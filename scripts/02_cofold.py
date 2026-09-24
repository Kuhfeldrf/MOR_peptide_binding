#!/usr/bin/env python3
"""Stage 2 - Chai-1 co-folding, single-sequence mode.

Runs one inference per (peptide, seed) and records the confidence metrics.
Runs in the `mor-chai` environment (numpy<2); cross-seed geometric agreement is
computed separately by 02_agreement.py in `mor-pilot`, because the two cannot
share an environment (see docs/decisions.md D5).

INTERPRETATION NOTE, to be carried into the README: these are *learned
confidence* scores, not binding affinity. Ranking is by the model's estimate of
its own accuracy. A confidently ranked top pose can still be confidently wrong.
This is the stated reason the pipeline continues into physics.

Usage:
    02_cofold.py --receptor results/00_receptor/mOR_clean.fasta \
                 --peptides results/01_library/peptides.fasta \
                 --outdir results/02_cofold --seeds 5 [--only ID,ID,...]
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import time

import numpy as np
import torch
from chai_lab.chai1 import run_inference


def read_fasta(p: pathlib.Path) -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    for line in p.read_text().splitlines():
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        elif line.strip():
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def interface_plddt(plddt, n_receptor: int) -> tuple[float, float]:
    """Mean pLDDT over the peptide chain, and over the whole complex.

    The peptide is the second chain, so its tokens follow the receptor's. This
    is the i-pLDDT proxy: confidence localised to the binding partner rather
    than averaged over a large, well-folded receptor that would swamp it.
    """
    a = np.asarray(plddt).squeeze()
    if a.ndim == 0:
        return float(a), float(a)
    flat = a.reshape(-1) if a.ndim == 1 else a.reshape(a.shape[0], -1)[0]
    if flat.size <= n_receptor:
        return float("nan"), float(np.mean(flat))
    return float(np.mean(flat[n_receptor:])), float(np.mean(flat))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True, type=pathlib.Path)
    ap.add_argument("--peptides", required=True, type=pathlib.Path)
    ap.add_argument("--outdir", required=True, type=pathlib.Path)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--only", default="", help="comma-separated peptide ids")
    a = ap.parse_args()

    receptor = read_fasta(a.receptor)[0][1]
    peptides = read_fasta(a.peptides)
    if a.only:
        want = {s.strip() for s in a.only.split(",") if s.strip()}
        peptides = [(n, s) for n, s in peptides if n in want]
        missing = want - {n for n, _ in peptides}
        if missing:
            raise SystemExit(f"FATAL: peptide ids not found: {sorted(missing)}")

    a.outdir.mkdir(parents=True, exist_ok=True)
    print(f"receptor: {len(receptor)} aa")
    print(f"peptides: {len(peptides)}  seeds: {a.seeds}")
    print(f"GPU: {torch.cuda.get_device_name(0)}  bf16={torch.cuda.is_bf16_supported()}")

    rows: list[dict] = []
    t_start = time.time()

    for pid, pseq in peptides:
        for seed in range(a.seeds):
            tag = f"{pid}_seed{seed}"
            rundir = a.outdir / pid / f"seed{seed}"
            # Chai-1 asserts its output_dir is empty, so inputs live elsewhere.
            indir = a.outdir / pid / "inputs"
            indir.mkdir(parents=True, exist_ok=True)
            fa = indir / f"seed{seed}.fasta"
            fa.write_text(f">protein|name=receptor\n{receptor}\n"
                          f">protein|name=peptide\n{pseq}\n")

            # Restartable: a seed that already produced scores is not redone.
            done = sorted(rundir.glob("scores.model_idx_*.npz")) if rundir.exists() else []
            if done:
                print(f"  {tag}: already complete ({len(done)} models), skipping")
                dt = float("nan")
            else:
                if rundir.exists():
                    for f in rundir.iterdir():
                        f.unlink()
                t0 = time.time()
                cand = run_inference(
                    fasta_file=fa,
                    output_dir=rundir,
                    use_msa_server=False,   # single-sequence mode, deliberate
                    use_esm_embeddings=True,
                    num_trunk_recycles=3,
                    num_diffn_timesteps=200,
                    seed=seed,
                    device="cuda:0",
                )
                dt = time.time() - t0
                np.savez(rundir / "_cand_meta.npz",
                         plddt=np.asarray(cand.plddt),
                         pae=np.asarray(cand.pae))

            meta = np.load(rundir / "_cand_meta.npz")
            ipl, cpl = interface_plddt(meta["plddt"], len(receptor))
            pae = meta["pae"]
            # Interface PAE: receptor-vs-peptide block, both directions.
            try:
                p2 = pae.reshape(-1, pae.shape[-2], pae.shape[-1])[0]
                nr = len(receptor)
                inter = np.concatenate([p2[:nr, nr:].reshape(-1),
                                        p2[nr:, :nr].reshape(-1)])
                pae_inter = float(np.mean(inter))
                pae_min = float(np.min(inter))
            except Exception:
                pae_inter = pae_min = float("nan")

            for idx, npz in enumerate(sorted(rundir.glob("scores.model_idx_*.npz"))):
                d = np.load(npz, allow_pickle=True)
                g = lambda k: float(np.ravel(d[k])[0]) if k in d.files else float("nan")
                pc = np.ravel(d["per_chain_ptm"]) if "per_chain_ptm" in d.files else [np.nan, np.nan]
                rows.append({
                    "peptide_id": pid,
                    "sequence": pseq,
                    "seed": seed,
                    "model_idx": idx,
                    "aggregate_score": round(g("aggregate_score"), 5),
                    "iptm": round(g("iptm"), 5),
                    "ptm": round(g("ptm"), 5),
                    "receptor_ptm": round(float(pc[0]), 5),
                    "peptide_ptm": round(float(pc[1]), 5) if len(pc) > 1 else "",
                    "i_plddt_peptide": round(ipl, 3),
                    "plddt_complex": round(cpl, 3),
                    "pae_interface_mean": round(pae_inter, 3),
                    "pae_interface_min": round(pae_min, 3),
                    # ipSAE is NOT produced by Chai-1. See docs/decisions.md D21.
                    "ipsae": "NOT_COMPUTED",
                    "has_inter_chain_clashes": bool(np.ravel(d["has_inter_chain_clashes"])[0])
                        if "has_inter_chain_clashes" in d.files else "",
                    "runtime_s": round(dt, 1),
                })
            print(f"  {tag}: {dt:.0f}s  iptm={rows[-1]['iptm']}  "
                  f"agg={rows[-1]['aggregate_score']}  i-pLDDT={rows[-1]['i_plddt_peptide']}")

    out = a.outdir / "scores.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    total = time.time() - t_start
    n_runs = len(peptides) * a.seeds
    (a.outdir / "cofold_stats.json").write_text(json.dumps({
        "n_peptides": len(peptides),
        "seeds_per_peptide": a.seeds,
        "inference_runs": n_runs,
        "structures_per_run": 5,
        "total_structures": len(rows),
        "total_runtime_s": round(total, 1),
        "mean_runtime_per_run_s": round(total / max(n_runs, 1), 1),
        "msa": False,
        "esm_embeddings": True,
        "num_trunk_recycles": 3,
        "num_diffn_timesteps": 200,
    }, indent=2, sort_keys=True) + "\n")

    print(f"\n{len(rows)} structures from {n_runs} runs in {total/60:.1f} min")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
