#!/usr/bin/env python3
"""Stage 2, arm B - Protenix co-folding, single-sequence mode.

Drop-in counterpart to 02_cofold.py (arm A, Chai-1). Same CLI, same one
(peptide, seed) per invocation, and the SAME scores.tsv schema, so
02_agreement.py and 07b_length_controls_analysis.py run unmodified against it.
That is what keeps arm A vs arm B a replication rather than a re-analysis.

Runs in the `mor-protenix` environment (numpy 2.4.1, torch 2.7.1), which cannot
be shared with `mor-chai` (numpy<2) - see docs/decisions.md D5.

MSA-free: no pairedMsaPath / unpairedMsaPath / templatesPath in the input JSON,
plus --use_msa false. Matches arm A's deliberate single-sequence mode.

INTERPRETATION NOTE, identical to arm A: these are *learned confidence* scores,
not binding affinity. Ranking is by the model's estimate of its own accuracy.
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import shutil
import subprocess
import time

# Protenix's summary_confidence JSON does NOT contain a PAE matrix, so the two
# interface-PAE columns arm A fills cannot be computed here. They are written
# empty rather than substituted with a different quantity. Protenix's own
# interface error metric, chain_pair_gpde, is recorded in its own clearly named
# column instead - never silently mapped onto the PAE column names.
PAE_UNAVAILABLE = ""


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


def lookup_sequence(library: pathlib.Path, peptide: str) -> str:
    with library.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("name") == peptide or row.get("peptide_id") == peptide:
                return row["sequence"].strip().upper()
    raise SystemExit(f"FATAL: peptide {peptide!r} not in {library}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True, type=pathlib.Path)
    ap.add_argument("--outdir", required=True, type=pathlib.Path)
    ap.add_argument("--peptide", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--library", required=True, type=pathlib.Path)
    ap.add_argument("--recycles", type=int, default=3)
    ap.add_argument("--timesteps", type=int, default=200)
    ap.add_argument("--samples", type=int, default=5)
    ap.add_argument("--model", default="protenix_base_default_v1.0.0",
                    help="2021-09-30 training cutoff, which holds 8F7Q "
                         "(released 2022-12-14) OUT of training. Changing this "
                         "to the 20250630 checkpoint would put it IN.")
    # Receptor MSA. Generated ONCE by jobs/ptx_msa.sbatch and shared by every
    # task, because the receptor sequence is identical across all 21x5 runs.
    #
    # WHY THIS IS NOT OPTIONAL-BY-DEFAULT: arm B v1 ran MSA-free and predicted
    # the receptor 16.2 A from 8F7Q chain R - the 7TM bundle never folded, so
    # there was no pocket and the whole arm was void. Chai-1 (arm A) got away
    # with single-sequence only because it supplies use_esm_embeddings=True as
    # an MSA surrogate; Protenix has no equivalent.
    #
    # The PEPTIDE deliberately gets no MSA: a 10-20mer fragment carries no
    # usable co-evolutionary signal, and it keeps this arm closer to arm A,
    # where neither chain had one.
    ap.add_argument("--paired-msa", type=pathlib.Path, default=None,
                    help="receptor pairing.a3m (absolute path recommended)")
    ap.add_argument("--unpaired-msa", type=pathlib.Path, default=None,
                    help="receptor non_pairing.a3m (absolute path recommended)")
    # Arm A accepts an m8 template file. Protenix templates need .a3m or .hhr,
    # so an m8 cannot be honoured. Reject it loudly rather than ignore it.
    ap.add_argument("--template-m8", type=pathlib.Path, default=None)
    a = ap.parse_args()

    if a.template_m8 is not None:
        raise SystemExit(
            "FATAL: --template-m8 is not supported by Protenix. It accepts "
            "templates as .a3m (hmmsearch) or .hhr via templatesPath, not m8. "
            "Refusing to silently drop a template that was asked for.")

    receptor = read_fasta(a.receptor)[0][1]
    seq = lookup_sequence(a.library, a.peptide)

    a.outdir.mkdir(parents=True, exist_ok=True)
    tag = f"{a.peptide}_seed{a.seed}"
    print(f"receptor: {len(receptor)} aa")
    print(f"peptide: {a.peptide} = {seq} ({len(seq)} aa)  seed: {a.seed}")

    pred = a.outdir / "pred"
    scores_out = a.outdir / "scores.tsv"

    # Restartable, on the same rule arm A uses: a run counts as complete only if
    # the parsed output AND the per-sample confidence JSONs are present.
    done = sorted(pred.rglob("*_summary_confidence_sample_*.json"))
    if scores_out.exists() and len(done) >= a.samples:
        print(f"  {tag}: already complete ({len(done)} samples), skipping")
        return

    if pred.exists():
        shutil.rmtree(pred)

    receptor_entity: dict = {"sequence": receptor, "count": 1}
    use_msa = False
    if a.paired_msa or a.unpaired_msa:
        for flag, path in (("--paired-msa", a.paired_msa),
                           ("--unpaired-msa", a.unpaired_msa)):
            if path is not None and not path.is_file():
                raise SystemExit(f"FATAL: {flag} {path} does not exist. "
                                 "Refusing to fall back to MSA-free, which is "
                                 "the condition that voided arm B v1.")
        if a.paired_msa:
            receptor_entity["pairedMsaPath"] = str(a.paired_msa.resolve())
        if a.unpaired_msa:
            receptor_entity["unpairedMsaPath"] = str(a.unpaired_msa.resolve())
        use_msa = True
        print(f"  receptor MSA: paired={a.paired_msa} unpaired={a.unpaired_msa}")
    else:
        print("  WARNING: no receptor MSA supplied; this reproduces the VOID "
              "arm B v1 condition and the receptor will not fold correctly.")

    inp = a.outdir / "input.json"
    inp.write_text(json.dumps([{
        "name": tag,
        "sequences": [
            {"proteinChain": receptor_entity},
            {"proteinChain": {"sequence": seq, "count": 1}},
        ],
    }], indent=2) + "\n")

    cmd = [
        "protenix", "pred",
        "-i", str(inp), "-o", str(pred),
        "-s", str(a.seed),
        "-n", a.model,
        "-e", str(a.samples),
        "--use_msa", "true" if use_msa else "false",
        "--use_default_params", "false",
        "--cycle", str(a.recycles),
        "--step", str(a.timesteps),
    ]
    print("  " + " ".join(cmd))
    t0 = time.time()
    rc = subprocess.call(cmd)
    dt = time.time() - t0
    if rc != 0:
        raise SystemExit(f"FATAL: protenix pred exited {rc} for {tag}")

    jsons = sorted(pred.rglob("*_summary_confidence_sample_*.json"),
                   key=lambda p: int(p.stem.rsplit("_", 1)[-1]))
    if not jsons:
        raise SystemExit(f"FATAL: protenix produced no confidence JSON for {tag}")

    rows: list[dict] = []
    for j in jsons:
        d = json.loads(j.read_text())
        idx = int(j.stem.rsplit("_", 1)[-1])
        chain_ptm = d.get("chain_ptm") or [float("nan")] * 2
        chain_plddt = d.get("chain_plddt") or [float("nan")] * 2
        # Interface gpde: receptor-vs-peptide off-diagonal element.
        try:
            gpde_if = float(d["chain_pair_gpde"][0][1])
        except (KeyError, IndexError, TypeError):
            gpde_if = float("nan")
        rows.append({
            "peptide_id": a.peptide,
            "sequence": seq,
            "seed": a.seed,
            "model_idx": idx,
            # Protenix's own ranking metric, the analogue of Chai-1's
            # aggregate_score. Named to match so downstream code is unchanged.
            "aggregate_score": round(float(d.get("ranking_score", float("nan"))), 5),
            "iptm": round(float(d.get("iptm", float("nan"))), 5),
            "ptm": round(float(d.get("ptm", float("nan"))), 5),
            "receptor_ptm": round(float(chain_ptm[0]), 5),
            "peptide_ptm": round(float(chain_ptm[1]), 5) if len(chain_ptm) > 1 else "",
            # Protenix reports complex plddt on a 0-100 scale and per-chain
            # plddt on 0-1. Arm A's columns are 0-1, so the complex value is
            # divided by 100 to keep the two arms on ONE scale.
            "i_plddt_peptide": round(float(chain_plddt[1]), 3) if len(chain_plddt) > 1 else "",
            "plddt_complex": round(float(d.get("plddt", float("nan"))) / 100.0, 3),
            "pae_interface_mean": PAE_UNAVAILABLE,
            "pae_interface_min": PAE_UNAVAILABLE,
            "ipsae": "NOT_COMPUTED",
            "has_inter_chain_clashes": bool(d.get("has_clash", False)),
            "runtime_s": round(dt, 1),
            # Protenix-only extras, appended so nothing upstream is displaced.
            "gpde_interface": round(gpde_if, 4),
            "disorder": round(float(d.get("disorder", float("nan"))), 4),
        })

    with scores_out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    (a.outdir / "cofold_stats.json").write_text(json.dumps({
        "engine": "protenix",
        "model": a.model,
        "n_peptides": 1,
        "seed": a.seed,
        "structures_per_run": a.samples,
        "total_structures": len(rows),
        "total_runtime_s": round(dt, 1),
        "msa": use_msa,
        "paired_msa": str(a.paired_msa) if a.paired_msa else None,
        "unpaired_msa": str(a.unpaired_msa) if a.unpaired_msa else None,
        "template": False,
        "num_trunk_recycles": a.recycles,
        "num_diffn_timesteps": a.timesteps,
    }, indent=2, sort_keys=True) + "\n")

    print(f"  {tag}: {dt:.0f}s  iptm={rows[0]['iptm']}  "
          f"rank={rows[0]['aggregate_score']}  gpde_if={rows[0]['gpde_interface']}")
    print(f"Wrote {scores_out} ({len(rows)} structures)")


if __name__ == "__main__":
    main()
