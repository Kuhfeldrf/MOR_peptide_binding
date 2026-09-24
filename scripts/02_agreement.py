#!/usr/bin/env python3
"""Stage 2 post-processing - cross-seed agreement.

Do the seeds converge on the same binding site and peptide orientation?

This is a cheap and genuinely informative diagnostic: disagreement across seeds
flags a prediction that should not be trusted, independently of how confident
the model claims to be. A high ipTM with poor cross-seed agreement is a
stronger warning than either number alone.

Method: superpose each seed's top-ranked model onto a reference seed using
RECEPTOR CA atoms only, then measure the peptide's displacement in that frame.
Superposing on the receptor rather than the whole complex is the point - it
asks "does the peptide land in the same place on the same receptor", not "are
the two structures similar overall", which a 281-residue receptor would
dominate.

Runs in `mor-pilot` (needs gemmi + numpy>=2), not `mor-chai`.

Usage:
    02_agreement.py --cofold results/02_cofold --outdir results/02_cofold
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import pathlib

import gemmi
import numpy as np


def load_top_model(seeddir: pathlib.Path) -> pathlib.Path | None:
    """Chai writes pred.model_idx_N.cif ranked by aggregate score; idx 0 is top."""
    c = sorted(seeddir.glob("pred.model_idx_*.cif"))
    return c[0] if c else None


def chain_ca(st: gemmi.Structure, chain_index: int) -> np.ndarray:
    ch = st[0][chain_index]
    pts = []
    for res in ch:
        at = res.find_atom("CA", "*")
        if at:
            pts.append([at.pos.x, at.pos.y, at.pos.z])
    return np.asarray(pts, dtype=float)


def kabsch(P: np.ndarray, Q: np.ndarray):
    """Rotation+translation mapping P onto Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cofold", required=True, type=pathlib.Path)
    ap.add_argument("--outdir", required=True, type=pathlib.Path)
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for pepdir in sorted(p for p in a.cofold.iterdir()
                         if p.is_dir() and not p.name.startswith("_")):
        seeds = sorted(d for d in pepdir.glob("seed*") if d.is_dir())
        models = {}
        for sd in seeds:
            cif = load_top_model(sd)
            if not cif:
                continue
            st = gemmi.read_structure(str(cif))
            st.setup_entities()
            if len(st[0]) < 2:
                continue
            models[sd.name] = (chain_ca(st, 0), chain_ca(st, 1))

        if len(models) < 2:
            print(f"{pepdir.name}: fewer than 2 usable seeds, skipping")
            continue

        names = sorted(models)
        pair_rmsd, pair_centroid = [], []
        for x, y in itertools.combinations(names, 2):
            rx, px = models[x]
            ry, py = models[y]
            n = min(len(rx), len(ry))
            if n < 3 or len(px) != len(py) or len(px) == 0:
                continue
            R, t = kabsch(rx[:n], ry[:n])
            px_al = (R @ px.T).T + t
            pair_rmsd.append(float(np.sqrt(np.mean(np.sum((px_al - py) ** 2, axis=1)))))
            pair_centroid.append(float(np.linalg.norm(px_al.mean(0) - py.mean(0))))

        if not pair_rmsd:
            continue
        mean_r, max_r = float(np.mean(pair_rmsd)), float(np.max(pair_rmsd))
        # Thresholds are reporting bands, not pass/fail: 2 A is roughly the
        # scale of a well-converged pose, 5 A means a different sub-site.
        verdict = ("CONVERGED" if max_r < 2.0
                   else "PARTIAL" if max_r < 5.0
                   else "DIVERGENT")
        rows.append({
            "peptide_id": pepdir.name,
            "n_seeds": len(names),
            "n_pairs": len(pair_rmsd),
            "peptide_rmsd_mean_A": round(mean_r, 3),
            "peptide_rmsd_max_A": round(max_r, 3),
            "centroid_shift_mean_A": round(float(np.mean(pair_centroid)), 3),
            "centroid_shift_max_A": round(float(np.max(pair_centroid)), 3),
            "cross_seed_agreement": verdict,
        })
        print(f"{pepdir.name}: {len(names)} seeds  "
              f"peptide RMSD mean={mean_r:.2f} max={max_r:.2f} A  -> {verdict}")

    if not rows:
        raise SystemExit("FATAL: no peptide had two usable seeds; "
                         "cross-seed agreement could not be computed.")

    out = a.outdir / "cross_seed_agreement.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    tally: dict[str, int] = {}
    for r in rows:
        tally[r["cross_seed_agreement"]] = tally.get(r["cross_seed_agreement"], 0) + 1
    (a.outdir / "agreement_stats.json").write_text(
        json.dumps({"n_peptides": len(rows), "by_verdict": tally,
                    "superposition": "receptor CA only",
                    "bands_A": {"CONVERGED": "<2", "PARTIAL": "2-5", "DIVERGENT": ">5"}},
                   indent=2, sort_keys=True) + "\n")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
