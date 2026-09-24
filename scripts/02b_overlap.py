#!/usr/bin/env python3
"""Stage 2b - training-overlap audit.

For every peptide, measure sequence similarity to peptide ligands present in
PDB mu-opioid receptor complexes, and emit a HIGH_OVERLAP / MODERATE / NOVEL
flag.

Why this stage exists: co-folding accuracy depends strongly on how much a
target resembles the model's training data. Reference peptides will predict
well partly through memorisation. Without this column, success on the reference
set would be mistaken for evidence of generalisation. Stage 7 stratifies by
this flag, and performance on NOVEL peptides is the honest measure.

CPU only, no GPU, seconds to run.

Usage:
    02b_overlap.py --peptides results/01_library/peptides.tsv \
                   --config config/config.yaml \
                   --outdir results/02b_overlap
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
from datetime import datetime, timezone

import yaml

# ---------------------------------------------------------------------------
# Known peptide ligands in deposited opioid-receptor structures.
#
# These are the sequences a co-folding model is most likely to have memorised.
# Each entry records the PDB entry it comes from and the deposition year, so
# the basis for the flag is auditable and can be extended as new structures
# appear. Curated 2026-09-24; DAMGO verified directly against data/raw/6ddf.cif.
# ---------------------------------------------------------------------------
PDB_OPIOID_PEPTIDES = [
    # (sequence, name, pdb, year, note)
    ("YGGFM",   "Met-enkephalin",   "8F7Q", 2023, "endogenous opioid peptide-receptor-Gi"),
    ("YGGFL",   "Leu-enkephalin",   "8F7R", 2023, "endogenous opioid peptide-receptor-Gi"),
    ("YGGFMRF", "Met-enkephalin-RF","8F7S", 2023, "endogenous opioid peptide-receptor-Gi"),
    ("YGGFLRRI","Dynorphin A(1-8)", "8F7W", 2023, "kappa/mu opioid peptide complexes"),
    ("YPWF",    "Endomorphin-1",    "-",    None, "widely studied mu agonist; in training corpora"),
    ("YPFF",    "Endomorphin-2",    "-",    None, "widely studied mu agonist; in training corpora"),
    ("YAGFMX",  "DAMGO",            "6DDF", 2018, "non-canonical: D-Ala2, N-MePhe4, Gly-ol5"),
]


def identity(a: str, b: str) -> float:
    """Best ungapped local identity, normalised by the shorter sequence.

    Deliberately simple: this is a memorisation screen, not an alignment study.
    A short peptide embedded in a longer one is the case that matters, and a
    sliding ungapped comparison captures it without introducing gap-penalty
    parameters that would need their own justification.
    """
    if not a or not b:
        return 0.0
    s, l = (a, b) if len(a) <= len(b) else (b, a)
    best = 0
    for off in range(len(l) - len(s) + 1):
        m = sum(1 for i, c in enumerate(s) if l[off + i] == c)
        best = max(best, m)
    return best / len(s)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--peptides", required=True, type=pathlib.Path)
    ap.add_argument("--config", required=True, type=pathlib.Path)
    ap.add_argument("--outdir", required=True, type=pathlib.Path)
    a = ap.parse_args()

    cfg = yaml.safe_load(a.config.read_text())["overlap"]
    hi = float(cfg["high_identity"])
    mod = float(cfg["moderate_identity"])
    a.outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    with a.peptides.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            seq = r["sequence"].strip().upper()
            if not seq:
                continue
            best_id, best_hit = 0.0, None
            for ref_seq, name, pdb, year, note in PDB_OPIOID_PEPTIDES:
                v = identity(seq, ref_seq)
                if v > best_id:
                    best_id, best_hit = v, (name, pdb, year)
            exact = any(seq == s for s, *_ in PDB_OPIOID_PEPTIDES)
            if exact or best_id >= hi:
                flag = "HIGH_OVERLAP"
            elif best_id >= mod:
                flag = "MODERATE"
            else:
                flag = "NOVEL"
            rows.append({
                "peptide_id": r["peptide_id"],
                "sequence": seq,
                "length": len(seq),
                "is_control": r.get("is_control", ""),
                "benchmark_include": r.get("benchmark_include", "NO"),
                "overlap_flag": flag,
                "best_identity": f"{best_id:.3f}",
                "exact_match_in_pdb": "yes" if exact else "no",
                "closest_known": best_hit[0] if best_hit else "",
                "closest_pdb": best_hit[1] if best_hit else "",
                "closest_year": best_hit[2] if best_hit and best_hit[2] else "",
            })

    out = a.outdir / "overlap.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["peptide_id"]):
            w.writerow(r)

    counts: dict[str, int] = {}
    bench: dict[str, int] = {}
    for r in rows:
        counts[r["overlap_flag"]] = counts.get(r["overlap_flag"], 0) + 1
        if r["benchmark_include"] == "YES":
            bench[r["overlap_flag"]] = bench.get(r["overlap_flag"], 0) + 1

    stats = {
        "n_peptides": len(rows),
        "by_flag": counts,
        "benchmark_by_flag": bench,
        "thresholds": {"high_identity": hi, "moderate_identity": mod},
        "reference_structures": len(PDB_OPIOID_PEPTIDES),
    }
    (a.outdir / "overlap_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n")

    print(f"Peptides audited: {len(rows)}")
    for k in ("HIGH_OVERLAP", "MODERATE", "NOVEL"):
        print(f"  {k:<13} {counts.get(k,0):>3}   "
              f"(benchmark set: {bench.get(k,0)})")
    print(f"\nWrote {out}")
    if bench.get("NOVEL", 0) == 0:
        print("\nWARNING: no benchmark peptide is flagged NOVEL. Stage 7 "
              "performance on this set therefore cannot distinguish "
              "generalisation from memorisation.")


if __name__ == "__main__":
    main()
