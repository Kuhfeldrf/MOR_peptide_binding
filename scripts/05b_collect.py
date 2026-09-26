#!/usr/bin/env python3
"""Collect per-peptide MM/GBSA scores into one table.

Also reports how many of the scored peptides are scrambled decoys, because a
triage score computed only over presumed binders cannot be validated - that
was the failure of the co-folding confidence tier, and this table is what
replaces it.
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import sys


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, type=pathlib.Path)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    a = ap.parse_args()

    rows: dict[str, dict[str, float]] = {}
    for f in sorted(a.dir.glob("*/score.tsv")):
        with f.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rows.setdefault(row["peptide"], {})[row["metric"]] = \
                    float(row["value_kcal_per_mol"])

    if not rows:
        sys.exit(f"FATAL: no score.tsv found under {a.dir}")

    metrics = sorted({m for v in rows.values() for m in v})
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w") as fh:
        fh.write("peptide\tis_decoy\t" + "\t".join(metrics) + "\n")
        for pep in sorted(rows):
            vals = "\t".join(
                f"{rows[pep][m]:.3f}" if m in rows[pep] else "NA"
                for m in metrics)
            fh.write(f"{pep}\t{str(pep.startswith('scr__')).lower()}\t{vals}\n")

    n_dec = sum(1 for p in rows if p.startswith("scr__"))
    print(f"{len(rows)} peptides scored, {n_dec} of them scrambled decoys")
    if n_dec == 0:
        print("WARNING: no decoys present. The score cannot be validated "
              "from this table - discrimination needs negatives.")
    print(f"Wrote {a.out}")


if __name__ == "__main__":
    main()
