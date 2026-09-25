#!/usr/bin/env python3
"""Stage 2e - collect co-folding scores across every peptide and seed.

One row per model, so the per-seed spread stays visible. Collapsing to a single
number per peptide would hide the variation that Stage 2 showed is the only
informative signal: the confidence scores themselves did not discriminate,
while the spread across seeds did.
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import statistics as st
import sys


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cofold", required=True, type=pathlib.Path)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    a = ap.parse_args()

    rows, fields = [], None
    for pepdir in sorted(d for d in a.cofold.iterdir()
                         if d.is_dir() and not d.name.startswith("_")):
        for seeddir in sorted(d for d in pepdir.glob("seed*") if d.is_dir()):
            tsv = seeddir / "scores.tsv"
            if not tsv.exists():
                continue
            with tsv.open() as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    row.setdefault("peptide_id", pepdir.name)
                    row["peptide_id"] = pepdir.name
                    row["seed"] = seeddir.name.replace("seed", "")
                    rows.append(row)
                    if fields is None:
                        fields = list(row.keys())

    if not rows:
        sys.exit(f"FATAL: no per-seed scores found under {a.cofold}")

    order = ["peptide_id", "seed"] + [f for f in fields
                                      if f not in ("peptide_id", "seed")]
    with a.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=order, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["peptide_id"], x["seed"])):
            w.writerow(r)

    peptides = sorted({r["peptide_id"] for r in rows})
    print(f"{len(rows)} models across {len(peptides)} peptides")
    print(f"\n{'peptide':<32}{'iptm range':>18}{'aggregate':>12}")
    for p in peptides:
        vals = [float(r["iptm"]) for r in rows
                if r["peptide_id"] == p and r.get("iptm") not in (None, "", "nan")]
        agg = [float(r["aggregate_score"]) for r in rows
               if r["peptide_id"] == p
               and r.get("aggregate_score") not in (None, "", "nan")]
        if vals:
            print(f"  {p:<30}{min(vals):.3f}-{max(vals):.3f}"
                  f"{st.mean(agg) if agg else float('nan'):>12.3f}")
    print(f"\nWrote {a.out}")


if __name__ == "__main__":
    main()
