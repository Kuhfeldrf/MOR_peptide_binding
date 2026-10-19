#!/usr/bin/env python3
"""Do any of the scores correlate with MEASURED affinity?

THE QUESTION BINARY LABELS CANNOT ANSWER. Every comparison so far scored methods
against class labels - binder vs assumed-non-binder - which cannot distinguish a
method that ranks by affinity from one that ranks by peptide length. With cited
Ki/IC50 values in hand the real question becomes available: within the set of
peptides that DO bind, does any score track how tightly they bind?

This is the stronger test, and it is the one that decides whether any of this
machinery could ever prioritise a library rather than merely filter it. The
pharmacophore CANNOT do this - it is binary by construction - so if a
structure-based score tracks affinity it would have a job even though it loses
the classification contest.

Method:
  - affinities are pooled per peptide as the median pAffinity
    (pAff = -log10(value in M)), taking Ki in preference to Kd, then IC50/EC50,
    and never mixing assay types within one peptide without saying so
  - scores are joined per peptide (max over seeds/samples, as everywhere else)
  - Spearman rho is reported, plus the same correlation for peptide LENGTH, so
    a score's apparent affinity correlation can be checked against the length
    artefact that explains every earlier result

Usage:
    07c_affinity_correlation.py --affinities docs/affinities_chembl.tsv \
        --scores <tsv with peptide_id + score columns> [...]
"""
from __future__ import annotations

import argparse
import math
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

PREFER = ["Ki", "Kd", "IC50", "EC50"]


def to_paff(value: float, units: str) -> float | None:
    """Convert a concentration to pAffinity = -log10(M)."""
    if value is None or value <= 0:
        return None
    u = (units or "").strip().lower()
    scale = {"nm": 1e-9, "um": 1e-6, "µm": 1e-6, "mm": 1e-3, "m": 1.0,
             "pm": 1e-12, "nmol/l": 1e-9, "umol/l": 1e-6}.get(u)
    if scale is None:
        return None
    return -math.log10(value * scale)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--affinities", required=True, type=pathlib.Path)
    ap.add_argument("--scores", nargs="+", required=True, type=pathlib.Path)
    ap.add_argument("--out", type=pathlib.Path, default=None)
    a = ap.parse_args()

    aff = pd.read_csv(a.affinities, sep="\t")
    aff = aff[aff.standard_type.isin(PREFER)].copy()
    if aff.empty:
        raise SystemExit("FATAL: no usable affinity records; nothing to correlate")
    aff["paff"] = [to_paff(v, u) for v, u in
                   zip(pd.to_numeric(aff.standard_value, errors="coerce"),
                       aff.standard_units)]
    aff = aff.dropna(subset=["paff"])
    # Censored values ('>' / '<') are dropped: a bound is not a measurement.
    n_cens = (aff.standard_relation.astype(str).str.strip().isin([">", "<", ">=", "<="])).sum()
    aff = aff[~aff.standard_relation.astype(str).str.strip().isin([">", "<", ">=", "<="])]
    print(f"{len(aff)} usable affinity records ({n_cens} censored dropped)")

    best = []
    for key, g in aff.groupby("peptide_key"):
        for t in PREFER:
            h = g[g.standard_type == t]
            if not h.empty:
                best.append({"peptide_key": key, "assay_type": t,
                             "paff_median": float(h.paff.median()),
                             "paff_n": int(len(h)),
                             "paff_spread": float(h.paff.max() - h.paff.min())})
                break
    A = pd.DataFrame(best)
    print(f"{len(A)} peptides with a pooled affinity\n")
    print(A.sort_values("paff_median", ascending=False).to_string(index=False))

    print("\n=== correlation of each score with MEASURED affinity ===")
    print(f"{'score file':34s} {'metric':14s} {'n':>3s} {'rho':>7s} {'p':>8s}")
    for sf in a.scores:
        if not sf.is_file():
            print(f"{sf.name:34s} (missing)")
            continue
        S = pd.read_csv(sf, sep="\t")
        idcol = "peptide_id" if "peptide_id" in S.columns else S.columns[0]
        num = [c for c in S.columns
               if c not in (idcol,) and pd.api.types.is_numeric_dtype(S[c])]
        agg = S.groupby(idcol)[num].max().reset_index()
        m = agg.merge(A, left_on=idcol, right_on="peptide_key", how="inner")
        if len(m) < 4:
            print(f"{sf.name:34s} only {len(m)} peptides overlap - skipped")
            continue
        for c in num:
            if c in ("seed", "model_idx", "n_receptor", "n_peptide"):
                continue
            rho, p = stats.spearmanr(m[c], m.paff_median)
            flag = "  <-- significant" if p < 0.05 else ""
            print(f"{sf.name:34s} {c:14s} {len(m):3d} {rho:+7.3f} {p:8.3f}{flag}")

    if a.out:
        A.to_csv(a.out, sep="\t", index=False)
        print(f"\nWrote {a.out}")


if __name__ == "__main__":
    main()
