#!/usr/bin/env python3
"""Stage 7a - does the co-folding triage separate binders from decoys?

Two independent questions, with very different statistical power:

  DISCRIMINATION  real peptides vs matched scrambled controls (17 vs 20).
                  Reported as AUROC, which is exactly the probability that a
                  randomly chosen real peptide scores above a randomly chosen
                  decoy, plus a Mann-Whitney U test. This asks whether the
                  triage works at all.

  CALIBRATION     Stage 2 score vs measured GPI IC50 for the benchmark
                  agonists (n=7). Spearman, because the relationship need only
                  be monotonic. This asks whether the triage ranks correctly
                  among things that do bind.

The second is deliberately reported with its limits attached: with n=7 a
Spearman rho needs roughly |rho| >= 0.79 to reach p < 0.05, and six of the
seven peptides sit inside a 9-fold IC50 window, so most of the available
signal is one peptide separating from the rest. Saying so is part of the
result.

Nothing here decides anything. It reports numbers for a human to judge.
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

# Higher is better for all of these except the PAE columns, where lower is.
LOWER_IS_BETTER = {"pae_interface_mean", "pae_interface_min",
                   "peptide_rmsd_mean_A"}


def auroc(pos: np.ndarray, neg: np.ndarray) -> float:
    """P(a random positive outranks a random negative), ties counted as half.

    This is the Mann-Whitney U statistic normalised, so it is computed from U
    rather than by sorting - identical value, and no tie-handling to get wrong.
    """
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return float(u / (len(pos) * len(neg)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cofold-dir", required=True, type=pathlib.Path)
    ap.add_argument("--reference", required=True, type=pathlib.Path)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    scores = pd.read_csv(a.cofold_dir / "scores.tsv", sep="\t")
    agree = pd.read_csv(a.cofold_dir / "cross_seed_agreement.tsv", sep="\t")
    ref = pd.read_csv(a.reference, sep="\t")

    # scores.tsv has one row per seed AND per model. Collapse to one row per
    # peptide by taking the BEST value of each metric across all of them, which
    # is what a screen would act on - the best pose it found, not the average
    # of poses it would have discarded.
    metrics = ["aggregate_score", "iptm", "i_plddt_peptide",
               "pae_interface_mean", "pae_interface_min"]
    metrics = [m for m in metrics if m in scores.columns]
    agg = {m: ("min" if m in LOWER_IS_BETTER else "max") for m in metrics}
    per_pep = scores.groupby("peptide_id").agg(agg).reset_index()

    if "peptide_rmsd_mean_A" in agree.columns:
        per_pep = per_pep.merge(
            agree[["peptide_id", "peptide_rmsd_mean_A", "cross_seed_agreement"]],
            on="peptide_id", how="left")
        metrics.append("peptide_rmsd_mean_A")

    ipsae_path = a.cofold_dir / "ipsae.tsv"
    if ipsae_path.exists():
        ip = pd.read_csv(ipsae_path, sep="\t")
        keep = [c for c in ("ipSAE", "pDockQ", "pDockQ2", "LIS")
                if c in ip.columns]
        if keep:
            ipa = ip.groupby("peptide_id")[keep].max().reset_index()
            per_pep = per_pep.merge(ipa, on="peptide_id", how="left")
            metrics += keep

    per_pep["is_decoy"] = per_pep["peptide_id"].str.startswith("scr__")
    per_pep.to_csv(a.out / "triage_per_peptide.tsv", sep="\t", index=False)

    n_real = int((~per_pep["is_decoy"]).sum())
    n_dec = int(per_pep["is_decoy"].sum())
    print(f"peptides: {n_real} real, {n_dec} scrambled decoys\n")

    if n_dec == 0:
        print("NO DECOYS PRESENT - discrimination cannot be assessed.")
        print("Set library.include_scrambles and re-run Stage 2.")

    # ------------------------------------------------ discrimination
    rows = []
    for m in metrics:
        pos = per_pep.loc[~per_pep["is_decoy"], m].dropna().to_numpy()
        neg = per_pep.loc[per_pep["is_decoy"], m].dropna().to_numpy()
        if len(pos) < 3 or len(neg) < 3:
            continue
        # Orient so that "higher is better" holds, making AUROC directly
        # readable: 0.5 is chance, 1.0 is perfect, below 0.5 is backwards.
        sign = -1.0 if m in LOWER_IS_BETTER else 1.0
        au = auroc(sign * pos, sign * neg)
        p = stats.mannwhitneyu(pos, neg, alternative="two-sided").pvalue
        rows.append({"metric": m, "n_real": len(pos), "n_decoy": len(neg),
                     "real_median": float(np.median(pos)),
                     "decoy_median": float(np.median(neg)),
                     "auroc": au, "mannwhitney_p": float(p)})

    disc = (pd.DataFrame(rows).sort_values("auroc", ascending=False)
            if rows else pd.DataFrame())
    if not disc.empty:
        disc.to_csv(a.out / "discrimination.tsv", sep="\t", index=False)
        print("DISCRIMINATION - real vs scrambled")
        print(f"{'metric':<22} {'real':>9} {'decoy':>9} {'AUROC':>7} {'p':>9}")
        for _, r in disc.iterrows():
            print(f"{r['metric']:<22} {r['real_median']:>9.3f} "
                  f"{r['decoy_median']:>9.3f} {r['auroc']:>7.3f} "
                  f"{r['mannwhitney_p']:>9.2e}")
        print("\nAUROC 0.5 = chance. Below 0.5 means the metric ranks decoys "
              "ABOVE real peptides.")

    # ------------------------------------------------ calibration
    bench = ref[(ref.get("benchmark_include") == "YES")
                & ref["ic50_um"].notna()][["name", "ic50_um"]]
    cal = per_pep.merge(bench, left_on="peptide_id", right_on="name")
    print(f"\nCALIBRATION - Stage 2 score vs GPI IC50 (n={len(cal)})")
    if len(cal) < 4:
        print("  too few peptides with IC50 to correlate")
    else:
        crows = []
        for m in metrics:
            v = cal[[m, "ic50_um"]].dropna()
            if len(v) < 4:
                continue
            # Potency is inverse to IC50, so a GOOD metric anti-correlates
            # with IC50. Sign is flipped so positive rho always means
            # "tracks potency".
            sign = -1.0 if m in LOWER_IS_BETTER else 1.0
            x = sign * v[m].to_numpy()
            y = -np.log10(v["ic50_um"].to_numpy())
            rho, p = stats.spearmanr(x, y)
            # Jackknife: with n=7 a single peptide can carry the whole
            # correlation. Dropping each in turn shows whether it does. The
            # worst case is what the claim is actually worth.
            jack = [stats.spearmanr(np.delete(x, i), np.delete(y, i)).statistic
                    for i in range(len(x))]
            crows.append({"metric": m, "n": len(v), "spearman_rho": float(rho),
                          "p": float(p),
                          "jackknife_min_rho": float(np.min(jack)),
                          "jackknife_max_rho": float(np.max(jack))})
        cdf = pd.DataFrame(crows).sort_values("spearman_rho", ascending=False)
        cdf.to_csv(a.out / "calibration.tsv", sep="\t", index=False)
        print(f"{'metric':<22} {'n':>3} {'rho':>7} {'p':>8} {'jackknife rho':>16}")
        for _, r in cdf.iterrows():
            print(f"{r['metric']:<22} {r['n']:>3.0f} {r['spearman_rho']:>7.3f} "
                  f"{r['p']:>8.3f}   "
                  f"{r['jackknife_min_rho']:>6.2f}..{r['jackknife_max_rho']:>5.2f}")

        # Two corrections that the headline number does not survive unaided.
        k = len(cdf)
        best = cdf.iloc[0]
        print(f"\n  n={len(cal)}: |rho| >= ~0.79 needed for p<0.05. "
              f"IC50 range {cal['ic50_um'].min():.2g}-"
              f"{cal['ic50_um'].max():.2g} uM.")
        print(f"  MULTIPLE COMPARISONS: {k} metrics were tested. Bonferroni "
              f"on the best (p={best['p']:.3f}) gives p={min(1.0, best['p']*k):.2f} "
              f"- not significant.")
        print(f"  JACKKNIFE: dropping one peptide moves the best metric's rho "
              f"to as low as {best['jackknife_min_rho']:.2f}.")
        print("  Treat this as a promising signal to test on a larger set, "
              "not an established correlation.")

    print(f"\nWrote {a.out}/")


if __name__ == "__main__":
    main()
