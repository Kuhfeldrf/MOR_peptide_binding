#!/usr/bin/env python3
import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("results/05_mmgbsa/scores.tsv", sep="\t")
d["is_decoy"] = d.is_decoy.astype(str).str.lower() == "true"
real, dec = d[~d.is_decoy], d[d.is_decoy]
print(f"real n={len(real)}   decoys n={len(dec)}\n")

print("DISCRIMINATION  (more negative = better binding; oriented so AUROC>0.5 = works)")
print(f"{'metric':<14}{'real med':>10}{'decoy med':>11}{'AUROC':>8}{'p':>9}")
for m in ["delta_total", "vdwaals", "eel", "egb", "esurf"]:
    a = real[m].dropna().to_numpy()
    b = dec[m].dropna().to_numpy()
    u = stats.mannwhitneyu(-a, -b, alternative="two-sided")
    print(f"{m:<14}{np.median(a):>10.2f}{np.median(b):>11.2f}"
          f"{u.statistic/(len(a)*len(b)):>8.3f}{u.pvalue:>9.3f}")

ref = pd.read_csv("docs/reference_peptides.tsv", sep="\t")
bench = ref[(ref.benchmark_include == "YES") & ref.ic50_um.notna()][["name", "ic50_um"]]
c = d.merge(bench, left_on="peptide", right_on="name")
print(f"\nCALIBRATION vs GPI IC50  (n={len(c)})")
if len(c) >= 4:
    for m in ["delta_total", "vdwaals", "eel"]:
        rho, p = stats.spearmanr(-c[m], -np.log10(c.ic50_um))
        print(f"  {m:<12} rho={rho:+.3f}  p={p:.3f}")
    print()
    print(c[["peptide", "ic50_um", "delta_total"]]
          .sort_values("ic50_um").to_string(index=False))
