#!/usr/bin/env python3
"""Is MM/GBSA measuring binding, or just ligand size?"""
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("results/05_mmgbsa/scores.tsv", sep="\t")
d["is_decoy"] = d.is_decoy.astype(str).str.lower() == "true"

# ligand heavy-atom count, straight from the parameterised ligand
n_at, n_res = [], []
for pep in d.peptide:
    f = pathlib.Path(f"results/03_membrane/{pep}/params/ligand_unique.pdb")
    lines = [l for l in f.read_text().splitlines()
             if l.startswith(("ATOM", "HETATM"))]
    n_at.append(len(lines))
    seq = pathlib.Path("results/01_library/peptides.tsv").read_text()
    n_res.append(np.nan)
d["n_ligand_atoms"] = n_at

for m in ["delta_total", "vdwaals", "eel"]:
    rho, p = stats.spearmanr(d[m], d.n_ligand_atoms)
    r = np.corrcoef(d[m], d.n_ligand_atoms)[0, 1]
    print(f"{m:<12} vs ligand atom count:  spearman rho={rho:+.3f} "
          f"(p={p:.2g})   pearson r={r:+.3f}")

print()
print("size-normalised delta_total (per ligand heavy atom):")
d["per_atom"] = d.delta_total / d.n_ligand_atoms
real, dec = d[~d.is_decoy], d[d.is_decoy]
for label, col in [("raw delta_total", "delta_total"),
                   ("per-atom", "per_atom")]:
    a = real[col].to_numpy()
    b = dec[col].to_numpy()
    u = stats.mannwhitneyu(-a, -b, alternative="two-sided")
    print(f"  {label:<16} real {np.median(a):>8.3f}  decoy {np.median(b):>8.3f}"
          f"   AUROC {u.statistic/(len(a)*len(b)):.3f}  p={u.pvalue:.3f}")

ref = pd.read_csv("docs/reference_peptides.tsv", sep="\t")
bench = ref[(ref.benchmark_include == "YES") & ref.ic50_um.notna()][["name", "ic50_um"]]
c = d.merge(bench, left_on="peptide", right_on="name")
print(f"\ncalibration vs IC50 (n={len(c)}):")
for col in ["delta_total", "per_atom"]:
    rho, p = stats.spearmanr(-c[col], -np.log10(c.ic50_um))
    print(f"  {col:<14} rho={rho:+.3f}  p={p:.3f}")
print()
print(c[["peptide", "ic50_um", "n_ligand_atoms", "delta_total", "per_atom"]]
      .sort_values("ic50_um").to_string(index=False))
