#!/usr/bin/env python3
"""Do the length controls separate message from length?"""
import numpy as np
import pandas as pd
from scipy import stats

s = pd.read_csv("results/02_cofold_lenctrl/scores.tsv", sep="\t")
ag = pd.read_csv("results/02_cofold_lenctrl/cross_seed_agreement.tsv", sep="\t")

per = s.groupby("peptide_id").agg(
    iptm=("iptm", "max"),
    i_plddt=("i_plddt_peptide", "max"),
    pae_min=("pae_interface_min", "min"),
    seq=("sequence", "first")).reset_index()
per = per.merge(ag[["peptide_id", "peptide_rmsd_mean_A", "cross_seed_agreement"]],
                on="peptide_id", how="left")
per["L"] = per.seq.str.len()
per["arm"] = np.where(per.peptide_id.str.startswith("pos__"), "positive", "negative")

print("POSITIVES (beta-endorphin N-terminal fragments, message YGGFM present)")
print(per[per.arm == "positive"][["peptide_id", "L", "iptm", "pae_min",
                                  "peptide_rmsd_mean_A", "cross_seed_agreement"]]
      .sort_values("L").to_string(index=False))

print("\nNEGATIVES by length (real milk fragments, no N-terminal Tyr)")
for L in (10, 15, 20):
    g = per[(per.arm == "negative") & (per.L == L)]
    p = per[(per.arm == "positive") & (per.L == L)]
    print(f"\n  {L} aa   negatives n={len(g)}")
    print(f"    iptm    negatives median {g.iptm.median():.4f}  "
          f"range {g.iptm.min():.4f}-{g.iptm.max():.4f}   "
          f"POSITIVE {p.iptm.iloc[0]:.4f}")
    rank = (g.iptm >= p.iptm.iloc[0]).sum()
    print(f"    negatives scoring >= the positive: {rank}/{len(g)}")
    print(f"    seed RMSD  negatives median {g.peptide_rmsd_mean_A.median():.2f} A"
          f"   POSITIVE {p.peptide_rmsd_mean_A.iloc[0]:.2f} A")

print("\nOVERALL discrimination, positives vs negatives (all lengths)")
a = per[per.arm == "positive"].iptm.to_numpy()
b = per[per.arm == "negative"].iptm.to_numpy()
u = stats.mannwhitneyu(a, b, alternative="two-sided")
print(f"  n=3 vs {len(b)}   AUROC {u.statistic/(len(a)*len(b)):.3f}   p={u.pvalue:.3f}")
print("  (n=3 positives - this is descriptive, not a powered test)")

print("\nDoes ipTM just track length?")
rho, p = stats.spearmanr(per.iptm, per.L)
print(f"  all controls:  rho={rho:+.3f}  p={p:.2g}")
g = per[per.arm == "negative"]
rho2, p2 = stats.spearmanr(g.iptm, g.L)
print(f"  negatives only: rho={rho2:+.3f}  p={p2:.2g}")

print("\nFull table")
print(per[["peptide_id", "L", "arm", "iptm", "pae_min",
           "peptide_rmsd_mean_A", "cross_seed_agreement"]]
      .sort_values(["L", "arm", "iptm"]).to_string(index=False))
