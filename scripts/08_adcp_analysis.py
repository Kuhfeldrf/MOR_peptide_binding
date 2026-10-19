#!/usr/bin/env python3
"""Does ADCP docking discriminate where co-folding did not?

Reads docs/adcp_energies.tsv (61 peptides x 30 replicas, fixed search effort)
and runs the SAME controls the co-folding arm failed, so the two are directly
comparable:

  1. length_positive vs length_negative at matched length (10/15/20 aa)
  2. real_peptide vs decoy_scramble
  3. does the score just track peptide length?
  4. the free baseline: "does the sequence start with Tyr?"
  5. replica convergence at fixed effort

Score convention: ADCP "best target energy", more negative = better. Per
peptide the statistic is the MINIMUM over 30 replicas (the best pose found),
with the replica spread kept as a convergence check.
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent.parent
df = pd.read_csv(HERE / "docs" / "adcp_energies.tsv", sep="\t")

per = df.groupby(["name", "class", "length", "sequence"]).agg(
    best=("best_target_energy", "min"),
    median_rep=("best_target_energy", "median"),
    worst=("best_target_energy", "max"),
    n_rep=("best_target_energy", "size"),
).reset_index()
per["spread"] = per.worst - per.best
per["starts_tyr"] = per.sequence.str.startswith("Y")

print(f"{len(df)} replica results over {len(per)} peptides "
      f"({per.n_rep.min()}-{per.n_rep.max()} replicas each)\n")
print("classes:")
print(per["class"].value_counts().to_string(), "\n")


def auroc(pos: np.ndarray, neg: np.ndarray) -> tuple[float, float]:
    """AUROC for 'more negative is better', plus Mann-Whitney p."""
    if len(pos) == 0 or len(neg) == 0:
        return float("nan"), float("nan")
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    # pos should be LOWER (better). AUROC in the "positive ranks better" sense:
    return 1.0 - u.statistic / (len(pos) * len(neg)), u.pvalue


print("=" * 72)
print("1. LENGTH CONTROLS - the test co-folding failed")
print("=" * 72)
lp = per[per["class"] == "length_positive"]
ln = per[per["class"] == "length_negative"]
for L in sorted(lp.length.unique()):
    p = lp[lp.length == L]
    g = ln[ln.length == L]
    if p.empty or g.empty:
        continue
    pv = p.best.iloc[0]
    better = (g.best <= pv).sum()
    print(f"\n  {L:2d} aa   negatives n={len(g)}")
    print(f"    ADCP best   negatives median {g.best.median():+.3f}  "
          f"range {g.best.min():+.3f}..{g.best.max():+.3f}")
    print(f"                POSITIVE {pv:+.3f}")
    print(f"    negatives scoring BETTER than the positive: {better}/{len(g)}")

a, b = lp.best.to_numpy(), ln.best.to_numpy()
roc, p = auroc(a, b)
print(f"\n  overall  n={len(a)} vs {len(b)}   AUROC {roc:.3f}   p={p:.3f}")
print("  (n=3 positives - descriptive, not a powered test)")

print("\n" + "=" * 72)
print("2. REAL PEPTIDES vs SCRAMBLE DECOYS")
print("=" * 72)
rp = per[per["class"] == "real_peptide"]
dc = per[per["class"] == "decoy_scramble"]
roc, p = auroc(rp.best.to_numpy(), dc.best.to_numpy())
print(f"  real n={len(rp)}  median {rp.best.median():+.3f}")
print(f"  scram n={len(dc)}  median {dc.best.median():+.3f}")
print(f"  AUROC {roc:.3f}   p={p:.3f}")
print("  NOTE: scrambles are assumed non-binders, not measured ones.")

print("\n" + "=" * 72)
print("3. DOES THE SCORE JUST TRACK LENGTH?")
print("=" * 72)
rho, p = stats.spearmanr(per.best, per.length)
print(f"  all peptides    rho={rho:+.3f}  p={p:.2g}  n={len(per)}")
for cls in ("real_peptide", "decoy_scramble", "length_negative"):
    g = per[per["class"] == cls]
    if len(g) > 2:
        r, pp = stats.spearmanr(g.best, g.length)
        print(f"  {cls:16s} rho={r:+.3f}  p={pp:.2g}  n={len(g)}")
print("\n  Interpretation: a strongly negative rho means longer peptides score")
print("  better regardless of sequence - the same size bias MM/GBSA showed.")

print("\n" + "=" * 72)
print("4. THE FREE BASELINE: does the sequence start with Tyr?")
print("=" * 72)
lib = per[per["class"].isin(["real_peptide", "decoy_scramble"])]
t = lib[lib.starts_tyr]
nt = lib[~lib.starts_tyr]
print(f"  starts-Y n={len(t)}   median {t.best.median():+.3f}")
print(f"  no-Y     n={len(nt)}  median {nt.best.median():+.3f}")
roc, p = auroc(t.best.to_numpy(), nt.best.to_numpy())
print(f"  ADCP separating Y-start from not: AUROC {roc:.3f}  p={p:.3f}")
tyr_auroc, _ = auroc(
    (~rp.starts_tyr).astype(float).to_numpy() * 0 - rp.starts_tyr.astype(float).to_numpy(),
    (~dc.starts_tyr).astype(float).to_numpy() * 0 - dc.starts_tyr.astype(float).to_numpy())
print(f"  the Tyr RULE ITSELF, real vs scramble: AUROC {tyr_auroc:.3f}")
print("  ADCP must beat that number to be worth its CPU time.")

print("\n" + "=" * 72)
print("5. CONVERGENCE AT FIXED EFFORT (30 replicas, 100k steps each)")
print("=" * 72)
print("  spread = worst-minus-best replica energy, per peptide")
for L in sorted(per.length.unique()):
    g = per[per.length == L]
    print(f"    {L:2d} aa  n={len(g):2d}  median spread {g.spread.median():5.2f}"
          f"   max {g.spread.max():5.2f}")
rho, p = stats.spearmanr(per.spread, per.length)
print(f"\n  spread vs length: rho={rho:+.3f}  p={p:.2g}")
print("  A rising spread means long peptides are NOT converged at this effort,")
print("  which is a statement about cost, not about binding.")

print("\n" + "=" * 72)
print("FULL TABLE (best over replicas)")
print("=" * 72)
print(per.sort_values(["class", "length", "best"])[
    ["name", "class", "length", "sequence", "best", "median_rep", "spread"]
].to_string(index=False))

out = HERE / "docs" / "adcp_per_peptide.tsv"
per.sort_values(["class", "length", "best"]).to_csv(out, sep="\t", index=False)
print(f"\nWrote {out}")
