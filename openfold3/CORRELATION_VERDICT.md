# Does anything from OpenFold3 correlate with IC50?

**Short answer: no — and dropping Met-enkephalin does not change that,
although it makes it look as though it might.**

---

## What was tested

OpenFold3 emits **no binding affinity**. It is a structure predictor. So every
quantity it *does* emit was tested, plus every geometric observable of the
predicted pose — **13 quantities in total**, each with its direction of
agreement declared in advance:

| Confidence (from the model) | Geometric (from the structure) |
|---|---|
| `sample_ranking_score` (OF3's own composite) | anchor distance, best sample |
| `iptm`, `ptm` | anchor distance, median |
| `chain_ptm` for the peptide chain | anchor reproducibility (fraction of 15 samples ≤ 5 Å) |
| `avg_plddt` | pocket contacts, median |
| `gpde`, `disorder` | depth, median |
| interface PAE (computed from the PAE matrix) | |

With 13 tests, the Bonferroni threshold is **α = 0.0038**.

---

## Result 1 — All nine peptides: nothing

| quantity | agrees if | ρ | p |
|---|---|---:|---:|
| sample_ranking_score | ρ < 0 | +0.033 | 0.93 ✗ wrong sign |
| **ipTM** | ρ < 0 | **0.000** | **1.00** |
| pTM | ρ < 0 | −0.017 | 0.97 |
| chain pTM (peptide) | ρ < 0 | −0.050 | 0.90 |
| avg pLDDT | ρ < 0 | −0.150 | 0.70 |
| gPDE | ρ > 0 | +0.350 | 0.36 |
| **interface PAE** | ρ > 0 | **0.000** | **1.00** |
| anchor distance, best | ρ > 0 | +0.433 | 0.24 |
| anchor reproducibility | ρ < 0 | −0.568 | 0.11 |
| pocket contacts | ρ < 0 | +0.391 | 0.30 ✗ wrong sign |
| depth | ρ > 0 | +0.333 | 0.38 |

**Correct-direction quantities with raw p < 0.05: zero.**
**Surviving Bonferroni: zero.**

Note ipTM and interface PAE both land at ρ = **exactly 0.000** — not weak,
absent. OpenFold3's own composite ranking score has the *wrong* sign.

---

## Result 2 — Dropping Met-enkephalin

The rationale is reasonable: Met-enkephalin is 0.2 µM against 6.5–160 µM for
the rest — **32× more potent than the next peptide** — and it is the endogenous
reference rather than a milk-derived peptide. It is a high-leverage point.

Excluding it (n = 8):

| quantity | ρ (n=9) | ρ (n=8) | p (n=8) |
|---|---:|---:|---:|
| **gPDE** | +0.350 | **+0.786** | **0.021** |
| avg pLDDT | −0.150 | −0.643 | 0.086 |
| **pocket contacts** ✗ *wrong sign* | +0.391 | **+0.866** | **0.005** |
| anchor reproducibility | −0.568 | −0.552 | 0.156 |

gPDE now reaches raw p = 0.021. **It still does not survive the
multiple-comparison threshold (0.0038).**

---

## Result 3 — Why that improvement is not real

Leave-one-out across *every* peptide, not just the one of interest:

**gPDE** — ρ with each peptide dropped in turn:

| dropped | ρ | p |
|---|---:|---:|
| *(none)* | +0.350 | 0.36 |
| **met_enkephalin** | **+0.786** | **0.021** |
| β-casomorphin-5 bov | +0.262 | 0.53 |
| β-casomorphin-5 hum | +0.262 | 0.53 |
| β-casomorphin-4 | +0.238 | 0.57 |
| β-casomorphin-7 hum | +0.238 | 0.57 |
| α-lactorphin | +0.357 | 0.39 |
| β-casomorphin-7 bov | +0.238 | 0.57 |
| neocasomorphin-6 | +0.381 | 0.35 |
| β-lactorphin | +0.381 | 0.35 |

**Only Met-enkephalin's removal produces the effect.** Every other exclusion
leaves ρ between 0.24 and 0.38, all p > 0.35. The same pattern holds for
avg pLDDT (−0.150 → −0.643 on Met-enkephalin's removal; −0.26 to +0.02 for
every other).

**The decisive tell:** removing Met-enkephalin *also* strengthens the
**pocket-contacts** correlation to ρ = +0.866, p = 0.005 — and that one runs
**backwards**, implying more orthosteric contacts go with *weaker* binding. If
excluding Met-enkephalin were revealing a masked true signal, it would not
simultaneously make a biologically inverted relationship the most significant
result in the table.

What is actually happening is ordinary high-leverage behaviour: Met-enkephalin
sits far from the others on log IC50, so removing it reshuffles ranks and
inflates whichever correlations happen to align — in both directions.

---

## Result 4 — The one quantity that is stable

**Anchor reproducibility** — the fraction of the 15 sampled models in which the
Tyr→D3.32 salt bridge forms — is the only quantity insensitive to which point
is dropped:

| dropped | ρ |
|---|---:|
| *(none)* | −0.568 |
| met_enkephalin | −0.552 |
| β-casomorphin-5 bov | −0.552 |
| β-casomorphin-5 hum | −0.498 |
| β-casomorphin-4 | −0.503 |
| β-casomorphin-7 hum | −0.520 |
| α-lactorphin | −0.520 |
| β-casomorphin-7 bov | −0.536 |
| neocasomorphin-6 | −0.536 |
| **β-lactorphin** | **−0.847** (p = 0.008) |

It holds between −0.50 and −0.57 under every exclusion but one, in the correct
direction, and **never reaches significance** (p = 0.11 at n = 9).

The β-lactorphin row is the mirror image of the Met-enkephalin temptation:
YLLF is the weakest binder (160 µM) yet has high reproducibility, so removing it
lifts ρ to −0.847. That is dropping the point that disagrees, and it should be
resisted for exactly the same reason.

---

## Verdict

1. **No OpenFold3 quantity — confidence or geometric — correlates with IC50 at
   n = 9.** Zero of thirteen, in the correct direction, at any threshold.
2. **Dropping Met-enkephalin improves two correlations, but the improvement is
   specific to that single exclusion and does not survive correction.** The same
   exclusion strengthens a wrong-signed correlation more than either of them.
3. **Anchor reproducibility is the only stable candidate** — correct direction,
   robust to exclusions, consistently ρ ≈ −0.55, never significant. It is worth
   a prospective test, not a claim.
4. n = 9 is simply too small. Under Spearman at n = 9, ρ must exceed ≈ 0.68 for
   raw significance; after correction for 13 tests, ≈ 0.88. Nothing here is
   close.

**Recommendation:** report the absence of correlation, with the full n = 9
table, as the honest result. If anchor reproducibility is to be pursued, it
needs a pre-registered test on an independent set of ≥ 20 peptides — which is
roughly 2 GPU-hours, and the only way the question gets answered rather than
argued.

---

Generated by `src/analyze_of3_correlations.py`; full tables in
`OF3_CORRELATIONS.md`, per-peptide values in `of3_quantities.csv`.
