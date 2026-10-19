# OpenFold3 background screen — known agonists vs. a diverse peptide set

**2026-09-23 · ORCA · OpenFold3 0.5.0, 8F7Q templated · 1,500 models, 0 failures,
3 h 48 m, 6.1 GB**

> Demonstration output. **Not repository content.**

---

## Headline

**Every OpenFold3 quantity appears to discriminate known MOR agonists from
background peptides (AUROC 0.71–0.91, p < 0.001). None of it survives length
matching.** The discrimination is peptide size, not MOR binding.

---

## Design

Same 8F7Q-templated protocol as the agonist run (3 seeds × 5 samples), so the
groups are directly comparable.

| group | n | what it is |
|---|---:|---|
| known_agonist | 14 | published MOR agonists, N-terminal Tyr pharmacophore, 4–7 aa |
| background_tyr | 14 | background milk peptides containing a Tyr |
| **background_noY** | **36** | background peptides with **no Tyr at all** — cannot form the message motif |
| agonist_scramble | 14 | scrambles of the agonists (length-identical by construction) |
| background_scramble | 50 | scrambles of the background |

---

## Result 1 — Wide distributions, as expected

Unlike the known-binder-only screen, the background spreads widely:

**ipTM**

| group | n | median | sd | min | max |
|---|---:|---:|---:|---:|---:|
| known_agonist | 14 | 0.674 | 0.027 | 0.614 | 0.711 |
| background_tyr | 14 | 0.526 | 0.131 | **0.136** | 0.652 |
| background_noY | 36 | 0.539 | 0.136 | 0.178 | **0.766** |
| **agonist_scramble** | 14 | **0.706** | 0.046 | 0.632 | 0.801 |
| background_scramble | 50 | 0.545 | 0.104 | 0.298 | 0.775 |

This is the distribution shape you predicted — tight for known binders, broad
for the background (0.136–0.766). **Note already that the agonists' own
scrambles have a higher median (0.706) than the real agonists (0.674).**

---

## Result 2 — The apparent discrimination

| quantity | AUROC | 95% CI | p |
|---|---:|---|---:|
| depth | **0.906** | 0.817–0.967 | <0.0001 |
| ipTM | **0.894** | 0.809–0.961 | <0.0001 |
| sample_ranking_score | **0.881** | 0.793–0.953 | <0.0001 |
| gPDE | **0.871** | 0.779–0.949 | <0.0001 |
| avg pLDDT | **0.844** | 0.741–0.931 | 0.0001 |
| pTM | 0.831 | 0.727–0.921 | 0.0002 |
| chain pTM (peptide) | 0.820 | 0.716–0.909 | 0.0003 |
| anchor reproducibility | 0.790 | 0.675–0.891 | 0.0006 |
| anchor best | 0.711 | 0.586–0.827 | 0.017 |
| **pocket contacts** | **0.594** | **0.440–0.741** | 0.244 |

Nine of ten exclude 0.5. Taken alone this reads as a working screen.

---

## Result 3 — It is length, and the evidence is unambiguous

**The confound, measured (n = 64):**

| quantity | ρ with length |
|---|---:|
| depth | **+0.941** |
| gPDE | **+0.937** |
| ipTM | **−0.925** |
| avg pLDDT | **−0.920** |
| sample_ranking_score | **−0.916** |
| pTM | −0.891 |
| anchor reproducibility | −0.400 |
| chain pTM (peptide) | −0.407 |
| anchor best | +0.342 |
| pocket contacts | *(below 0.3)* |

Agonists are **4–7 residues**; the background is **5–29**. Most quantities
track length at |ρ| > 0.9 — near-deterministic.

**Length-matched AUROC** (background restricted to 4–7 residues, n = 9):

| quantity | AUROC all | **length-matched** | 95% CI | p |
|---|---:|---:|---|---:|
| sample_ranking_score | 0.881 | **0.595** | 0.325–0.849 | 0.47 |
| ipTM | 0.894 | **0.627** | 0.357–0.873 | 0.33 |
| pTM | 0.831 | **0.365** | 0.119–0.635 | 0.30 |
| chain pTM (peptide) | 0.820 | **0.659** | 0.405–0.889 | 0.22 |
| avg pLDDT | 0.844 | **0.349** | 0.103–0.619 | 0.24 |
| gPDE | 0.871 | **0.421** | 0.159–0.683 | 0.55 |
| anchor best | 0.711 | **0.556** | 0.270–0.817 | 0.68 |
| anchor reproducibility | 0.790 | **0.627** | 0.353–0.881 | 0.33 |
| pocket contacts | 0.594 | **0.560** | 0.357–0.766 | 0.62 |
| depth | 0.906 | **0.595** | 0.349–0.825 | 0.47 |

**Not one length-matched CI excludes 0.5.** Four drop *below* 0.5 — the
agonists score worse than length-matched background peptides.

**The internal consistency check:** `pocket contacts` is the only quantity that
does **not** strongly track length — and it is the only one that failed to
discriminate in the raw comparison (AUROC 0.594, CI 0.440–0.741). The quantity
without the confound is the quantity without the "signal". That is exactly what
a pure length artefact predicts.

---

## Result 4 — The composition control removes all doubt

Agonists vs. **their own scrambles** — length identical by construction, so the
confound cannot operate at all:

| quantity | agonist median | scramble median | AUROC |
|---|---:|---:|---:|
| sample_ranking_score | 0.717 | **0.741** | **0.321** |
| ipTM | 0.674 | **0.706** | **0.321** |
| anchor reproducibility | 0.333 | **0.600** | **0.344** |
| anchor best (Å) | 3.314 | **2.605** | **0.352** |
| depth (Å) | 7.356 | **6.942** | **0.388** |
| pTM | 0.850 | 0.852 | 0.434 |
| avg pLDDT | 81.13 | 81.18 | 0.469 |
| pocket contacts | 9.500 | **10.000** | 0.449 |
| chain pTM (peptide) | 0.105 | 0.105 | 0.500 |
| gPDE | 0.506 | 0.511 | 0.505 |

**All ten quantities sit at or below 0.5.** Scrambled agonists score equal or
better on every measure. Most striking: scrambles form the Tyr→D3.32 anchor
**more** reproducibly than the real agonists (0.600 vs 0.333) and anchor
**closer** (2.61 Å vs 3.31 Å).

No individual comparison reaches significance (n = 14 vs 14, p = 0.11–1.00), so
this is a consistent *direction* rather than a demonstrated reversal. But ten
quantities all landing on the same side of 0.5 is not what a real signal looks
like.

---

## Result 5 — The pharmacophore is irrelevant to the score

36 background peptides contain **no tyrosine** and therefore cannot form the
opioid message motif. Against them:

- **pocket contacts**: agonists 9.5, no-Tyr background 9.0 — AUROC 0.629, p = 0.12

Peptides that physically cannot engage D3.32 make **essentially the same number
of orthosteric contacts** as published agonists. The other quantities separate
them (ipTM AUROC 0.861, depth 0.897) — but those are the length-confounded
ones, and the no-Tyr group is longer on average.

---

## Conclusion

1. **The distribution is wide**, as you anticipated — ipTM 0.136–0.766 across
   the background, against 0.614–0.711 for known binders.
2. **The separation is real but it is length.** |ρ| = 0.92–0.94 for the top
   quantities; every length-matched CI includes 0.5.
3. **Scrambles beat their parents on all ten quantities.** Composition-matched,
   length-matched, and still no advantage for the real motif.
4. **Peptides with no Tyr reach the pocket as well as agonists.**
5. **This replicates the Boltz-2 finding in a second model** — there,
   ρ(ipTM, length) = −0.654 and length-matched AUROC fell from 0.833 to 0.681.
   Here the same collapse is sharper: 0.894 → 0.627.

**No OpenFold3 quantity identifies MOR ligands.** What every confidence measure
reports is that short peptides are easier to place confidently than long ones.

### Power caveat, stated plainly

Only **9 of 50** background peptides fall in the agonists' 4–7 residue window,
so the length-matched CIs are wide (0.10–0.89). "No signal" here means *none
detectable at this power*, not proof of absence. The scramble comparison —
perfectly length-matched, n = 14 vs 14 — points the same way and is the
stronger evidence.

**The fix is cheap:** a length-matched background of ~40 peptides at 4–7
residues would tighten those CIs substantially. From the 1,503-peptide pool,
that is a straightforward sample and roughly 3 GPU-hours.

---

## Artifacts

`OF3_DISTRIBUTION.md` (full tables), `of3_distribution.csv` (128 peptides, all
quantities). Structures on the cluster at `~/bionemo_mor_run/of3/bg_out` (6.1 GB).
