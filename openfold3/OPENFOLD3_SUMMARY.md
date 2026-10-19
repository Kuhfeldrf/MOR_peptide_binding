# OpenFold3 experiment — MOR agonists, 8F7Q templated

**Date:** 2026-09-22 · ORCA · OpenFold3 0.5.0 (`of3-ob-2025-06-30-174k.pt`)
**84/84 queries succeeded · 420 models · 0 failures · 1 h 05 m wall, ~1.1 GPU-hours**

> Demonstration output. **Not repository content.**

---

## Design

14 MOR **agonist** peptides (antagonists excluded per instruction, since a single
receptor cannot support an active/inactive call) plus 14 composition-matched
scrambles. Receptor chain templated on **8F7Q** via `template_cif_paths`;
3 seeds × 5 diffusion samples = **15 models per peptide**.

**Two constraints worth stating up front:**

1. **Direct pocket targeting is unavailable for peptides.** OpenFold3's
   `PocketConstraint` requires `ligand_chain_id` to name a chain of type
   `LIGAND`; a peptide supplied as `PROTEIN` is rejected — the same restriction
   that blocks Boltz-2's affinity head. Conditioning was therefore achieved by
   **templating the receptor** on 8F7Q, not by constraining the peptide.
2. **OpenFold3 produces no binding affinity.** It is a structure predictor. What
   is correlated with IC50 below is (a) its confidence and (b) *geometric
   measurements of the predicted pose* — the latter being a physical observable
   rather than a model self-assessment.

**Numbering.** OpenFold3 renumbers output chains 1..N, so canonical Asp149 (8F7Q
frame) appears at residue **84**. A first measurement read residue 149 — which is
**GLN** — and reported "0–1 of 5 pocket contacts, anchor n/a". Every number here
comes from `src/measure_pose.py`, which auto-detects the offset by residue
identity (found: **65**) and refuses to measure if no unique offset validates.

---

## Result 1 — The pose works. Very well.

Every one of the 14 real agonists formed the Tyr→D3.32 salt bridge (≤5.0 Å) in
at least one sample, with 8–10 of 10 orthosteric pocket contacts.

| peptide | seq | len | IC50 (µM) | anchor best | anchor med | pocket med | anchor formed |
|---|---|---:|---:|---:|---:|---:|---:|
| met_enkephalin | YGGFM | 5 | 0.2 | **2.03 Å** | 5.49 | 10 | 7/15 |
| beta_casomorphin_5_bov | YPFPG | 5 | 6.5 | 3.26 | 5.46 | 8 | 7/15 |
| beta_casomorphin_5_hum | YPFVE | 5 | 14 | 2.18 | 4.59 | 9 | 8/15 |
| beta_casomorphin_4_bov | YPFP | 4 | 22 | 3.37 | 5.18 | 9 | 7/15 |
| beta_casomorphin_7_hum | YPFVEPI | 7 | 25 | 2.62 | 5.66 | 9 | 5/15 |
| alpha_lactorphin | YGLF | 4 | 50 | 3.06 | 6.94 | 9 | 5/15 |
| beta_casomorphin_7_bov | YPFPGPI | 7 | 57 | 3.99 | 8.23 | 9 | 1/15 |
| neocasomorphin_6 | YPVEPF | 6 | 59 | 3.52 | 12.10 | 10 | 2/15 |
| beta_lactorphin | YLLF | 4 | 160 | 2.54 | 5.01 | 10 | 7/15 |

For context: the **experimental** Met-enkephalin pose measures 4.6 Å at this
contact. OpenFold3's best sample reaches **2.03 Å**, and it does so from
sequence in ~12 s per model — against 4000 rigid-body models and a 1-in-3 seed
success rate in the prior HADDOCK3 work.

**This is a genuine technical success and the strongest result in the project
so far.** It is a statement about *pose generation*, not about potency.

---

## Result 2 — It does not rank potency (n = 9)

| quantity | agreement direction | Spearman ρ | p | partial ρ (length controlled) | p |
|---|---|---:|---:|---:|---:|
| anchor distance, best | positive | +0.433 | 0.244 | +0.447 | 0.228 |
| anchor distance, median | positive | +0.350 | 0.356 | +0.404 | 0.281 |
| pocket contacts, median | negative | **+0.391** ✗ | 0.298 | +0.394 | 0.294 |
| **ipTM, median** | negative | **0.000** | **1.000** | +0.019 | 0.962 |
| depth, median | positive | +0.333 | 0.381 | +0.539 | 0.134 |
| **anchor reproducibility** | negative | **−0.568** | 0.110 | **−0.665** | **0.051** |

Nothing reaches significance at n = 9.

- **ipTM has a Spearman ρ of exactly 0.000 with IC50.** Not weak — absent.
  Consistent with the Boltz-2 result, in a different model with a different
  confidence calibration.
- **Pocket contacts run the wrong way** (positive where agreement needs
  negative).
- **Anchor reproducibility** — the fraction of the 15 samples that form the
  salt bridge — is the only quantity with a correct-direction, appreciable
  effect: ρ = −0.568, rising to **−0.665 (p = 0.051)** once length is
  controlled.

**That last row must not be oversold.** Six quantities were tested; a
Bonferroni-corrected threshold would be α ≈ 0.008, and p = 0.051 does not
survive it. It is a **hypothesis worth testing prospectively**, not a finding.
It is also plausible mechanistically — how *reliably* a model re-forms a
pharmacophore contact is a different thing from how well it can do so once.

---

## Result 3 — Scrambles pose just as well

Best-sample anchor distance, paired:

- **Real peptides anchor closer in only 5 of 14 pairs** — worse than chance
- Wilcoxon signed-rank **p = 0.626**
- Mean difference +0.62 Å, driven entirely by two outliers
  (β-casomorphin-7 hum +9.30 Å, soymorphin-6 +5.20 Å); the median is negative

Scrambled controls reach the orthosteric site with the salt bridge formed as
readily as the real motifs. Examples: scrambled soymorphin-7 anchors at
**1.67 Å**, scrambled β-lactorphin at **1.91 Å** — both tighter than
Met-enkephalin's real 2.03 Å.

**Interpretation:** the pose is achievable by essentially any Tyr-containing
peptide of this size. Templating the receptor on a peptide-bound structure
supplies a pocket that will accept almost anything presented to it. Pose
formation is therefore **not evidence of binding**, and cannot be used as one.

---

## Result 4 — The length confound, replicated independently

- ρ(ipTM, length) = **−0.736 (p = 0.024)** — shorter peptides score higher

This reproduces the background screen's ρ = −0.654 in a **different model**, on
a **different peptide set**. The confound is a property of co-folding confidence
generally, not a Boltz-2 artefact.

Usefully, on this particular subset ρ(length, log IC50) = **0.017** — length and
potency are unrelated here. So the failure of ipTM to track IC50 is *not*
explained away by length; ipTM tracks length, and length tracks nothing.

**The 4-mer test your IC50 table made possible:** β-lactorphin (YLLF) is the
weakest binder in the set at **160 µM** — 800× weaker than Met-enkephalin — yet
it achieves a **2.54 Å** anchor in 7 of 15 samples, out-posing peptides that
bind 20× more tightly. That single comparison makes the point more plainly than
any correlation coefficient.

---

## Conclusions

**What OpenFold3 + 8F7Q templating does well:** generate physically credible
orthosteric poses with the correct pharmacophore contact, fast, from sequence,
reproducibly enough to be useful as a *structure generator*.

**What it does not do:** distinguish strong binders from weak ones, or real
motifs from scrambles. No available scalar — confidence or geometric —
correlates with potency at this sample size.

This is the same dissociation the prior HADDOCK3 diagnostic work identified and
scoped precisely: *search finds the pose; scoring cannot see the register.* A
different stack, a different scoring function, the same conclusion. That
consistency across independent methods is itself the most defensible finding
here.

---

## Recommended next steps

1. **Test anchor reproducibility prospectively.** It is the only correct-direction
   candidate. A pre-registered test on an independent peptide set with n ≥ 20
   would settle it. ~2 GPU-hours.
2. **Do not use pose formation as evidence of binding.** Scrambles pose equally
   well; any screen built on "does it reach D3.32" will not enrich.
3. **Physics-based rescoring** (MM-GBSA/FEP) on these generated poses is now
   more attractive than before: OpenFold3 supplies high-quality starting
   structures cheaply, which is exactly what such a rescore needs. This is the
   natural division of labour — generation by the neural model, ranking by
   physics.
4. **Trace the IC50 values to primary literature.** Nine values now carry the
   entire correlation analysis; they currently trace to a prior analysis script
   and a user-supplied table.

---

## Artifacts

`OF3_ANALYSIS.md` (generated), `of3_per_peptide.csv` (28 peptides),
`of3_samples.csv` (420 models, every measurement). Structures remain on the
cluster at `~/bionemo_mor_run/of3/agonist_out` (1.7 GB).
