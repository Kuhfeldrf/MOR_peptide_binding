# Does the co-folding triage separate binders from non-binders?

**No. Not with sequence alone, and not with the receptor pinned to its own
experimental structure.** This is the central negative result of the pilot, and
it is what makes the physics stages mandatory rather than optional.

## The test

Two questions with very different statistical power:

| | comparison | n | power |
|---|---|---|---|
| **Discrimination** | 17 real peptides vs 20 scrambled controls | 37 | good |
| **Calibration** | Stage 2 score vs measured GPI IC50 | 7 | poor |

Discrimination is the one that matters. The decoys are **permutations** of their
real counterparts - identical amino acid composition, identical length - so
composition and length bias cancel exactly, and what remains is sequence
**order**, which is what determines binding.

AUROC is reported because it is literally the probability that a randomly
chosen real peptide outranks a randomly chosen decoy. 0.5 is chance.

## Two arms

| arm | receptor information |
|---|---|
| **A** | sequence only (single-sequence mode, ESM embeddings) |
| **B** | + 6DDF chain R as template, so the receptor is pinned to experiment |

Both arms: same 37 peptides, same 5 seeds, same metrics. The peptide never gets
a template - a 4-7mer has no homologue, and supplying one would hand the model
the answer.

## Result: discrimination absent in both arms

| metric | arm A | arm B |
|---|---|---|
| ipTM | 0.376 | **0.415** |
| aggregate_score | 0.397 | 0.418 |
| pae_interface_min | 0.362 | 0.484 |
| pae_interface_mean | 0.391 | 0.479 |
| i_plddt_peptide | 0.400 | 0.476 |
| peptide_rmsd_mean | 0.512 | 0.503 |

Every value moved **toward** 0.5 under templating and **none crossed it**. The
decoys still score marginally higher than the real peptides (ipTM 0.352 vs
0.332).

Nothing here is individually significant (p = 0.16-0.99), so the honest
statement is **"failed to demonstrate discrimination"** rather than "reliably
inverted". But the direction is consistent across all six metrics and both
arms, and 17-vs-20 is far better powered than the n=7 calibration.

**Giving the model the true receptor did not give it the ability to tell a
binder from a scramble.** That is a capability it lacks, not an information
deficit.

## Calibration improved, but rests on a noisy substrate

| | arm A | arm B |
|---|---|---|
| best rho vs IC50 (n=7) | 0.821 | 0.857 |
| p | 0.023 | 0.014 |
| Bonferroni (6 metrics) | 0.14 | **0.08** |
| jackknife floor | 0.71 | 0.77 |

Still not significant after correction. And see below - the underlying
per-peptide measurements are unstable, so this should not be read as a real
improvement.

## Cross-seed agreement at 5 seeds is too noisy to rank peptides

Aggregate barely moved (real 4.01 → 4.07 A; decoys 4.11 → 4.22 A), but
individual peptides swung violently in **both** directions:

| peptide | arm A | arm B |
|---|---|---|
| met-enkephalin | 1.67 | **1.01** (CONVERGED) |
| beta-casomorphin-7-bov | 2.03 | 1.35 |
| **neocasomorphin-6** | **1.08** | **7.89** |

Neocasomorphin-6 went from the best-converged peptide in the set to badly
divergent, from a change that pins only the **receptor**. A 7 A swing from that
is not a real effect - it means the 5-seed estimate of pose reproducibility has
variance comparable to its range.

**Consequence:** `cross_seed_agreement` cannot be used to rank peptides at
n_seeds = 5, and any metric derived from it inherits that noise.

## What this does and does not say

**Does not say co-folding is useless here.** It is usable for **pose
generation** - met-enkephalin converged to 1.01 A with the template, and
Stage 2 reproduces the orthosteric pocket to 0.80 A against 8F7Q. Generating a
starting structure is a different task from judging whether a peptide binds.

**Does say the confidence scores cannot triage.** ipTM answers *"is my pose
geometrically right?"*, not *"does this bind?"* There is no physical theory
connecting the two, and empirically here there is no relationship either.

**Does say the physics stages are not optional.** The pilot's stated premise
was that learned confidence is not affinity. That is now measured rather than
asserted.

## Consequences for the design

1. **Do not use ipTM or any Chai-1 confidence metric as a screening filter.**
   It would discard real binders at the same rate as decoys.
2. **The triage tier must be physics** - MM/GBSA (Stage 5) at minimum, since
   an approximate energy is at least the right kind of quantity.
3. **Arm C (MSA) is not worth running.** Arm B gave the receptor *perfect*
   information and discrimination did not appear; an MSA supplies strictly less
   than the experimental structure.
4. **Keep co-folding for pose generation**, with templating on, since it
   improves the best cases and costs nothing.
5. **Raise n_seeds or stop using cross-seed agreement as a ranking metric.**

## Caveat on the met-enkephalin row

Met-enkephalin is flagged `HIGH_OVERLAP`: its sequence is the N-terminal
message sequence of beta-endorphin in 8F7Q (deposited 2023), so the model has
very likely seen it. It is the only peptide that converges. Its scores should
not be treated as a prediction.

This is also why **6DDF and not 8F7Q** was used as the template - 8F7Q contains
met-enkephalin's own sequence, and templating on it would pre-shape the pocket
for the one peptide whose result is already contaminated.

Discrimination is robust to this: it compares 17 real against 20 decoys, and
scrambles appear in no deposited structure.

## Data

| file | contents |
|---|---|
| `triage_discrimination_armA_seqonly.tsv` | AUROC, arm A |
| `triage_discrimination_armB_template.tsv` | AUROC, arm B |
| `triage_calibration_armA_seqonly.tsv` | rho vs IC50, arm A |
| `triage_calibration_armB_template.tsv` | rho vs IC50, arm B |
| `triage_per_peptide_arm*.tsv` | per-peptide metric values |
| `stage2_cross_seed_agreement_armB_template.tsv` | arm B convergence |

Generated by `scripts/07_validate_triage.py`; arm B by `jobs/armB.sbatch`.
