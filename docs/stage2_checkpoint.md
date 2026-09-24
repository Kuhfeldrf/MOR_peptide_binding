# Stage 2 checkpoint - co-folding confidence and cross-seed agreement

Three test peptides, five seeds each, Chai-1 single-sequence mode.
15 inference runs producing 75 structures (5 diffusion samples per run).

## Headline finding

**ipTM could not distinguish the three peptides. Cross-seed agreement
separated them cleanly.**

Across a 32-fold affinity range and two pharmacologies, every ipTM value
fell in a narrow band of 0.27-0.35. The same three peptides gave peptide
RMSD across seeds of 1.6 A, 4.2 A and 7.1 A - converged, divergent and very
divergent. **The confidence score was blind to a difference the geometry
made obvious.**

## Per-peptide results

| Peptide | Seq | IC50 (uM) | Role | Overlap | ipTM (min-max) | aggregate | i-pLDDT | seed RMSD mean / max | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| met_enkephalin | YGGFM | 0.2 | agonist | HIGH_OVERLAP | 0.265-0.308 | 0.405 | 0.537 | 1.637 / 1.984 A | **CONVERGED** |
| beta_casomorphin_5_bov | YPFPG | 6.5 | agonist | MODERATE | 0.270-0.310 | 0.400 | 0.532 | 4.173 / 6.803 A | **DIVERGENT** |
| casoxin_C | YIPIQYVLSR | 5 | antagonist | MODERATE/NOVEL | 0.270-0.358 | 0.409 | 0.475 | 7.083 / 11.231 A | **DIVERGENT** |

## Reading this honestly

**1. The one peptide that converged is the one the model has seen.**
Met-enkephalin (YGGFM) is an exact sequence match to 8F7Q, deposited 2023,
and is flagged HIGH_OVERLAP by Stage 2b. It is the only peptide whose five
seeds agree (1.64 A mean, 1.98 A max). The two peptides with no deposited
mu-opioid complex diverge by 4-11 A - different sub-sites, not just
different rotamers.

The most economical reading is that the model reproduces a memorised pose
consistently, and genuinely does not know where the other two bind. That is
exactly the failure mode Stage 2b exists to expose, and it appeared on the
first three peptides tested.

**2. ipTM did not detect this.** Met-enkephalin's ipTM (0.286-0.304) is
indistinguishable from casoxin C's (0.277-0.352), despite one being
converged and the other scattered across 11 A. **A confidence score that
cannot separate a memorised answer from no answer cannot be used to rank
candidates.** This is direct evidence for the pipeline's stated premise:
learned confidence is not affinity, and the physics stages are not optional.

**3. The absolute ipTM values should not be over-interpreted yet.** All are
low (0.27-0.35), but ipTM is known to be diluted when one chain is large and
mostly not at the interface - a 281-residue receptor against a 5-residue
peptide is close to the worst case. **ipSAE exists precisely to correct
this, and is not yet computed (decision D21).** Until it is, these numbers
cannot be read as evidence that the poses are bad - only that ipTM is
uninformative here.

**4. i-pLDDT tracks length, not affinity.** The 10-mer casoxin C scores
lowest (0.459-0.497) and the two 5-mers higher (~0.52-0.54). The ordering
follows peptide length, not the 32-fold affinity range. No evidence of
affinity discrimination.

## Runtime and budget

| Quantity | Value |
|---|---|
| Per inference run (cached weights) | **~78 s** |
| Per run including one-off weight download | 349 s |
| This checkpoint: 3 peptides x 5 seeds | 19.7 min |
| Structures produced | 75 (5 per run) |
| **Projected full library: 40 peptides x 5 seeds** | **~4.3 h on one L40S** |

Comfortably inside budget. The instructions ask for a flag if runtime
exceeds the stated budget by more than ~2x; it does not.

## What this does not show

- Nothing here is validated. Three peptides, one receptor, one model.
- Cross-seed convergence is **not** evidence of correctness. A model can
  converge confidently on a wrong pose; convergence only shows consistency.
- No comparison against an experimental pose has been made. Met-enkephalin's
  converged pose could be checked against 8F7Q, which would say whether it is
  recalling the right answer - that is worth doing and has not been done.

