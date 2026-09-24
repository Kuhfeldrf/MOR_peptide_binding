# Prior work: what this replaces, and what carries forward

## What is being replaced

This pipeline replaces a **HADDOCK3 / Uni-Dock / Rosetta** docking stack with a
GPU structure-prediction stack built on the open-source NVIDIA BioNeMo
framework (MSA-Search, OpenFold2/3, Boltz-2).

Prior public repositories:

* [`Kuhfeldrf/hpc.cqls`](https://github.com/Kuhfeldrf/hpc.cqls) — Slurm arrays,
  dependency chaining, GROMACS MD, plotting conventions. The orchestration
  patterns in `slurm/` follow it.
* [`Kuhfeldrf/peptide-md-docking`](https://github.com/Kuhfeldrf/peptide-md-docking)
  — the Rosetta peptide docking pipeline.

**A caveat on those repositories, recorded because it caused a problem before:**
neither archives the MOR opioid work. They archive a Rosetta pipeline docking
two milk peptides against TLR4 (3FXI). Earlier proposal text cited them as the
archive for the MOR workflow; they are not. **This repository exists in part to
be the archive that claim actually needs.**

## Baseline timings — what the new stack is measured against

From the prior Rosetta stack, per peptide:

| Stage | Wall time |
|---|---|
| Global docking | ~2–4 h |
| FlexPepDock refinement | ~1–2 h |

HADDOCK3, from the completed July 2026 run: rigidbody sampling alone
(4000 models) took **~2 h 48 min on 16 cores**, and the whole-node
`--exclusive` 64-core configuration was adopted to make the array finish
inside a walltime at all.

Measured timings for the new stack are reported from the demo run, not
estimated here.

## What carries forward unchanged

So that the two stacks are **directly comparable**, the following are reused
verbatim rather than regenerated:

| Carried forward | Source | Note |
|---|---|---|
| Peptide set (20) | `preliminary_agonist_binding_data_v2/peptides.py` | sequences, names, agonist/antagonist class |
| Scrambled controls (20) | `.../scrambles.tsv` | reused, not regenerated; composition verified against each parent |
| Seeds | `917`, `424242`, `20260715` | same three replicates |
| Receptor structures | `updated_docking_Met_Enkephalin_to_Agonist_MOR/` | 5C1M active, 4DKL inactive, already prepared |
| Pose-validation logic | `validate_pose.py` / `validate_batch.py` | the D3.32-anchor / pocket-contact / depth criteria |
| IC50 values (6) | `anganost.py` | the only published potencies available |

## Lessons from the prior runs, and where each one landed in this code

The prior work is unusually well documented about its own failures. Each
lesson below is either implemented here or recorded in `ASSUMPTIONS.md`.

**1. Always use `auth_seq_id`, never `label_seq_id`.**
An earlier preparation renumbered MOR by −51, so canonical Asp147 became
"Asp96" and file residue 147 was really Ile198. Every canonical-numbered
restraint then pointed at the wrong residue, silently.
→ `run_msa.sequence_from_pdb()` reads author numbering only, and
`verify_residues()` **aborts the run** if a declared binding-site residue is
not the expected type. Every receptor declares its own `d332`; nothing
hardcodes 147.

**2. A run with no ground truth cannot tell you it failed.**
The original Met-enkephalin docking ran `caprieval` with no
`reference_fname`, so no RMSD to a real pose was ever computed. The shipped
pose sat **14.3 Å** from its target while the caption claimed the opposite.
→ The experimental Met-enkephalin pose (8F7Q chain B, residues 1–5) is
retained in `config/receptors.yaml` as a positive control, and Met-enkephalin
is in the screened set as the endogenous reference agonist.

**3. Restraints cannot demonstrate discovery.**
The original run used restraints naming D3.32 as the target — it told the
peptide where to sit, and *still* missed by 14.3 Å. A pipeline that supplies
the answer cannot demonstrate that it can find the answer.
→ Structure prediction from sequence uses no orthosteric restraints at all.
Nothing in `config/` names the pocket as a target; `binding_site_residues` is
used only to *check* a predicted pose, never to steer one.

**4. Probe the pocket before docking into any new receptor.**
Rigid 5C1M with BU72 removed is sealed to anything larger than water.
→ Recorded as `rigid_docking_competent: false` on the active receptor. Its
scope — it constrains rigid docking, not structure prediction — is discussed
in `ASSUMPTIONS.md` §2.2, including what remains unverified.

**5. Measure the noise floor before claiming a class difference.**
This is the prior stack's central negative result. Three seeds of the same
peptide against the same receptor gave a **score SD of ≈29** and an
anchor-distance SD of **8.4 Å**. Every agonist-vs-antagonist contrast measured
(Cohen's d from −0.22 to +0.64) had a 95% CI spanning zero and was **an order
of magnitude smaller than that noise floor**, at n = 14 vs 6. The honest
conclusion was that *no classification potential was demonstrated*.
→ `rank.py` computes the seed-to-seed noise floor from the data, warns when
`differential_margin` sits below it, and flags every affected peptide
`below_noise_floor`. `validate_against_known.py` reports the same ruler
alongside every effect size. The margin in the config is a **placeholder to be
set from the measured floor**, not a constant.

**6. Selection, not sampling, was the bottleneck.**
Near-native poses were found (i-RMSD 0.99 Å, DockQ 0.64) but not reliably
*selected*: FCC clustering rank #1 disagreed with best-by-score in **20/20**
runs and had 0/10 pocket contacts every time. More sampling would have
reproduced the same silent failure at higher cost.
→ A structure predictor returns ranked samples with confidence directly, so
the clustering failure mode does not arise in the same form. Whether its
confidence ranking is any better at resolving the D3.32 register is an open
question this demo can begin to answer, and is not assumed.

**7. Count what survives filtering before buying compute.**
The ≥5-residue filter silently dropped four 4-mers — all with published
activity, and one of them (casoxin B) from a six-member antagonist class.
→ The filter is applied but **never silent**: excluded peptides are listed by
name in every `--dry-run`, and `validate_against_known.py --dry-run` reports
the surviving class counts (11 agonists, 5 antagonists) and warns that this is
underpowered. See `ASSUMPTIONS.md` §3.1 for the filter's true provenance.

**8. The scrambled controls were never run.**
The prior run dropped them for time, so the claim "real peptides beat
length-matched nonsense" was never tested.
→ Carried forward and **screened in this demo**. Test 3 of
`validate_against_known.py` is the control the prior stack never ran.

## Relationship to the classification result

`classification_RESULTS.md` on ORCA asked whether agonist/antagonist status
could be read off HADDOCK binding scores. Its answer was **no** — no signal
beyond peptide length and run-to-run noise, underpowered at 14 vs 6.

This repository asks the same question of a different stack. It does **not**
assume the answer will be different. The scientific point of screening both
conformational states is preserved in the config schema and the output table
because it is the hypothesis under test, not because it is established.
