# uOR peptide screening pilot

A seven-stage pipeline that takes a list of peptide sequences and a mu-opioid
receptor structure, predicts each peptide-receptor complex with a co-folding
model, and then carries the complexes through membrane molecular dynamics and
physics-based free-energy scoring. The final stage compares both families of
score against published binding affinities for a curated reference set. The
purpose is to show that every stage exists, executes, and is reproducible - not
to produce a screening result.

> **This repository is evidence of technical readiness, not a results campaign.**
> Read the "Validated versus demonstrated" section before citing any number in it.

---

## Stage status

Runtime and hardware are recorded only for stages that have actually run.

| Stage | Description | Status | Runtime | Hardware |
|-------|-------------|--------|---------|----------|
| 0 | Receptor preparation (6DDF) | `STUBBED` | - | - |
| 1 | Peptide library ingestion | `STUBBED` | - | - |
| 2 | Co-folding (Chai-1, 5 seeds) | `STUBBED` | - | - |
| 2b | Training-overlap audit | `STUBBED` | - | - |
| 3 | Membrane system build | `STUBBED` | - | - |
| 4 | Molecular dynamics (GROMACS) | `STUBBED` | - | - |
| 5 | MM/GBSA triage | `STUBBED` | - | - |
| 6 | ABFE (DAMGO first) | `STUBBED` | - | - |
| 7 | Benchmark figure | `STUBBED` | - | - |

**Every stage is currently stubbed.** The rule graph, inputs, and outputs are
real and wired together; each stub script exits non-zero with a `STUBBED`
marker rather than writing an empty or fabricated output, so a stubbed stage
cannot be mistaken for one that ran.

---

## Validated versus demonstrated

Nothing in this repository is validated. Nothing in it has yet been
demonstrated either - no stage has run.

When stages do run, the distinction that will apply is:

- **Demonstrated** - the stage executed end to end on the three-peptide test set
  and produced output of the expected shape. This says the plumbing works.
- **Validated** - the stage's output has been shown to agree with an independent
  measurement across a set large enough to support the claim.

A stage that has run once on three peptides is **demonstrated, not validated**.
The pilot as scoped can reach "demonstrated" for all seven stages and
"validated" for none of them, because the reference set carries affinity values
for only six peptides (see below). Any README revision that blurs this line
should be treated as a regression.

---

## Reference data and its limits

Stage 1 reads a real, curated peptide set copied verbatim from prior ORCA work
(`data/reference/`): 20 food-derived opioid peptides plus 20 scrambled
composition-matched negative controls.

Two constraints govern how Stage 7 may use it:

1. **The affinity column is IC50 in micromolar, not Ki in nanomolar**, and only
   **6 of 20** peptides carry any value. No Cheng-Prusoff conversion is applied,
   because the source does not record the radioligand Kd or concentration that
   such a conversion requires. `ki_nM` and `ic50_um` are separate columns.
2. **Those six values are currently marked `UNSOURCED`.** The source attributes
   them to an analysis script, not a paper. A script reference does not satisfy
   the sourcing requirement. They must be traced to primary literature - with
   radioligand, assay system, and species - before Stage 7 treats them as a
   validation target.

With n = 6, a Spearman coefficient is weakly determined. That belongs in the
figure caption, not only in the methods.

---

## Reproducing the test path

The test profile shortens MD production and reduces ABFE lambda windows. It
never alters production parameters: overrides live in `test_overrides` in
`config/config.yaml` and are applied once, in the Snakefile.

```bash
export PILOT_ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$PILOT_ROOT"

source /etc/profile.d/z00_lmod.sh      # module is undefined otherwise
module load apptainer cuda/12.9.0

conda activate mor-pilot

snakemake --profile config/slurm --config profile=test -n   # dry run
snakemake --profile config/slurm --config profile=test
```

Target: three peptides through all seven stages in under one hour on a single
GPU node. **Not yet achieved - every stage is stubbed.**

---

## Known limitations

- **No experimental validation is in scope.** Nothing here is measured.
- **Co-folding accuracy depends on training-set overlap.** Reference peptides
  will predict well partly through memorisation. Stage 2b exists to flag this,
  and Stage 7 stratifies by it; performance on `NOVEL` peptides is the honest
  measure.
- **Co-folding confidence is not affinity.** ipTM, i-pLDDT, and ipSAE are the
  model's estimate of its own accuracy. A confidently ranked top pose can be
  confidently wrong. This is the stated reason the pipeline continues into
  physics.
- **ABFE convergence is the known weak point** for flexible, often charged
  peptides. Stage 6 reports a convergence assessment, not a bare number.
- **Only 6 reference peptides carry an affinity value, and it is IC50 not Ki,
  and it is currently unsourced.** This is the binding constraint on Stage 7.
- **DAMGO is extracted from experiment, not predicted.** It is
  Tyr-D-Ala-Gly-N-MePhe-Gly-ol - a D-amino acid, an N-methylated residue and a
  C-terminal alcohol - so plain FASTA input to Chai-1 cannot express it. It is
  instead taken from 6DDF, where it is experimentally resolved as entity 5, and
  is excluded from sequence-based co-folding. **The consequence is that DAMGO's
  ABFE result validates the free-energy machinery only - decoupling, restraints,
  MBAR, convergence - and does not test pose prediction**, because it starts
  from the experimental complex. See `docs/damgo_notes.md`.
- **`MEA` and `ETA` (N-methyl-Phe and the Gly-ol cap) have no standard protein
  force-field parameters** and need GAFF2 treatment via antechamber. The charge
  derivation scheme is not yet decided. This is the live technical risk on the
  DAMGO path.
- **ref2015-based tools are absent by design.** FlexPepDock and related Rosetta
  refinement are excluded because ref2015 has no membrane term; membrane MD
  performs the equivalent role in the correct physical environment.
- **Amber-to-GROMACS topology conversion at the Stage 3/4 boundary** is a real
  step and a genuine place for silent error. Stage 3 must verify atom count and
  total charge across the conversion.
- **The OSU target architecture is unconfirmed**, and the engine choice depends
  on it. See `docs/arch_notes.md`.

---

## Scaling beyond the pilot

What a full receptor panel and peptidome would require, to be quantified once
the pilot has produced real timings:

- **Storage.** The trajectory policy is the governing constraint. Stage 4 writes
  a reduced trajectory (receptor + peptide + proximal lipids) at analysis
  frequency and full-system frames sparsely. Stage 6 saves energies and
  subsampled coordinates only: full trajectories across all lambda windows would
  be 30-50 GB *per peptide*, versus 2-5 GB for what MBAR actually needs.
  Uncontrolled trajectory writing is the single most likely way to fill the
  filesystem.
- **GPU hours.** Chai-1 cost scales roughly quadratically with total token
  count. Stage 6 dominates: ~36 windows x several ns each, per peptide.
- **Porting work.** See `docs/arch_notes.md`. If the OSU system is aarch64, every
  Python dependency needs a wheel-availability check; GROMACS was chosen to keep
  the MD engine itself out of that risk category.

---

## Layout

```
config/      global parameters + Snakemake SLURM profile
workflow/    Snakefile and one .smk rule per stage
scripts/     one script per rule; no logic lives in the Snakefile
data/        raw PDBs, curated reference set, MS input
test/        three-peptide test set
docs/        provenance.md (versions, builds) + arch_notes.md (porting)
results/     gitignored
logs/        gitignored
```

## Provenance

Cluster details, module versions, and build flags: `docs/provenance.md`.
Architecture and engine rationale: `docs/arch_notes.md`.

**Note:** the build instructions describe ORCA as an Oregon State University
cluster. The host in use is `login.orca.pdx.edu` (Portland State University).
The discrepancy is recorded rather than resolved.
