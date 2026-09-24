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
| 0 | Receptor preparation (6DDF) | `WORKING` | ~40 s | login node, CPU |
| 1 | Peptide library ingestion | `WORKING` | < 1 s | login node, CPU |
| 2 | Co-folding (Chai-1, 5 seeds) | `WORKING` | ~78 s/run | 1x L40S |
| 2b | Training-overlap audit | `WORKING` | < 1 s | CPU |
| 3 | Membrane system build | `IN PROGRESS` | ~20-40 min | CPU, 8 cores |
| 4 | Molecular dynamics (GROMACS) | `STUBBED` | - | - |
| 5 | MM/GBSA triage | `STUBBED` | - | - |
| 6 | ABFE (DAMGO first) | `STUBBED` | - | - |
| 7 | Benchmark figure | `STUBBED` | - | - |

**Stages 0, 1, 2 and 2b run. Stage 3 is in progress. Stages 4-7 are stubbed.** The rule graph, inputs, and
outputs are real and wired together throughout; each remaining stub exits
non-zero with a `STUBBED` marker rather than writing an empty or fabricated
output, so a stubbed stage cannot be mistaken for one that ran.

### Stage 0 findings

Two expectations in the build instructions did not survive contact with the
structure, and are recorded rather than worked around:

- **There is no nanobody and no scFv16 in 6DDF.** Only the Gi heterotrimer
  (chains A, B, C) is present alongside the receptor (R) and DAMGO (D). The
  instruction to strip them is a no-op here.
- **The receptor has no chain breaks.** Residues 65-345 are continuous, so
  **ECL2 is fully resolved** and no loop modelling is performed. The
  instructions call for modelling missing residues "particularly ECL2"; there
  are none. This is a property of the structure, not a skipped step.

Protonation at pH 7.4 (PROPKA via pdb2pqr, AMBER naming), for the two residues
the instructions name explicitly:

| Residue | State | Evidence |
|---------|-------|----------|
| **D147** | deprotonated, -1 | no `HD2` on OD2 |
| **H297** | **HID** - neutral, delta-protonated | `HD1` present, `HE2` absent |

Determined from the hydrogens actually placed rather than from residue names,
because `pdb2pqr --pdb-output` retains the generic `HIS`/`ASP` names and cannot
distinguish HID/HIE/HIP on the name alone.

**Determinism verified:** two runs from a clean directory produce byte-identical
`mOR_clean.pdb`, `damgo_ref.pdb`, `mOR_clean.pqr` and `mOR_chainR_raw.pdb`.
SHA-256 digests are written to `results/00_receptor/checksums.json`. Timestamps
are confined to the run log and excluded from the hashed outputs.

### Stage 1 output

40 peptides ingested: 20 real food-derived opioid peptides and 20
composition-matched scrambled controls. 13 carry a literature-sourced IC50.
The GPI benchmark subset is **n = 10** (7 agonists, 3 antagonists), spanning
0.2-200 uM. No peptide was dropped by the filters.

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

## Stage 2 findings: confidence is not discrimination

Three test peptides, five seeds each. Full report in
`docs/stage2_checkpoint.md`.

**ipTM could not tell the peptides apart; cross-seed agreement could.** Across
a 32-fold affinity range and both pharmacologies, every ipTM fell in a narrow
0.27-0.35 band. The same peptides gave cross-seed peptide RMSDs of 1.6, 4.2 and
7.1 A - converged, divergent and very divergent.

The one peptide whose seeds converged, Met-enkephalin, is the only
`HIGH_OVERLAP` case: its sequence is the N-terminal message sequence of
beta-endorphin in 8F7Q. Validated against that structure, Chai-1 reproduces the
orthosteric pocket at **0.80 A** and places the peptide within **2.32 A** of
experiment. So it is both consistent *and* correct - for the one peptide it has
seen. The two with no deposited muOR complex scatter by 4-11 A, which is a
different sub-site rather than a different rotamer.

**This is the pipeline's premise demonstrated rather than asserted:** a learned
confidence score that cannot separate a memorised answer from no answer cannot
rank candidates, and the physics stages are not optional.

Caveat carried forward: ipTM is known to be diluted when one chain is large and
mostly off-interface, and 281 residues against 5 is close to the worst case.
**ipSAE corrects exactly this and is not yet computed** (decision D21), so the
low absolute values show that ipTM is uninformative here, not that the poses
are bad.

## Membrane model: what it represents, and what it does not

The receptor is muOR in **enteroendocrine cells of the intestine**, so the
bilayer is chosen for that environment. Full reasoning in decisions D26/D27.

**Intestinal epithelial membranes are strongly asymmetric**, so "intestinal"
does not by itself specify a composition:

| Face | cholesterol : phospholipid : glycolipid | Cholesterol |
|------|------------------------------------------|-------------|
| Apical (brush border) | ~1 : 1 : 1 | ~50 mol%, glycosphingolipid-rich |
| Basolateral | ~1 : 2.5 : 0.3 | ~28-29 mol% |

**The basolateral/neuronal case is the right target.** Intestinal muOR is
reported mainly on **enteric neurons** of the myenteric and submucosal plexus,
with epithelial expression tied to basolateral function. Decisively, **every
IC50 in the benchmark comes from the GPI assay**, which is guinea-pig ileum
longitudinal muscle / myenteric plexus.

**"Enteric neurons" means intestinal tissue, not brain.** The enteric nervous
system is roughly 500 million neurons embedded in the gut wall. The GPI
preparation is gut throughout; the distinction is which cell type within it -
neurons in the muscle layer rather than epithelium lining the lumen.

**Composition used: POPC : cholesterol at 7:3 (30 mol%).**

### What this membrane is not

- **Not the enterocyte apical membrane.** At ~50 mol% cholesterol with heavy
  glycosphingolipid content, that is a different physical environment.
- **No leaflet asymmetry.** Real plasma membranes keep PS and PE inner-facing.
- **No PE, PS, sphingomyelin or glycosphingolipid.** POPC stands in for the
  entire phospholipid fraction. This follows the build instructions, which
  specify "POPC with cholesterol"; a multi-component bilayer would substitute a
  specified component and is not done silently.
- **No glycocalyx.**
- **The 30 mol% figure is measured for intestinal epithelial basolateral
  membrane, not for enteric neurons.** No enteric-neuron lipidomics was found.
  It is plausible for neuronal plasma membrane generally, but it is not a
  measured value for the tissue the benchmark affinities come from.

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

## Decisions

Every tool and parameter choice, with its reasoning, the alternatives
considered, the cost accepted, and the condition that should trigger revisiting
it: **`docs/decisions.md`**. Reversals are recorded rather than overwritten.

Open items currently recorded there: the LICENSE copyright holder, the OSU
target architecture, and two provisional choices (AM1-BCC charges, dropping
human casoxin B as a duplicate).

## Provenance

Cluster details, module versions, and build flags: `docs/provenance.md`.
Architecture and engine rationale: `docs/arch_notes.md`.

**Note:** the build instructions describe ORCA as an Oregon State University
cluster. The host in use is `login.orca.pdx.edu` (Portland State University).
The discrepancy is recorded rather than resolved.
