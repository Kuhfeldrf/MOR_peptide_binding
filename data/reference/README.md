# Curated reference peptide set

`reference_peptides.tsv` carries the published Ki values that make Stage 7 a
validation rather than a plumbing check.

**Rule: do not invent affinity values.** Every non-null `ki_nM` must carry a
`citation` resolving to a specific paper. Where a value cannot be sourced it is
left null and `ki_status` is set to `UNSOURCED`.

| Column       | Notes                                                       |
|--------------|-------------------------------------------------------------|
| `name`       | Common peptide name.                                        |
| `sequence`   | One-letter code, or `NON_CANONICAL` (see below).             |
| `ki_nM`      | Published Ki at uOR, nanomolar. Null if unsourced.           |
| `ki_status`  | `SOURCED` / `UNSOURCED`.                                     |
| `citation`   | DOI or full reference.                                       |
| `role`       | `agonist` / `negative_control`.                              |
| `notes`      | Assay conditions, species, caveats.                          |

## Open item: non-canonical residues

DAMGO is Tyr-D-Ala-Gly-N-MePhe-Gly-ol. It contains a D-amino acid, an
N-methylated residue, and a C-terminal alcohol. It therefore has **no valid
canonical one-letter FASTA representation**, which has consequences for
Stage 2 (Chai-1 sequence input) and Stage 6 (force-field parameters).
This is flagged, not silently resolved. See README "Known limitations".

## Source data

`known_opioid_peptides.csv` and `background_peptides.csv` are copied **verbatim**
from `~/bionemo_mor/data/peptides/` on ORCA. Provenance for those files is in
`peptides_source.yaml` and `PRIOR_ASSUMPTIONS.md` / `PRIOR_WORK.md`.

Contents: 20 real food-derived opioid peptides (beta-casomorphins, casoxins,
lactorphins, soymorphins, alphas1-casein exorphins, Met-enkephalin) and 20
**scrambled controls**. The scrambles serve as the negative controls the build
instructions require, and are stronger than arbitrary non-binders because their
amino-acid composition is matched to their parent sequence.

## Critical: IC50 is not Ki

The source file carries **`ic50_um` (micromolar IC50)**, not Ki. The build
instructions ask for Ki. These are **not interchangeable** and are kept in
separate columns. No Cheng-Prusoff conversion is applied: it requires the
radioligand Kd and concentration, which the source does not record.

Only 6 of 20 peptides carry any affinity value at all:

| Peptide | IC50 (uM) |
|---------|-----------|
| Met-enkephalin | 0.2 |
| beta-casomorphin-5 (bovine) | 6.5 |
| beta-casomorphin-5 (human) | 14 |
| beta-casomorphin-7 (human) | 25 |
| beta-casomorphin-7 (bovine) | 57 |
| neocasomorphin-6 | 59 |

Stage 7 reports **n** alongside every correlation. With n = 6 a Spearman
coefficient is weakly determined; that caveat belongs in the figure caption,
not just the methods.

## Citation status: UNSOURCED pending literature check

The source attributes these values to `prior_analysis:anganost.py` - **a script,
not a paper**. A script reference does not satisfy the sourcing requirement, so
every value above is currently marked `UNSOURCED` and must be traced back to
primary literature before Stage 7 treats it as a validation target.

A sourced value records the radioligand, the tissue or cell system, and the
species alongside the number. Affinities for these peptides vary by roughly an
order of magnitude across assay formats, so a bare number without conditions is
not a sourced value.

## Status

Real sequences: **present**. Affinity values: **present but UNSOURCED**.
Ki column: **empty**, pending the literature check in step 4.
