# The positive control fails: co-folding does not recognise the opioid message

This is the strongest negative result in the pilot, and it is stronger than the
scrambled-decoy result because it does not depend on assuming the decoys are
non-binders.

## Why this test was built

The scrambled decoys control composition, which makes them a good control for
size bias and a weak one for binding. The opioid message is an N-terminal
tyrosine with a free protonated alpha-amine - the amine makes the Asp147 salt
bridge seen in 100% of frames of the 50 ns DAMGO run. Scrambling `YPFP` to
`PYPF` keeps the tyrosine, the free amine, both aromatics and the identical
hydrophobic surface. The scrambles may therefore be **weak binders rather than
non-binders**, which makes an AUROC near 0.5 ambiguous.

The fix is a control set with a **real positive**: if a peptide that genuinely
carries the pharmacophore also fails to score well, the ambiguity disappears.

## Design

| arm | n | content |
|---|---|---|
| **Positive** | 3 | beta-endorphin(1-10), (1-15), (1-20) - begins `YGGFM` |
| **Negative** | 18 | real casein / beta-lactoglobulin / alpha-lactalbumin fragments, 6 per length, **none beginning with Tyr** |

Length is varied in both arms and the message in only one, so "too long to
bind" and "lacks the pharmacophore" are separable. All sequences are traceable
to a UniProt accession and residue range; beta-endorphin was located at POMC
residue 237, matching the known 237-267 span.

## Result

**All 105 co-folding jobs completed.** There is no computational failure at 10,
15 or 20 residues - the pipeline handles the length. What it does not do is
discriminate.

| length | positive ipTM | negatives median | negatives scoring **above** the positive |
|---|---|---|---|
| 10 aa | 0.349 | 0.308 | 2 / 6 |
| 15 aa | 0.331 | 0.280 | 1 / 6 |
| **20 aa** | **0.215** | 0.276 | **5 / 6** |

**At 20 residues the known binder ranks 6th of 7**, below five random
non-opioid milk fragments.

The positive's score also degrades monotonically with length - 0.349 -> 0.331
-> 0.215 - even though full-length beta-endorphin (31 aa) is a high-affinity
muOR ligand and every one of these fragments contains met-enkephalin at its
N-terminus.

Overall AUROC is 0.574 (p = 0.74), but with n = 3 positives that number is
descriptive only. The per-length ranking above is the informative part.

### Everything diverges

All 21 controls came back `DIVERGENT`, cross-seed peptide RMSD 2.7-23.6 A, with
beta-endorphin(1-20) at 11.6 A. The model does not know where these peptides
go. That is at least honest uncertainty rather than confident error - but it
also means the poses fed to any downstream scoring are unreliable.

### ipTM does not simply track length here

Spearman rho = -0.308 (p = 0.17) across all controls, -0.223 (p = 0.37) within
the negatives. The size-bias mechanism measured for MM/GBSA
([`mmgbsa_validation.md`](mmgbsa_validation.md), vdW vs atom count r = -0.930)
does **not** extend to ipTM in this length regime. The failure here is not a
size artefact; it is a failure to see the pharmacophore.

## What this establishes

**Co-folding confidence does not recognise the opioid message.** A real peptide
carrying the real pharmacophore scores below random milk fragments. That cannot
be explained by decoys secretly binding, which was the loophole in the scramble
result.

Combined with the earlier tests, three independent scoring approaches now fail
on the same control problem:

| tier | test | result |
|---|---|---|
| ipTM, sequence | 17 real vs 20 scrambles | AUROC 0.376 |
| ipTM, 6DDF template | same | AUROC 0.415 |
| MM/GBSA | 16 real vs 18 scrambles | AUROC 0.483 |
| **ipTM, real positive** | **beta-endorphin vs 18 milk fragments** | **positive ranks 6/7 at 20 aa** |

The common factor is that **co-folding places every peptide in the orthosteric
pocket**, so every downstream score ranks structures that are all already
"bound". The information identifying a non-binder - that it fails to adopt a
favourable pose at all - is discarded before scoring begins.

**The bottleneck is pose generation, not scoring.**

## A baseline that has to be beaten

Nearly every opioid peptide begins with tyrosine. On the real library, "does
the sequence start with Tyr?" is a free, one-line predictor. **None of the
three scoring tiers tested here would beat it.** Any future method must be
measured against that baseline, or it is not adding anything.

(On this control set the rule separates the arms perfectly by construction, so
it is circular here - it is a baseline for the real library, not for these
controls.)

## Caveat on the positives

The beta-endorphin fragments are positives by **pharmacophore reasoning, not by
a cited measurement for these specific constructs**. Fragments 1-9, 1-16 and
1-17 are reported active, and 1-5 IS met-enkephalin, so substantial affinity is
a strong inference - but under this project's own sourcing rule that is an
inference, not a source. Obtaining measured affinities for the exact fragments
would make this a validated positive control rather than a reasoned one.

## Force-field dependence of the physics score

Scoring met-enkephalin's MM/GBSA interaction energy under both parameterisation
paths, same pose:

| parameterisation | DELTA TOTAL |
|---|---|
| GAFF2 / AM1-BCC | **-40.71** kcal/mol |
| ff19SB | **-26.44** kcal/mol |

A 14 kcal/mol spread for the same peptide in the same pose. This is a further
reason to treat these as ranking scores within one parameterisation and never
as affinities - and a reason to keep the force field fixed across any
comparison.

## Data

| file | contents |
|---|---|
| `../data/reference/length_controls.csv` | the 21 control sequences with provenance |
| `lengthctrl_cofold_scores.tsv` | per-seed, per-model confidence metrics |
| `lengthctrl_cross_seed.tsv` | pose reproducibility |

Generated by `scripts/01b_length_controls.py` (seeded), co-folded by
`jobs/controls.sbatch`, analysed by `scripts/07b_length_controls_analysis.py`.
