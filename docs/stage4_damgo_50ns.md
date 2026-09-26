# Stage 4 result: DAMGO, 50 ns

First completed production run. DAMGO in the mu-opioid receptor, POPC with
30 mol% cholesterol, 310 K, semi-isotropic NPT, starting from the experimental
6DDF pose.

| Quantity | Value |
|----------|-------|
| Duration | 50 ns, 501 frames |
| Wall time | 8 h 43 m on one A30 |
| Performance | 137.7 ns/day |
| Ligand RMSD, mean | **1.30 A** |
| Ligand RMSD, second half | **1.26 A** |
| Ligand centre-of-mass displacement, final | **0.50 A** |
| Receptor CA RMSD, final | 1.58 A |
| Residues in contact >50% of frames | 19 |
| Remains in pocket | yes |

## The complex is stable

A ligand RMSD of 1.3 A held over 50 ns, with the second half no higher than the
mean, is a ligand that is not going anywhere. The centre of mass finishes
0.50 A from where it started. Receptor CA RMSD of 1.58 A is unremarkable for a
GPCR on this timescale.

## The contacts are the canonical pocket

| Residue | % of frames |
|---------|-------------|
| **ASP147** | **100.0** |
| TRP318 | 98.4 |
| TYR148 | 98.2 |
| THR218 | 97.8 |
| LEU219 | 96.8 |
| ILE296 | 95.4 |
| ILE322 | 95.0 |
| VAL236 | 94.6 |
| **CYX217** | 93.8 |
| **HID297** | 93.6 |
| ILE144 | 89.8 |
| PHE221 | 89.4 |
| MET151 | 87.6 |
| VAL300 | 84.0 |

**ASP147 is in contact in every one of 501 frames.** That is the conserved TM3
aspartate whose salt bridge to the ligand's protonated amine is the defining
interaction of opioid ligand binding. The remaining residues are the
orthosteric pocket lining across TM3, TM5, TM6, TM7 and ECL2.

## Three earlier decisions are validated by this

The contact list is not just chemically sensible, it exercises three choices
that were flagged at the time as places where a silent error would have been
invisible:

1. **ASP147 deprotonated.** Stage 0 determined this from the hydrogens PROPKA
   placed, not from the residue name. Had it been protonated, the salt bridge
   could not form. It is present in 100% of frames.
2. **DAMGO's N-terminus protonated, net charge +1.** obabel silently ignores pH
   when its protonate and add-hydrogens flags are combined and produced a
   neutral amine; the assertion added afterwards caught it. A neutral amine has
   nothing to donate to ASP147.
3. **HID297 and the CYX217 disulfide.** Both appear at >93%. tleap infers
   histidine protonation from the residue NAME, which pdb2pqr does not write,
   and forms a disulfide only when told to. Both were set explicitly in
   Stage 3.

Each of those would have produced a system that ran. This trajectory is the
evidence that they were set correctly.

## What this does NOT show

**It does not show the pipeline can predict binding.** The run starts from the
experimental pose. It demonstrates that the force field, membrane, protonation
and parameters hold a known-correct structure in place over 50 ns - which is a
prerequisite for trusting a free energy, not a result about affinity.

**It is one peptide and one run.** No replicate, no statistics. A single
trajectory cannot distinguish a stable complex from one that had not yet found
its way out.

**Stability is not affinity.** A tightly held pose says nothing about how
tightly the peptide binds relative to another peptide. That comparison is
Stage 5 and Stage 6, which are still stubbed.

## Trajectory policy held

| Output | Size | Frames |
|--------|------|--------|
| Reduced group (receptor + DAMGO + proximal lipids) | 27.2 MB | 500 |
| Full system, sparse | 12.8 MB | 10 |

40 MB total for 50 ns. The instructions name uncontrolled trajectory writing as
the single most likely way to fill the filesystem; writing the full system at
analysis frequency instead would have produced roughly 2 GB for the same run.

## Correction, 2026-09-26

This report originally said the run took 8 h 43 m **on one L40S**. It was an
**A30** (job 182895, node orcaga20). The `md_production` rule was pinned to the
`long` partition, which contains only A30 nodes, so the card was a consequence
of that setting rather than a choice - and 137.7 ns/day is therefore an A30
number.

Production has since moved to `normal`, which holds the L40S nodes and still
allows 24 h. Whether the L40S is materially faster here is an expectation, not
yet a measurement.
