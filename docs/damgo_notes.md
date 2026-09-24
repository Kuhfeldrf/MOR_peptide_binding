# DAMGO: representation decision

**Decided 2026-09-24. Resolved by extraction from experiment, not prediction.**

## The problem

The build instructions name DAMGO as the first ABFE target (Stage 6). DAMGO is
**Tyr-D-Ala-Gly-N-MePhe-Gly-ol**: a D-amino acid, an N-methylated residue, and a
C-terminal alcohol instead of a carboxylate. It therefore has **no valid
canonical one-letter FASTA representation**, so it cannot be handed to Chai-1
as a sequence the way every other peptide in the library can.

## The resolution

**DAMGO is already present, experimentally resolved, in 6DDF** - which is the
Stage 0 primary receptor structure. It is downloaded regardless. Verified
directly from `data/raw/6ddf.cif`:

```
_entity 5  polymer syn  DAMGO  513.587  'analogue of enkephalin'
chain 5 residues: TYR DAL GLY MEA ETA
```

where the PDB chemical components map exactly onto DAMGO's chemistry:

| Component | Identity | DAMGO position |
|-----------|----------|----------------|
| `TYR` | L-tyrosine | Tyr1 |
| `DAL` | **D-alanine** | D-Ala2 |
| `GLY` | glycine | Gly3 |
| `MEA` | **N-methyl-phenylalanine** | N-MePhe4 |
| `ETA` | **ethanolamine** | Gly5-ol (C-terminal alcohol) |

So the correct move is to **extract DAMGO's coordinates from 6DDF** rather than
predict them.

## Why this is better, not merely easier

1. **ABFE should start from the experimental pose.** Stage 6 measures a binding
   free energy; feeding it a predicted pose would compound pose error into the
   free energy. The cryo-EM pose is the best available starting structure.
2. **Receptor and ligand come from the same structure.** Stage 0 builds the
   receptor from 6DDF, so the extracted DAMGO is already correctly positioned in
   that receptor's frame. No docking or alignment step is required, and no
   opportunity for silent misplacement is introduced.
3. **It removes the non-canonical input problem from Stage 2 entirely** instead
   of working around it.

## What this decision costs - state this in the README

**DAMGO's ABFE result does not test pose prediction.** It starts from the
experimental complex, so it validates the free-energy machinery
(decoupling, restraints, MBAR, convergence) and nothing upstream of it. That is
the appropriate first validation target - it isolates one thing - but it must
not be presented as end-to-end validation of the pipeline.

## Consequences per stage

- **Stage 0** - additionally emit `results/00_receptor/damgo_ref.pdb` from 6DDF
  entity 5, alongside the cleaned receptor. Same determinism requirement.
- **Stage 2 (Chai-1)** - DAMGO is **excluded from sequence-based co-folding**.
  It is a validation standard, not a screening candidate.
  *Optional diagnostic:* supply DAMGO to Chai-1 as a CCD/SMILES ligand and ask
  whether co-folding recovers the known pose. Since 6DDF long predates any
  current model's training cutoff this is a **memorisation check, not a
  generalisation test**, and Stage 2b must flag it `HIGH_OVERLAP`. Informative,
  but it must be labelled honestly.
- **Stage 3/4 (parameters)** - `TYR` and `GLY` are standard. `DAL` is D-alanine,
  obtainable by mirroring the L-Ala parameters in ff14SB/ff19SB. **`MEA` and
  `ETA` have no standard protein parameters** and require GAFF2 via
  `antechamber`/`parmchk2` (AmberTools, already in `environment.yml`), with
  charges derived consistently with the rest of the peptide.
  **This is the remaining technical risk on the DAMGO path** and should be
  verified early, not discovered during Stage 6.
- **Stage 7** - DAMGO's published affinity must still be sourced to primary
  literature under the Stage 1 rules; it is not exempt.

## Resolved: parameterisation

**Settled 2026-09-24 - see `docs/forcefield_decision.md`.** DAMGO is treated as
a single **CGenFF** molecule under **CHARMM36m**, matching published DAMGO-mu-OR
MD, which uses exactly this combination in GROMACS. This removes the
mixed-force-field junction that a per-residue Amber treatment would have
created mid-peptide.

The RESP-versus-AM1-BCC question does not arise: CGenFF assigns charges by
analogy rather than using either scheme.

**One item remains open:** CGenFF parameter generation needs the licensed
`cgenff` binary or the interactive ParamChem web service, and the latter would
break the pipeline's scriptability requirement. Stage 3 stays `STUBBED` until
that is resolved. See the FLAG section of `docs/forcefield_decision.md`.
