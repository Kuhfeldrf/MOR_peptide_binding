#!/usr/bin/env bash
set -eo pipefail
cd /scratch/kuhfeldr-Kuhfeld_temp

cat > docs/forcefield_decision.md <<'DOC'
# Force field decision: DAMGO and the membrane system

**Decided 2026-09-24: AMBER stack — ff19SB (protein) + LIPID21 (membrane) +
GAFF2 with AM1-BCC charges (DAMGO), built with `packmol-memgen` in its default
Amber mode, run in GROMACS.**

This document records a reversal. An earlier decision the same day selected
CHARMM36m + CGenFF on the strength of literature precedent. That decision was
**superseded** after the reasoning behind it was found to be incomplete. Both
the original argument and its flaw are kept below, because a decision record
that hides its own error is not evidence of anything.

## The problem restated

DAMGO is Tyr-D-Ala-Gly-N-MePhe-Gly-ol. Three features have no standard protein
force-field representation:

| Feature | PDB component | Why it is a problem |
|---------|---------------|---------------------|
| D-alanine | `DAL` | Mirror-image backbone; L-residue parameters do not transfer directly |
| N-methyl-phenylalanine | `MEA` | N-methylation removes the backbone amide H and alters omega/phi propensities |
| Glycinol (C-terminal alcohol) | `ETA` | Replaces the C-terminal carboxylate; removes a formal negative charge |

**Structure is not the problem.** DAMGO's bound coordinates come directly from
6DDF entity 5 (see `docs/damgo_notes.md`); no modelling or prediction is
involved. What is missing is the **energy function** - partial charges,
Lennard-Jones terms, bonded parameters. A structure records where atoms sat in
one experiment; it says nothing about the forces acting on them, and parameters
cannot be derived from coordinates.

## The reversal, and why

The superseded decision argued that the Amber route required "custom
parameterisation with no published precedent," because there is no
ff14SB/ff19SB-compatible parameter library for N-methylated amino acids. **That
is true only of a residue-based treatment** - DAMGO as ff19SB residues plus two
custom ones.

The same decision then proposed treating DAMGO as a **single whole molecule**
under CGenFF, to avoid the mixed-force-field junction. That reframing was never
applied back to the Amber side. It should have been: **`antechamber` + GAFF2
parameterises arbitrary organic molecules, N-methylated amides and primary
alcohols included. That is what GAFF exists for.** Treating DAMGO as one GAFF2
molecule sidesteps precisely the problem cited as disqualifying, and needs no
custom residue library at all.

Once that is corrected, the comparison is close, and the deciding factor is not
precedent:

| | GAFF2 + AM1-BCC | CGenFF |
|---|---|---|
| Licence | **None** (AmberTools, open source) | Academic registration; binary needs an institutional signatory, ~weeks |
| Scriptable | **Yes - already installed and verified** | Web app per molecule, or the licensed binary |
| Precedent for DAMGO-muOR | Not this specific system | Yes, the published studies |
| System force field | ff19SB + LIPID21 | CHARMM36m |
| ABFE track record | Very heavily used and validated | Also standard |

**Chosen on the grounds of what can actually be run and reproduced.** For a
repository whose entire purpose is demonstrating reproducibility, a
licence-gated interactive web step in the middle of parameter generation is a
worse defect than using a well-validated open force field with less direct
precedent for this exact receptor.

A hard constraint drove the rest: **the ligand force field must match the
system force field.** GAFF2 with CHARMM36m is invalid - different
Lennard-Jones combination rules. Choosing the ligand force field therefore
chooses the whole stack.

## The decision

| Component | Choice |
|-----------|--------|
| Protein | ff19SB |
| Membrane | LIPID21 (POPC + cholesterol) |
| Water / ions | OPC (ff19SB's matched water model), 150 mM |
| DAMGO | **GAFF2, as a single molecule**, AM1-BCC charges via `antechamber`/`sqm` |
| Builder | `packmol-memgen`, default Amber output - **unchanged from the instructions** |
| Engine | GROMACS, via ParmEd topology conversion |

DAMGO is treated as one GAFF2 molecule rather than as residues. This avoids a
mixed-force-field junction mid-peptide, which is a known source of quiet error
where charges must be reconciled across the boundary.

## Costs accepted

1. **Amber-to-GROMACS topology conversion returns.** `packmol-memgen` emits
   Amber topologies; GROMACS needs them via ParmEd. This was briefly removed by
   the CHARMM decision and is now back. **Stage 3 must verify atom count and
   total system charge across the conversion** and fail loudly on a mismatch.
2. **AM1-BCC rather than RESP.** Standard and well validated for ligand work,
   and far cheaper. If Stage 6 convergence looks suspect, RESP at HF/6-31G(d,p)
   is the documented escalation - it is the level Amber's own charges were
   derived at.
3. **Less direct precedent for this exact receptor.** Published DAMGO-muOR MD
   used CHARMM36m + CGenFF. This is a real cost, accepted deliberately in
   exchange for a pipeline that runs without a licence gate.

## To revisit

- **If ABFE convergence is poor** or DAMGO's DeltaG is badly off the published
  affinity: re-derive charges with RESP before questioning the sampling.
- **If the pipeline scales to many novel peptides**: every new sequence needs
  fresh ligand parameters. `antechamber` scales fine; CGenFF's web app would
  not, which is a further argument for the route chosen here.
- **If a future comparison against the published CHARMM36m work is wanted**:
  `packmol-memgen --charmm` exists, so the CHARMM path stays open as a
  cross-check rather than being foreclosed.

## Sources

- 6DDF, DAMGO-bound mu-opioid receptor-Gi complex:
  Koehl et al., Nature (2018) - https://www.nature.com/articles/s41586-018-0219-7
- Published DAMGO-muOR MD using CHARMM36m + CGenFF (the precedent not followed):
  Cell Research (2025) - https://www.nature.com/articles/s41422-025-01191-8
- Amber-based mu-opioid ligand parameterisation with GAFF + RESP:
  Frontiers Mol Neurosci (2019) -
  https://www.frontiersin.org/journals/molecular-neuroscience/articles/10.3389/fnmol.2019.00207/full
- Amber force-field models: https://ambermd.org/AmberModels_proteins.php
- ff19SB: Tian et al., JCTC (2020) -
  https://pubs.acs.org/doi/abs/10.1021/acs.jctc.9b00591
DOC

python3 - <<'PY'
import pathlib
p = pathlib.Path("config/config.yaml"); t = p.read_text()
t = t.replace("""  # CHARMM36m + CGenFF, per docs/forcefield_decision.md. packmol-memgen is
  # unchanged as the builder; only its output format differs (--charmm).
  output_format: charmm
  forcefield:
    protein: charmm36m
    lipid: charmm36
    ligand: cgenff
    water: tip3p_charmm""",
"""  # AMBER stack, per docs/forcefield_decision.md. packmol-memgen default mode.
  output_format: amber
  forcefield:
    protein: ff19SB
    lipid: lipid21
    ligand: gaff2
    water: opc
  # Amber -> GROMACS conversion happens here. Stage 3 MUST verify atom count
  # and total charge across it and fail loudly on mismatch.
  verify_conversion: true""")
t = t.replace("""  # DAMGO is decoupled as a single CGenFF molecule, avoiding a mixed
  # force-field junction mid-peptide. See docs/forcefield_decision.md.
  ligand_forcefield: cgenff""",
"""  # DAMGO is decoupled as a single GAFF2 molecule, avoiding a mixed
  # force-field junction mid-peptide. See docs/forcefield_decision.md.
  ligand_forcefield: gaff2
  ligand_charges: am1bcc      # escalate to RESP if convergence is poor""")
p.write_text(t); print("config updated")
PY

python3 - <<'PY'
import pathlib
p = pathlib.Path("docs/damgo_notes.md"); t = p.read_text()
i = t.find("## Resolved: parameterisation")
assert i > 0
p.write_text(t[:i] + """## Resolved: parameterisation

**Settled 2026-09-24 - see `docs/forcefield_decision.md`.** DAMGO is treated as
a single **GAFF2** molecule with **AM1-BCC** charges, generated by
`antechamber`/`parmchk2` (AmberTools, already installed and verified), under
**ff19SB + LIPID21**.

Treating it as one molecule rather than as residues avoids a mixed-force-field
junction mid-peptide, and means **no N-methyl amino acid residue library is
needed** - GAFF2 parameterises arbitrary organic molecules, N-methylated amides
and primary alcohols included.

An earlier decision that day chose CHARMM36m + CGenFF on literature precedent;
it was reversed because CGenFF parameter generation is licence-gated and
interactive, which conflicts with the reproducibility this repository exists to
demonstrate. The reversal and its reasoning are recorded in
`docs/forcefield_decision.md` rather than silently overwritten.

**Nothing here blocks Stage 3.** The parameterisation runs entirely from tools
already present in `mor-pilot`.
""")
print("damgo_notes updated")
PY

python3 - <<'PY'
import pathlib
p = pathlib.Path("docs/arch_notes.md"); t = p.read_text()
t = t.replace("""- ~~**Topology conversion at the Stage 3 boundary.**~~ **No longer applies.**
  The force-field decision of 2026-09-24 (`docs/forcefield_decision.md`) moved
  the stack to CHARMM36m, which GROMACS reads natively via
  `packmol-memgen --charmm`. The Amber-to-GROMACS conversion step, previously
  flagged here as a silent-error risk, is removed entirely. This was not the
  reason for the force-field choice, but it is a real secondary benefit.""",
"""- **Topology conversion at the Stage 3 boundary.** `packmol-memgen` emits Amber
  topologies; GROMACS needs them converted via ParmEd. A routine but real extra
  step, and a genuine place for silent error. **Stage 3 must verify atom count
  and total system charge across the conversion and fail loudly on mismatch.**
  (This risk was briefly removed by a CHARMM36m decision on 2026-09-24 and
  returned when that decision was reversed - see `docs/forcefield_decision.md`.)""")
p.write_text(t); print("arch_notes updated")
PY

git add -A
git commit -q -m "Revert to AMBER stack: ff19SB + LIPID21 + GAFF2/AM1-BCC

Supersedes the CHARMM36m + CGenFF decision made earlier today. That decision
rested on an incomplete argument: it ruled out the Amber route because no
N-methyl amino acid residue library exists for ff19SB, then proposed treating
DAMGO as a single whole molecule to avoid a force-field junction - without
applying that same reframing back to Amber. antechamber and GAFF2 parameterise
arbitrary organic molecules, N-methylated amides and alcohols included, so the
cited blocker does not apply to a whole-molecule treatment.

With that corrected the comparison is close, and the deciding factor is that
CGenFF generation is licence-gated and interactive while antechamber is already
installed and verified. For a repository whose purpose is demonstrating
reproducibility, a licence gate in the middle of parameter generation is the
worse defect.

Accepts the costs: the Amber-to-GROMACS topology conversion returns as a
silent-error risk and Stage 3 must verify atom count and total charge across
it; AM1-BCC rather than RESP, with RESP documented as the escalation if Stage 6
convergence is poor; and less direct precedent for this exact receptor.

Stage 3 is no longer blocked on an external licence.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
git --no-pager log --oneline | head -3
