#!/usr/bin/env bash
set -eo pipefail
cd /scratch/kuhfeldr-Kuhfeld_temp

cat > docs/forcefield_decision.md <<'DOC'
# Force field decision: DAMGO and the membrane system

**Decided 2026-09-24 after a literature review. Chosen: CHARMM36m (protein,
lipids) + CGenFF (DAMGO), in GROMACS.**

This supersedes the earlier working assumption of ff19SB + LIPID21 + GAFF2.

## The problem restated

DAMGO is Tyr-D-Ala-Gly-N-MePhe-Gly-ol. Three features have no standard protein
force-field representation:

| Feature | PDB component | Why it is a problem |
|---------|---------------|---------------------|
| D-alanine | `DAL` | Mirror-image backbone; L-residue parameters do not transfer directly |
| N-methyl-phenylalanine | `MEA` | N-methylation removes the backbone amide H and alters omega/phi propensities |
| Glycinol (C-terminal alcohol) | `ETA` | Replaces the C-terminal carboxylate; removes a formal negative charge |

Stage 6 computes an absolute binding free energy for DAMGO. **Ligand parameter
error propagates directly into that DeltaG and does not announce itself.**

## What the literature actually does

Published MD of the DAMGO-bound mu-opioid receptor uses **CHARMM36m for the
protein, lipids and ions, with CGenFF for the ligand, run in GROMACS**. That is
the same engine this pilot selected, so the precedent transfers without
reinterpretation.

On the AMBER side the picture is weaker. There is **no published,
ff14SB/ff19SB-compatible parameter library for N-methylated amino acids**.
Amber-family coverage of protein mimetics exists in `ff15ipq-m` (D-alpha and
C-alpha-methylated residues), but that is a different force-field family and
would mean abandoning ff19SB for the receptor. Published Amber-based mu-opioid
work parameterises *small-molecule* ligands (morphine, fentanyl, herkinorin)
with GAFF/GAFF2 and RESP charges - not a modified peptide like DAMGO.

So the Amber route for DAMGO would require **custom parameterisation with no
published precedent to check it against**, on the one molecule whose free
energy is the pilot's validation target.

## Why CGenFF suits this molecule specifically

DAMGO is 513 Da - small-molecule sized. Treating it as **a single CGenFF
molecule** avoids the mixed-force-field junction entirely, rather than joining
ff19SB residues to GAFF2 residues mid-chain and having to reconcile charges
across the boundary. That junction is a well-known source of quiet error, and
removing it is worth more here than force-field purity.

CGenFF and CHARMM36m are designed to be combined: shared Lennard-Jones
combination rules and a common parameterisation philosophy.

CGenFF also emits a **penalty score** per assigned parameter. That is an
auditable, quantitative statement of how far each parameter was extrapolated by
analogy - exactly the kind of honest uncertainty signal this pilot is supposed
to surface rather than hide. Convention: penalties **< 10** are acceptable,
**10-50** warrant inspection, **> 50** require QM refinement. The penalties
must be recorded in `results/` and reported, not just consulted and discarded.

## Effect on the rest of the pipeline

| Stage | Change |
|-------|--------|
| 3 | `packmol-memgen --charmm` emits CHARMM-format output. **The mandated Stage 3 tool is unchanged** - only its output format. |
| 3 | Lipids become CHARMM36 rather than LIPID21. POPC and cholesterol are both well validated in C36. |
| 4 | GROMACS ships CHARMM36m natively; no conversion step. **This removes the Amber-to-GROMACS topology conversion previously flagged as a silent-error risk** in arch_notes. |
| 6 | DAMGO decoupled as one CGenFF molecule - simpler alchemical setup than a multi-residue mixed-FF ligand. |

Water: CHARMM-modified TIP3P, per C36 convention.

## FLAG: CGenFF parameter generation is not yet scriptable here

Generating CGenFF parameters requires either the **`cgenff` binary** (free for
academic use, but requires registration) or the **ParamChem web service**.
The web service is interactive and **would break the "scripts cleanly and
reproducibly" requirement** that motivated choosing packmol-memgen over
CHARMM-GUI in the first place.

**This needs resolving before Stage 3 runs.** Options, in preference order:

1. Obtain the `cgenff` binary under an academic licence and vendor it into the
   container. Keeps the pipeline scriptable. **Preferred.**
2. Generate DAMGO's stream file once via ParamChem, commit it to
   `data/reference/damgo_cgenff.str`, and document it as a **manual,
   non-reproducible step** with the date and input hash recorded. Honest but
   weaker, and it must be labelled as such in the README rather than glossed.
3. Fall back to the Amber route (ff19SB + GAFF2/RESP for `MEA`/`ETA`). Keeps
   everything scriptable but discards the literature precedent and takes on the
   junction problem.

Until this is resolved, Stage 3 stays `STUBBED`. **No substitute force field
will be selected silently.**

## Charge derivation, resolved

The earlier open question - RESP vs AM1-BCC - **does not arise under this
decision**. CGenFF assigns charges by analogy from its own parameterised
training set; it does not use either scheme. If option 3 above is taken
instead, the question returns, and the answer would then be RESP at
HF/6-31G(d,p), which is the level Amber's own charges were derived at and the
level used in published Amber-based mu-opioid ligand work.

## Sources

- CHARMM36m + CGenFF for DAMGO-mu-opioid receptor MD in GROMACS:
  Cell Research (2025), "The molecular basis of mu-opioid receptor signaling
  plasticity" - https://www.nature.com/articles/s41422-025-01191-8
- 6DDF, DAMGO-bound mu-opioid receptor-Gi complex:
  Koehl et al., Nature (2018) - https://www.nature.com/articles/s41586-018-0219-7
- Amber-based mu-opioid ligand parameterisation (GAFF + RESP, small molecules):
  Frontiers Mol Neurosci (2019) -
  https://www.frontiersin.org/journals/molecular-neuroscience/articles/10.3389/fnmol.2019.00207/full
- Amber protein force-field coverage, incl. ff15ipq-m for protein mimetics:
  https://ambermd.org/AmberModels_proteins.php
- ff19SB: Tian et al., JCTC (2020) -
  https://pubs.acs.org/doi/abs/10.1021/acs.jctc.9b00591
DOC

python3 - <<'PY'
import pathlib
p = pathlib.Path("config/config.yaml"); t = p.read_text()
t = t.replace("""membrane:
  builder: packmol-memgen""","""membrane:
  builder: packmol-memgen
  # CHARMM36m + CGenFF, per docs/forcefield_decision.md. packmol-memgen is
  # unchanged as the builder; only its output format differs (--charmm).
  output_format: charmm
  forcefield:
    protein: charmm36m
    lipid: charmm36
    ligand: cgenff
    water: tip3p_charmm""")
t = t.replace("""  estimator: MBAR""","""  estimator: MBAR
  # DAMGO is decoupled as a single CGenFF molecule, avoiding a mixed
  # force-field junction mid-peptide. See docs/forcefield_decision.md.
  ligand_forcefield: cgenff""")
p.write_text(t); print("config updated")
PY

python3 - <<'PY'
import pathlib
p = pathlib.Path("docs/damgo_notes.md"); t = p.read_text()
t = t.replace("""## Still open

The **charge-derivation scheme for `MEA` and `ETA`** is not yet decided
(RESP vs AM1-BCC, and whether to parameterise DAMGO as one GAFF2 unit or as
standard residues plus two custom ones). Flagged rather than assumed.""",
"""## Resolved: parameterisation

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
that is resolved. See the FLAG section of `docs/forcefield_decision.md`.""")
p.write_text(t); print("damgo_notes updated")
PY

python3 - <<'PY'
import pathlib
p = pathlib.Path("docs/arch_notes.md"); t = p.read_text()
t = t.replace("""- **Topology conversion at the Stage 3 boundary.** `packmol-memgen` emits Amber
  topologies; GROMACS needs them converted via ParmEd or acpype. This is a
  routine but real extra step, and a genuine place for silent error. Stage 3
  must verify atom counts and total charge across the conversion.""",
"""- ~~**Topology conversion at the Stage 3 boundary.**~~ **No longer applies.**
  The force-field decision of 2026-09-24 (`docs/forcefield_decision.md`) moved
  the stack to CHARMM36m, which GROMACS reads natively via
  `packmol-memgen --charmm`. The Amber-to-GROMACS conversion step, previously
  flagged here as a silent-error risk, is removed entirely. This was not the
  reason for the force-field choice, but it is a real secondary benefit.""")
p.write_text(t); print("arch_notes updated")
PY

git add -A
git commit -q -m "Choose CHARMM36m + CGenFF for DAMGO, per literature precedent

Published MD of the DAMGO-bound mu-opioid receptor uses CHARMM36m with CGenFF
for the ligand, run in GROMACS - the same engine this pilot selected, so the
precedent transfers directly. The Amber alternative has no published
ff14SB/ff19SB-compatible parameter library for N-methylated amino acids, which
would have meant custom parameterisation with nothing to check it against, on
the one molecule whose free energy is the validation target.

DAMGO is 513 Da, so treating it as a single CGenFF molecule also removes the
mixed-force-field junction that a per-residue Amber treatment would create
mid-peptide. CGenFF penalty scores give an auditable measure of how far each
parameter was extrapolated, which will be recorded rather than discarded.

packmol-memgen is retained as the Stage 3 builder as mandated; only its output
format changes, via --charmm. A secondary benefit is that the Amber-to-GROMACS
topology conversion previously flagged as a silent-error risk disappears.

Flags that CGenFF generation needs a licensed binary or the interactive
ParamChem service, the latter breaking scriptability. Stage 3 stays STUBBED
until that is resolved; no substitute force field is selected silently.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
git --no-pager log --oneline | head -3
