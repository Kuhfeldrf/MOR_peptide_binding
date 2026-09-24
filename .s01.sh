#!/usr/bin/env bash
set -eo pipefail
cd /scratch/kuhfeldr-Kuhfeld_temp

cat > workflow/rules/00_receptor.smk <<'RULE'
# Stage 0 - receptor preparation from 6DDF.
#
# STATUS: WORKING. Output verified byte-identical across runs.
#
# Also emits the experimentally resolved DAMGO ligand (chain D, entity 5),
# which Stage 6 uses as its ABFE starting structure rather than a predicted
# pose. See docs/damgo_notes.md.

rule receptor_prep:
    input:
        cif="data/raw/6ddf.cif",
    output:
        receptor="results/00_receptor/mOR_clean.pdb",
        damgo="results/00_receptor/damgo_ref.pdb",
        checksums="results/00_receptor/checksums.json",
        log="results/00_receptor/receptor_prep.log",
    params:
        ph=config["receptor"]["ph"],
    conda:
        "../../environment.yml"
    shell:
        "python3 scripts/00_receptor.py "
        "--cif {input.cif} --outdir results/00_receptor --ph {params.ph}"
RULE

cat > workflow/rules/01_library.smk <<'RULE'
# Stage 1 - peptide library ingestion.
#
# STATUS: WORKING.
#
# Joins the literature-sourced affinity data onto the screening library by
# sequence. Nothing is silently dropped: every input row appears in the output
# with an explicit `included` flag and, where excluded, a reason.

rule library_ingest:
    input:
        library="data/reference/known_opioid_peptides.csv",
        reference="data/reference/reference_peptides.tsv",
        config="config/config.yaml",
    output:
        table="results/01_library/peptides.tsv",
        fasta="results/01_library/peptides.fasta",
        stats="results/01_library/library_stats.json",
    conda:
        "../../environment.yml"
    shell:
        "python3 scripts/01_library.py "
        "--library {input.library} --reference {input.reference} "
        "--config {input.config} --outdir results/01_library"
RULE

python3 - <<'PY'
import pathlib
p=pathlib.Path("README.md"); t=p.read_text()
t=t.replace(
"""| 0 | Receptor preparation (6DDF) | `STUBBED` | - | - |
| 1 | Peptide library ingestion | `STUBBED` | - | - |""",
"""| 0 | Receptor preparation (6DDF) | `WORKING` | ~40 s | login node, CPU |
| 1 | Peptide library ingestion | `WORKING` | < 1 s | login node, CPU |""")

t=t.replace(
"""**Every stage is currently stubbed.** The rule graph, inputs, and outputs are
real and wired together; each stub script exits non-zero with a `STUBBED`
marker rather than writing an empty or fabricated output, so a stubbed stage
cannot be mistaken for one that ran.""",
"""**Stages 0 and 1 run. Stages 2-7 are stubbed.** The rule graph, inputs, and
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
0.2-200 uM. No peptide was dropped by the filters.""")
p.write_text(t); print("README updated")
PY

git add -A
git commit -q -m "Implement Stages 0 and 1

Stage 0 extracts the receptor and the experimentally resolved DAMGO from 6DDF,
assigns protonation at pH 7.4 via PROPKA, and verifies its own output. Two
instruction expectations did not hold and are recorded rather than worked
around: 6DDF contains no nanobody or scFv16, and the receptor has no chain
breaks, so ECL2 is fully resolved and no loop modelling is performed.

Protonation states are determined from the hydrogens actually placed rather
than from residue names, because pdb2pqr --pdb-output keeps the generic
HIS/ASP names and cannot distinguish HID/HIE/HIP. D147 is deprotonated (-1);
H297 is HID, neutral and delta-protonated.

Determinism verified as the instructions require: two runs from a clean
directory give byte-identical structural outputs, with SHA-256 digests written
to checksums.json. Timestamps are confined to the run log.

Stage 1 joins the literature-sourced affinities onto the screening library by
sequence, matching only unmodified C-termini so an amide or methyl ester cannot
be joined onto the wrong peptide. 40 peptides in, 13 with affinity, benchmark
subset n=10. Every input row is reported with an explicit included flag and an
exclusion reason where applicable, so nothing is dropped silently.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
git --no-pager log --oneline | head -3
