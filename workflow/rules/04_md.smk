# Stage 4 - membrane MD with GROMACS. See docs/arch_notes.md for engine choice.
#
# STATUS: STUBBED. The rule structure, inputs, and outputs are real; the
# script it calls exits non-zero with a STUBBED marker until implemented.

rule md_production:
    input:
        "results/03_membrane/system.parm7"
    output:
        "results/04_md/stability.tsv"
    log:
        "logs/md_production.log"
    script:
        "../../scripts/04_md.py"
