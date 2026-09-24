# Stage 6 - ABFE by alchemical decoupling, MBAR. Energies + subsampled coords only.
#
# STATUS: STUBBED. The rule structure, inputs, and outputs are real; the
# script it calls exits non-zero with a STUBBED marker until implemented.

rule abfe_windows:
    input:
        "results/04_md/stability.tsv"
    output:
        "results/06_abfe/dg.tsv"
    log:
        "logs/abfe_windows.log"
    script:
        "../../scripts/06_abfe.py"
