# Stage 5 - MM/GBSA triage. Rank ordering only, not absolute affinity.
#
# STATUS: STUBBED. The rule structure, inputs, and outputs are real; the
# script it calls exits non-zero with a STUBBED marker until implemented.

rule mmgbsa_triage:
    input:
        "results/04_md/stability.tsv"
    output:
        "results/05_mmgbsa/scores.tsv"
    log:
        "logs/mmgbsa_triage.log"
    script:
        "../../scripts/05_mmgbsa.py"
