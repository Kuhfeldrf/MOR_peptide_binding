# Stage 0 - receptor preparation. Output must be byte-identical across runs.
#
# STATUS: STUBBED. The rule structure, inputs, and outputs are real; the
# script it calls exits non-zero with a STUBBED marker until implemented.

rule receptor_prep:
    input:
        config["receptor"]["primary_pdb"]
    output:
        "results/00_receptor/mOR_clean.pdb"
    log:
        "logs/receptor_prep.log"
    script:
        "../../scripts/00_receptor.py"
