# Stage 2 - Chai-1 co-folding, single-sequence mode, 5 seeds per peptide.
#
# STATUS: STUBBED. The rule structure, inputs, and outputs are real; the
# script it calls exits non-zero with a STUBBED marker until implemented.

rule cofold_chai:
    input:
        "results/01_library/peptides.tsv", "results/00_receptor/mOR_clean.pdb"
    output:
        "results/02_cofold/scores.tsv"
    log:
        "logs/cofold_chai.log"
    conda:
        "../../envs/chai.yml"
    script:
        "../../scripts/02_cofold.py"
