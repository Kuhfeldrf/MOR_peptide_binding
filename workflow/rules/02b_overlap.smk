# Stage 2b - training-overlap audit. CPU only. Emits HIGH_OVERLAP/MODERATE/NOVEL.
#
# STATUS: STUBBED. The rule structure, inputs, and outputs are real; the
# script it calls exits non-zero with a STUBBED marker until implemented.

rule overlap_audit:
    input:
        "results/01_library/peptides.tsv"
    output:
        "results/02b_overlap/overlap.tsv"
    log:
        "logs/overlap_audit.log"
    conda:
        "../../environment.yml"
    script:
        "../../scripts/02b_overlap.py"
