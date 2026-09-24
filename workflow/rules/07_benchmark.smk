# Stage 7 - Spearman + MAE against published affinity, stratified by overlap flag.
#
# STATUS: STUBBED. The rule structure, inputs, and outputs are real; the
# script it calls exits non-zero with a STUBBED marker until implemented.

rule benchmark:
    input:
        "results/02_cofold/scores.tsv", "results/05_mmgbsa/scores.tsv", "results/02b_overlap/overlap.tsv"
    output:
        "results/07_benchmark/benchmark.tsv"
    log:
        "logs/benchmark.log"
    script:
        "../../scripts/07_benchmark.py"
