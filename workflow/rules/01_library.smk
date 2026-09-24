# Stage 1 - peptide library ingestion and filtering.
#
# STATUS: STUBBED. The rule structure, inputs, and outputs are real; the
# script it calls exits non-zero with a STUBBED marker until implemented.

rule library_ingest:
    input:
        "data/reference/known_opioid_peptides.csv"
    output:
        "results/01_library/peptides.tsv"
    log:
        "logs/library_ingest.log"
    conda:
        "../../environment.yml"
    script:
        "../../scripts/01_library.py"
