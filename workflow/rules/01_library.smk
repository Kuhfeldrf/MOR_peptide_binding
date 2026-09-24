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
