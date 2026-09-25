# Stage 2b - training-overlap audit. CPU only, seconds to run.
#
# STATUS: WORKING.
#
# Co-folding accuracy depends strongly on how much a target resembles the
# model's training data, so reference peptides predict well partly through
# memorisation. Without this column, success on the reference set would be
# mistaken for evidence of generalisation.
#
# Stage 2 made that concrete: the one peptide whose seeds converged was
# Met-enkephalin, whose sequence is the N-terminal message sequence of
# beta-endorphin in 8F7Q. Performance on NOVEL peptides is the honest measure,
# which is why Stage 7 stratifies on this flag.


rule overlap_audit:
    input:
        peptides=f"{RESULTS}/01_library/peptides.tsv",
        config="config/config.yaml",
    output:
        table=f"{RESULTS}/02b_overlap/overlap.tsv",
        stats=f"{RESULTS}/02b_overlap/overlap_stats.json",
    log:
        f"{LOGS}/02b_overlap.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/02b_overlap.py "
        "--peptides {input.peptides} --config {input.config} "
        "--outdir $(dirname {output.table}) > {log} 2>&1"
