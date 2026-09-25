# Stage 7 - benchmark against published affinity.
#
# STATUS: STUBBED.
#
# THE SCIENTIFIC POINT OF THIS STAGE is the contrast between learned-confidence
# scores and physics-based scores, reported SEPARATELY. That contrast is the
# evidence for why the expensive downstream stages exist, so collapsing them
# into one ranking would destroy the figure's argument.
#
# Stage 2 already showed the confidence side fails here: ipTM sat flat at
# 0.27-0.35 across a 32-fold affinity range and ipSAE could not separate a
# correct pose from a failed one either (D21). Cross-seed agreement was the
# only Stage 2 diagnostic that discriminated, so it is carried in as a column
# rather than left in a log.
#
# STRATIFY BY OVERLAP FLAG. Reference peptides predict well partly through
# memorisation, so performance on NOVEL peptides is the honest measure.
#
# CAVEATS THE FIGURE MUST CARRY, not just the methods text:
#   * n is small. Report it on the plot, not only in the caption.
#   * affinities are IC50 in micromolar, not Ki; no conversion is applied.
#   * the set mixes laboratories. Koch 1985 reports bovine casomorphins an
#     order of magnitude from the Brantl 1981 values used here, same peptides
#     and same assay type, which bounds what any correlation can mean.
#   * the receptor is active-state only, while some benchmark peptides are
#     antagonists that prefer the inactive state.


rule benchmark:
    input:
        cofold=f"{RESULTS}/02_cofold/scores.tsv",
        agreement=f"{RESULTS}/02_cofold/cross_seed_agreement.tsv",
        ipsae=f"{RESULTS}/02_cofold/ipsae.tsv",
        overlap=f"{RESULTS}/02b_overlap/overlap.tsv",
        mmgbsa=f"{RESULTS}/05_mmgbsa/scores.tsv",
        reference=config["library"]["reference_set"],
    output:
        table=f"{RESULTS}/07_benchmark/benchmark.tsv",
        figure=f"{RESULTS}/07_benchmark/benchmark.svg",
        summary=f"{RESULTS}/07_benchmark/summary.md",
    params:
        stratify=config["benchmark"]["stratify_by"],
        assay=config["benchmark"]["assay_filter"],
    log:
        f"{LOGS}/07_benchmark.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/07_benchmark.py "
        "--cofold {input.cofold} --agreement {input.agreement} "
        "--ipsae {input.ipsae} --overlap {input.overlap} "
        "--mmgbsa {input.mmgbsa} --reference {input.reference} "
        "--stratify-by {params.stratify} --assay {params.assay} "
        "--table {output.table} --figure {output.figure} "
        "--summary {output.summary} > {log} 2>&1"
