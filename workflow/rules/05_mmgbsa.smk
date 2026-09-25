# Stage 5 - MM/GBSA triage, PER PEPTIDE.
#
# STATUS: STUBBED.
#
# PURPOSE IS RANK ORDERING, NOT ABSOLUTE AFFINITY. MM/GBSA omits explicit
# solvent, treats the membrane implicitly and neglects configurational entropy;
# its numbers are not free energies. In the scaled-up design this is the cheap
# filter that decides which peptides are worth Stage 6, which costs orders of
# magnitude more per peptide.
#
# Read over the production trajectory at a stride, so the cost scales with the
# number of frames analysed rather than the length of the run.


rule mmgbsa:
    input:
        xtc=f"{RESULTS}/04_md/{{pep}}/prod.xtc",
        tpr=f"{RESULTS}/04_md/{{pep}}/prod.tpr",
        parm=f"{RESULTS}/03_membrane/{{pep}}/system.parm7",
    output:
        score=f"{RESULTS}/05_mmgbsa/{{pep}}/score.tsv",
    params:
        stride=config["mmgbsa"]["stride"],
    log:
        f"{LOGS}/05_mmgbsa_{{pep}}.log",
    resources:
        runtime=240,
        mem_mb=32000,
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/05_mmgbsa.py "
        "--traj {input.xtc} --topology {input.parm} "
        "--stride {params.stride} --out {output.score} > {log} 2>&1"


rule mmgbsa_collect:
    """One table across every peptide that reached Stage 5."""
    input:
        expand(f"{RESULTS}/05_mmgbsa/{{pep}}/score.tsv", pep=MD_PEPTIDES),
    output:
        table=f"{RESULTS}/05_mmgbsa/scores.tsv",
    log:
        f"{LOGS}/05_collect.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/05b_collect.py "
        "--dir $(dirname {output.table}) --out {output.table} > {log} 2>&1"
