# Stage 5 - MM/GBSA triage, PER PEPTIDE.
#
# STATUS: IMPLEMENTED, NOT YET VALIDATED. The score is produced; whether it
# separates binders from decoys is the open question this stage exists to
# answer, and it is not answered until the AUROC is in.
#
# PURPOSE IS RANK ORDERING, NOT ABSOLUTE AFFINITY. MM/GBSA omits explicit
# solvent, treats the membrane implicitly and neglects configurational entropy;
# its numbers are not free energies.
#
# THIS RUNS WITHOUT MD, AND THAT IS THE POINT.
# The rule used to take prod.xtc, a 50 ns production trajectory. An ensemble
# average is the better calculation, but it makes this post-hoc analysis of the
# expensive step rather than the filter that decides whether to spend it - and
# it could never be run over the 20 scrambled decoys, so the score could never
# be validated. Scoring a minimised complex costs minutes, so the whole library
# is rankable and the decoys are included.
#
# The trajectory-ensemble variant belongs downstream, on survivors, and is not
# implemented.


rule mmgbsa:
    """Interaction energy for one peptide, from a minimised complex."""
    input:
        complex=f"{RESULTS}/03_membrane/{{pep}}/complex_opm.pdb",
        mol2=f"{RESULTS}/03_membrane/{{pep}}/params/ligand.mol2",
        frcmod=f"{RESULTS}/03_membrane/{{pep}}/params/ligand.frcmod",
        ligpdb=f"{RESULTS}/03_membrane/{{pep}}/params/ligand_unique.pdb",
    output:
        score=f"{RESULTS}/05_mmgbsa/{{pep}}/score.tsv",
    params:
        igb=config["mmgbsa"]["igb"],
        salt=config["membrane"]["ionic_strength_mM"] / 1000.0,
        minsteps=config["mmgbsa"]["minsteps"],
        ss=config["mmgbsa"]["disulfide"],
    log:
        f"{LOGS}/05_mmgbsa_{{pep}}.log",
    resources:
        runtime=180,
        mem_mb=16000,
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/05_mmgbsa.py "
        "--complex {input.complex} --mol2 {input.mol2} "
        "--frcmod {input.frcmod} --ligand-pdb {input.ligpdb} "
        "--outdir $(dirname {output.score}) --out {output.score} "
        "--peptide {wildcards.pep} --igb {params.igb} --salt {params.salt} "
        "--minsteps {params.minsteps} --disulfide {params.ss} > {log} 2>&1"


rule mmgbsa_collect:
    """One table across every peptide scored, decoys included.

    Over TRIAGE_PEPTIDES rather than MD_PEPTIDES: a triage score that is only
    computed for the peptides already selected for MD cannot be said to have
    selected anything, and cannot be tested against the decoys.
    """
    input:
        expand(f"{RESULTS}/05_mmgbsa/{{pep}}/score.tsv", pep=TRIAGE_PEPTIDES),
    output:
        table=f"{RESULTS}/05_mmgbsa/scores.tsv",
    log:
        f"{LOGS}/05_collect.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/05b_collect.py "
        "--dir $(dirname {output.table}) --out {output.table} > {log} 2>&1"
