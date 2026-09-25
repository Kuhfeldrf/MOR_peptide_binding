# Stage 6 - absolute binding free energy, PER PEPTIDE PER LAMBDA WINDOW.
#
# STATUS: STUBBED.
#
# The pilot runs DAMGO first: its affinity is well established and its pose is
# experimental rather than predicted, so a disagreement points at the method
# rather than at the starting structure.
#
# STORAGE POLICY, set before anything runs. Full trajectories across every
# lambda window would be 30-50 GB per peptide; energies and subsampled
# coordinates bring that to 2-5 GB, and MBAR needs only the latter. The window
# rule therefore writes dhdl and a strided trajectory, never full frames.
#
# The lambda windows fan out as wildcards, so a 36-window calculation is 36
# independent jobs rather than one long serial chain.

LAMBDAS = list(range(config["abfe"]["n_lambda_windows"]))


rule abfe_window:
    """One lambda window. Independent of every other window."""
    input:
        gro=f"{RESULTS}/04_md/{{pep}}/npt_final.gro",
        top=f"{RESULTS}/03_membrane/{{pep}}/system.top",
        ndx=f"{RESULTS}/04_md/{{pep}}/index.ndx",
    output:
        dhdl=f"{RESULTS}/06_abfe/{{pep}}/lambda{{lam}}/dhdl.xvg",
    params:
        ns=config["abfe"]["ns_per_window"],
        nwin=config["abfe"]["n_lambda_windows"],
        stride=config["abfe"]["subsample_stride"],
    log:
        f"{LOGS}/06_abfe_{{pep}}_lambda{{lam}}.log",
    resources:
        runtime=720,
        slurm_partition="long",
        gres="gpu:1",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/06_abfe.py window "
        "--pep {wildcards.pep} --lambda-index {wildcards.lam} "
        "--n-windows {params.nwin} --ns {params.ns} "
        "--subsample {params.stride} "
        "--outdir $(dirname {output.dhdl}) > {log} 2>&1"


rule abfe_mbar:
    """MBAR over every window, with a convergence assessment.

    The instructions require a free energy reported WITH its convergence, not a
    bare number: peptides are flexible and often charged, convergence is the
    known weak point, and an honest uncertainty is worth more here than a
    confident value.
    """
    input:
        lambda w: expand(
            f"{RESULTS}/06_abfe/{w.pep}/lambda{{lam}}/dhdl.xvg", lam=LAMBDAS),
    output:
        dg=f"{RESULTS}/06_abfe/{{pep}}/dg.tsv",
        convergence=f"{RESULTS}/06_abfe/{{pep}}/convergence.txt",
    log:
        f"{LOGS}/06_mbar_{{pep}}.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/06_abfe.py mbar "
        "--dir $(dirname {output.dg}) --out {output.dg} "
        "--convergence {output.convergence} > {log} 2>&1"


rule abfe_collect:
    input:
        expand(f"{RESULTS}/06_abfe/{{pep}}/dg.tsv", pep=ABFE_PEPTIDES),
    output:
        table=f"{RESULTS}/06_abfe/dg.tsv",
    log:
        f"{LOGS}/06_collect.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/06b_collect.py "
        "--dir $(dirname {output.table}) --out {output.table} > {log} 2>&1"
