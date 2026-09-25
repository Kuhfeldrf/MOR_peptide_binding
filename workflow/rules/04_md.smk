# Stage 4 - membrane MD, PER PEPTIDE.
#
# STATUS: WORKING.
#
# minimise -> heat -> staged restraint release -> production, as the build
# instructions specify. Each peptide is an independent chain of jobs, so the
# whole set runs concurrently under the SLURM executor.


rule md_setup:
    """Index groups for thermostat coupling and the reduced trajectory group."""
    input:
        parm=f"{RESULTS}/03_membrane/{{pep}}/system.parm7",
        gro=f"{RESULTS}/03_membrane/{{pep}}/system.gro",
    output:
        ndx=f"{RESULTS}/04_md/{{pep}}/index.ndx",
        posre=f"{RESULTS}/04_md/{{pep}}/posre_system1.itp",
    params:
        proximal=config["md"]["proximal_cutoff"],
    log:
        f"{LOGS}/04_setup_{{pep}}.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/04a_setup_md.py "
        "--membrane $(dirname {input.parm}) --outdir $(dirname {output.ndx}) "
        "--proximal {params.proximal} > {log} 2>&1"


rule md_minimise:
    input:
        gro=f"{RESULTS}/03_membrane/{{pep}}/system.gro",
        top=f"{RESULTS}/03_membrane/{{pep}}/system.top",
        ndx=f"{RESULTS}/04_md/{{pep}}/index.ndx",
        mdp="config/min.mdp",
    output:
        gro=f"{RESULTS}/04_md/{{pep}}/min.gro",
    log:
        f"{LOGS}/04_min_{{pep}}.log",
    resources:
        runtime=90,
        slurm_extra="'--gres=gpu:1'",
    conda:
        "../../environment.yml"
    shell:
        "bash {SCRIPTS}/04_gmx_run.sh min {input.mdp} {input.gro} {input.gro} "
        "{input.top} {input.ndx} $(dirname {output.gro}) > {log} 2>&1"


rule md_heat:
    """NVT to 310 K, body temperature, with the solute restrained."""
    input:
        gro=f"{RESULTS}/04_md/{{pep}}/min.gro",
        top=f"{RESULTS}/03_membrane/{{pep}}/system.top",
        ndx=f"{RESULTS}/04_md/{{pep}}/index.ndx",
        mdp="config/nvt.mdp",
    output:
        gro=f"{RESULTS}/04_md/{{pep}}/nvt.gro",
        cpt=f"{RESULTS}/04_md/{{pep}}/nvt.cpt",
    log:
        f"{LOGS}/04_nvt_{{pep}}.log",
    resources:
        runtime=120,
        slurm_extra="'--gres=gpu:1'",
    conda:
        "../../environment.yml"
    shell:
        "bash {SCRIPTS}/04_gmx_run.sh nvt {input.mdp} {input.gro} {input.gro} "
        "{input.top} {input.ndx} $(dirname {output.gro}) > {log} 2>&1"


rule md_equilibrate:
    """Staged restraint release, 1000 -> 500 -> 100 -> 10 kJ/mol/nm^2.

    Releasing in one step lets the solute relax into a lipid and water
    arrangement that has not yet adapted to it.
    """
    input:
        gro=f"{RESULTS}/04_md/{{pep}}/nvt.gro",
        cpt=f"{RESULTS}/04_md/{{pep}}/nvt.cpt",
        top=f"{RESULTS}/03_membrane/{{pep}}/system.top",
        ndx=f"{RESULTS}/04_md/{{pep}}/index.ndx",
        mdp="config/npt.mdp",
    output:
        gro=f"{RESULTS}/04_md/{{pep}}/npt_final.gro",
        cpt=f"{RESULTS}/04_md/{{pep}}/npt_final.cpt",
        summary=f"{RESULTS}/04_md/{{pep}}/equilibration.txt",
    params:
        stages=" ".join(str(x) for x in config["md"]["restraint_stages"]),
    log:
        f"{LOGS}/04_npt_{{pep}}.log",
    resources:
        runtime=240,
        slurm_extra="'--gres=gpu:1'",
    conda:
        "../../environment.yml"
    shell:
        "bash {SCRIPTS}/04_gmx_equilibrate.sh $(dirname {output.gro}) "
        "{input.top} {input.ndx} {input.mdp} {params.stages} > {log} 2>&1"


rule md_production:
    input:
        gro=f"{RESULTS}/04_md/{{pep}}/npt_final.gro",
        cpt=f"{RESULTS}/04_md/{{pep}}/npt_final.cpt",
        top=f"{RESULTS}/03_membrane/{{pep}}/system.top",
        ndx=f"{RESULTS}/04_md/{{pep}}/index.ndx",
        mdp="config/prod.mdp",
    output:
        xtc=f"{RESULTS}/04_md/{{pep}}/prod.xtc",
        tpr=f"{RESULTS}/04_md/{{pep}}/prod.tpr",
    log:
        f"{LOGS}/04_prod_{{pep}}.log",
    resources:
        runtime=lambda wc, attempt: 60 * config["md"]["production_ns"],
        slurm_partition="long",
        slurm_extra="'--gres=gpu:1'",
    conda:
        "../../environment.yml"
    shell:
        "bash {SCRIPTS}/04_gmx_run.sh prod {input.mdp} {input.gro} {input.gro} "
        "{input.top} {input.ndx} $(dirname {output.xtc}) {input.cpt} > {log} 2>&1"


rule md_analysis:
    """Peptide RMSD, contact persistence, and whether it stayed in the pocket."""
    input:
        xtc=f"{RESULTS}/04_md/{{pep}}/prod.xtc",
        tpr=f"{RESULTS}/04_md/{{pep}}/prod.tpr",
        ndx=f"{RESULTS}/04_md/{{pep}}/index.ndx",
    output:
        stats=f"{RESULTS}/04_md/{{pep}}/md_stats.json",
        stability=f"{RESULTS}/04_md/{{pep}}/md_stability.tsv",
        contacts=f"{RESULTS}/04_md/{{pep}}/contact_persistence.tsv",
    params:
        offset=config["receptor"]["resid_offset"],
    log:
        f"{LOGS}/04_analysis_{{pep}}.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/04c_md_analysis.py "
        "--dir $(dirname {input.xtc}) --out $(dirname {output.stats}) "
        "--resid-offset {params.offset} > {log} 2>&1"
