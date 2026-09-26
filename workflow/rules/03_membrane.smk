# Stage 3 - membrane system build, PER PEPTIDE.
#
# STATUS: WORKING.
#
# This is where the workflow's fan-out lives: each peptide gets its own bilayer,
# parameterisation and topology, and Snakemake schedules them independently
# through the SLURM executor. Scaling from 3 peptides to 300 is a config change.
#
# The starting pose comes from one of two places. DAMGO and any other ligand
# resolved in the receptor structure are EXTRACTED from experiment; everything
# else uses its best-ranked Chai-1 pose. Which applies is decided by config, not
# by a branch buried in a script.
#
# Both rules below want the PEPTIDE ALONE, so both use ligand_only. There used
# to be a second function, complex_source, returning Chai-1's best_pose.pdb -
# the receptor AND peptide together - and membrane_complex passed it as
# --ligand. 03a_prep_complex.py relabels every atom it is handed as residue
# LIG, so each co-folded system got the receptor twice: once properly as chain
# R, and again as a 2294-atom "ligand" sitting at Chai-1's coordinates ~200 A
# away. The bounding box grew to 251 x 251 x 303 A and packmol dutifully filled
# it with 1.9M atoms instead of 84k. Nothing failed; the jobs just ran for
# hours. DAMGO was unaffected because for an experimental ligand both functions
# returned the same peptide-only file, so the bug could not show up in the one
# system that had been checked by hand.
#
# complex_source is deleted rather than fixed: two nearly identical functions
# whose names differ by what they include is the footgun that caused this.


def ligand_only(wildcards):
    """The PEPTIDE ALONE - never receptor+peptide. See the note above."""
    if wildcards.pep in EXPERIMENTAL:
        return f"{RESULTS}/00_receptor/{wildcards.pep}_ref.pdb"
    return f"{RESULTS}/02_cofold/{wildcards.pep}/best_ligand.pdb"


rule membrane_complex:
    """Receptor + peptide in the OPM membrane frame, AMBER-named."""
    input:
        receptor=f"{RESULTS}/00_receptor/mOR_clean.pdb",
        ligand=ligand_only,
        opm=f"{DATA}/raw/6ddf_opm.pdb",
    output:
        complex=f"{RESULTS}/03_membrane/{{pep}}/complex_opm.pdb",
        renames=f"{RESULTS}/03_membrane/{{pep}}/prep_renames.log",
    log:
        f"{LOGS}/03_complex_{{pep}}.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/03a_prep_complex.py "
        "--receptor {input.receptor} --ligand {input.ligand} "
        "--opm {input.opm} --outdir $(dirname {output.complex}) > {log} 2>&1"


rule ligand_params:
    """GAFF2 / AM1-BCC parameters for the peptide as a single molecule (D10)."""
    input:
        ligand=ligand_only,
    output:
        mol2=f"{RESULTS}/03_membrane/{{pep}}/params/ligand.mol2",
        frcmod=f"{RESULTS}/03_membrane/{{pep}}/params/ligand.frcmod",

    log:
        f"{LOGS}/03_params_{{pep}}.log",
    resources:
        runtime=60,
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/03b_ligand_params.py "
        "--ligand {input.ligand} --outdir $(dirname {output.mol2}) "
        "> {log} 2>&1"


rule membrane_pack:
    """Pack the bilayer with packmol-memgen.

    --pbc is mandatory: without it the patch overhangs the box and atoms
    overlap their own periodic images (D29). --apl_offset leaves room for the
    receptor footprint, which the default lipid count does not (D27).
    """
    input:
        complex=f"{RESULTS}/03_membrane/{{pep}}/complex_opm.pdb",
    output:
        packed=f"{RESULTS}/03_membrane/{{pep}}/membrane_system.pdb",
        packlog=f"{RESULTS}/03_membrane/{{pep}}/packmol.log",
    params:
        lipids=":".join(config["membrane"]["lipids"]),
        ratio=config["membrane"]["lipid_ratio"],
        apl=config["membrane"]["apl_offset"],
        dist=config["membrane"]["boundary_dist"],
        wat=config["membrane"]["water_layer"],
        salt=config["membrane"]["ionic_strength_mM"] / 1000.0,
        cation=config["membrane"]["cation"],
    log:
        f"{LOGS}/03_pack_{{pep}}.log",
    resources:
        runtime=720,
        mem_mb=64000,
    conda:
        "../../environment.yml"
    shell:
        "bash {SCRIPTS}/03_pack_membrane.sh "
        "{input.complex} {output.packed} "
        "{params.lipids} {params.ratio} {params.apl} "
        "{params.dist} {params.wat} {params.salt} {params.cation} "
        "> {log} 2>&1"


rule membrane_qc:
    """Reject a pack that did not converge or is not periodic (D28, D29)."""
    input:
        packed=f"{RESULTS}/03_membrane/{{pep}}/membrane_system.pdb",
    output:
        report=f"{RESULTS}/03_membrane/{{pep}}/qc.txt",
    params:
        chol=config["membrane"]["cholesterol_fraction"],
    log:
        f"{LOGS}/03_qc_{{pep}}.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/03c_membrane_qc.py "
        "--dir $(dirname {input.packed}) --expect-chol-frac {params.chol} "
        "> {output.report} 2>&1"


rule membrane_topology:
    """tleap: ff19SB + LIPID21 + OPC + GAFF2 ligand, then verified conversion."""
    input:
        packed=f"{RESULTS}/03_membrane/{{pep}}/membrane_system.pdb",
        qc=f"{RESULTS}/03_membrane/{{pep}}/qc.txt",
        mol2=f"{RESULTS}/03_membrane/{{pep}}/params/ligand.mol2",
        frcmod=f"{RESULTS}/03_membrane/{{pep}}/params/ligand.frcmod",
    output:
        parm=f"{RESULTS}/03_membrane/{{pep}}/system.parm7",
        rst=f"{RESULTS}/03_membrane/{{pep}}/system.rst7",
        top=f"{RESULTS}/03_membrane/{{pep}}/system.top",
        gro=f"{RESULTS}/03_membrane/{{pep}}/system.gro",
    log:
        f"{LOGS}/03_topology_{{pep}}.log",
    resources:
        runtime=120,
        mem_mb=32000,
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/03d_build_topology.py "
        "--dir $(dirname {output.parm}) > {log} 2>&1"
