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
# away. The bounding box grew from 86 x 86 x 119 A to 251 x 251 x 303 A and
# packmol dutifully filled it with 1.9M atoms instead of 84k. Nothing failed;
# the jobs just ran for hours. DAMGO was unaffected because for an experimental
# ligand both functions returned the same peptide-only file, so the bug could
# not show up in the one system that had been checked by hand.
#
# complex_source is deleted rather than fixed: two nearly identical functions
# whose names differ by what they include is the footgun that caused this.


def ligand_only(wildcards):
    """The PEPTIDE ALONE - never receptor+peptide. See the note above."""
    if wildcards.pep in EXPERIMENTAL:
        return f"{RESULTS}/00_receptor/{wildcards.pep}_ref.pdb"
    return f"{RESULTS}/02_cofold/{wildcards.pep}/best_ligand.pdb"


def pose_reference(wildcards):
    """The co-folded receptor+peptide, used ONLY to relate two coordinate
    frames - never as the ligand itself.

    A co-folding model predicts in its own frame near the origin; the
    experimental receptor sits at its crystal coordinates. The prediction's own
    receptor is the only thing that relates the two, so it is superposed onto
    the experimental receptor and the resulting RIGID transform is applied to
    the peptide. The peptide therefore keeps exactly the pose that was
    predicted, in the correct frame.

    An experimental ligand was extracted from the receptor structure and is
    already in frame, so it needs no reference.
    """
    if wildcards.pep in EXPERIMENTAL:
        return []
    return f"{RESULTS}/02_cofold/{wildcards.pep}/best_pose.pdb"


def pose_flags(wildcards, input):
    if wildcards.pep in EXPERIMENTAL:
        return ""
    # The prediction numbers the construct from 1; the crystal starts at 65.
    # Superposing without that offset silently matches the wrong residues -
    # it gave a 24.6 A fit where the correct mapping gives 2.7 A.
    return (f"--pose {input.pose} "
            f"--pose-chain {config['cofold']['receptor_chain']} "
            f"--pose-offset {config['receptor']['resid_offset']}")


rule membrane_complex:
    """Receptor + peptide in the OPM membrane frame, AMBER-named."""
    input:
        receptor=f"{RESULTS}/00_receptor/mOR_clean.pdb",
        ligand=ligand_only,
        pose=pose_reference,
        opm=f"{DATA}/raw/6ddf_opm.pdb",
    params:
        pose=pose_flags,
    output:
        complex=f"{RESULTS}/03_membrane/{{pep}}/complex_opm.pdb",
        renames=f"{RESULTS}/03_membrane/{{pep}}/prep_renames.log",
        geometry=f"{RESULTS}/03_membrane/{{pep}}/pose_geometry.tsv",
    log:
        f"{LOGS}/03_complex_{{pep}}.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/03a_prep_complex.py "
        "--receptor {input.receptor} --ligand {input.ligand} {params.pose} "
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
        nloop_all=config["membrane"]["nloop_all"],
        nloop=config["membrane"]["nloop"],
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
        "{params.nloop_all} {params.nloop} "
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
    # The report is written to the LOG first, then copied to the output on
    # success. It used to be written straight to {output.report}, which meant a
    # FAILED pack destroyed its own explanation: Snakemake deletes the outputs
    # of a failed job, {log} stayed empty because nothing was redirected there,
    # and the only way to learn why QC rejected a system was to re-run the
    # script by hand. A rule whose whole purpose is to explain a rejection must
    # keep its reasoning on the failure path.
    shell:
        "python3 {SCRIPTS}/03c_membrane_qc.py "
        "--dir $(dirname {input.packed}) --expect-chol-frac {params.chol} "
        "> {log} 2>&1 && cp {log} {output.report}"


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
