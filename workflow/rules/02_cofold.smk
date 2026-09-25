# Stage 2 - Chai-1 co-folding, single-sequence mode, PER PEPTIDE PER SEED.
#
# STATUS: WORKING.
#
# Peptides are keyed by the `name` column of the curated reference table
# throughout the workflow. The screening library uses its own `peptide_id`
# scheme with a `real__` prefix; joining the two on name rather than carrying
# both conventions is what keeps the wildcards consistent from here to Stage 7.
#
# INTERPRETATION NOTE for the README: these are learned-confidence scores, not
# binding affinity. Stage 2 found ipTM flat at 0.27-0.35 across a 32-fold
# affinity range, and ipSAE no better (D21). Cross-seed agreement was the only
# diagnostic that separated a correct pose from a failed one, which is the
# stated reason the pipeline continues into physics.


rule cofold_seed:
    """One Chai-1 inference: one peptide, one seed, five diffusion samples."""
    input:
        receptor=f"{RESULTS}/00_receptor/mOR_clean.fasta",
        library=f"{RESULTS}/01_library/peptides.tsv",
    output:
        scores=f"{RESULTS}/02_cofold/{{pep}}/seed{{seed}}/scores.tsv",
        meta=f"{RESULTS}/02_cofold/{{pep}}/seed{{seed}}/_cand_meta.npz",
    params:
        msa=config["cofold"]["use_msa"],
        recycles=config["cofold"]["num_trunk_recycles"],
        steps=config["cofold"]["num_diffn_timesteps"],
    log:
        f"{LOGS}/02_cofold_{{pep}}_seed{{seed}}.log",
    resources:
        runtime=90,
        mem_mb=64000,
        gres="gpu:1",
    conda:
        "../../envs/chai.yml"
    shell:
        "python3 {SCRIPTS}/02_cofold.py "
        "--receptor {input.receptor} --library {input.library} "
        "--peptide {wildcards.pep} --seed {wildcards.seed} "
        "--outdir $(dirname {output.scores}) "
        "--recycles {params.recycles} --timesteps {params.steps} > {log} 2>&1"


rule cofold_best_pose:
    """Top-ranked model for this peptide, as the Stage 3 starting structure."""
    input:
        seeds=lambda w: expand(
            f"{RESULTS}/02_cofold/{w.pep}/seed{{seed}}/scores.tsv", seed=SEEDS),
    output:
        pose=f"{RESULTS}/02_cofold/{{pep}}/best_pose.pdb",
    log:
        f"{LOGS}/02_bestpose_{{pep}}.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/02d_best_pose.py "
        "--dir $(dirname {output.pose}) --out {output.pose} > {log} 2>&1"


rule cofold_scores:
    """Aggregate confidence metrics across every peptide and seed."""
    input:
        expand(f"{RESULTS}/02_cofold/{{pep}}/seed{{seed}}/scores.tsv",
               pep=COFOLD, seed=SEEDS),
    output:
        table=f"{RESULTS}/02_cofold/scores.tsv",
    log:
        f"{LOGS}/02_scores.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/02e_aggregate.py "
        "--cofold $(dirname {output.table}) --out {output.table} > {log} 2>&1"


rule cofold_agreement:
    """Cross-seed agreement: do the seeds converge on the same pose?

    Superposes on RECEPTOR CA only, then measures the peptide in that frame -
    the question is whether the peptide lands in the same place on the same
    receptor, which a whole-complex superposition would not answer (D23).
    """
    input:
        expand(f"{RESULTS}/02_cofold/{{pep}}/seed{{seed}}/scores.tsv",
               pep=COFOLD, seed=SEEDS),
    output:
        table=f"{RESULTS}/02_cofold/cross_seed_agreement.tsv",
    log:
        f"{LOGS}/02_agreement.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/02_agreement.py "
        "--cofold $(dirname {output.table}) --outdir $(dirname {output.table}) "
        "> {log} 2>&1"


rule cofold_ipsae:
    """ipSAE via the vendored reference implementation (D21).

    Reported because the instructions require it, and because its failure at
    this peptide length is informative: d0 is clamped to 1.0 A for interfaces
    of 27 residues or fewer, so a 4-10 residue peptide cannot score well
    however correct the pose. It is not used to rank.
    """
    input:
        expand(f"{RESULTS}/02_cofold/{{pep}}/seed{{seed}}/_cand_meta.npz",
               pep=COFOLD, seed=SEEDS),
    output:
        table=f"{RESULTS}/02_cofold/ipsae.tsv",
    log:
        f"{LOGS}/02_ipsae.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/02c_ipsae.py "
        "--cofold $(dirname {output.table}) --out {output.table} > {log} 2>&1"
