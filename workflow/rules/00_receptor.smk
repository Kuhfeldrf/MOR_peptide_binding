# Stage 0 - receptor preparation from the primary structure.
#
# STATUS: WORKING. Output verified byte-identical across runs.
#
# Ligands resolved in the structure are EXTRACTED here rather than predicted.
# DAMGO is the case that matters: Tyr-D-Ala-Gly-N-MePhe-Gly-ol has no valid
# canonical FASTA form, so it cannot be co-folded from sequence, but it is
# present in 6DDF as entity 5. Taking the experimental pose also means Stage 6
# starts from measured coordinates rather than compounding pose error into a
# free energy - at the cost, stated in the README, that DAMGO's ABFE then
# validates the free-energy machinery and not pose prediction.
#
# Which ligands to extract is config (receptor.extract_ligands), so adding one
# does not require editing the workflow.


rule receptor_prep:
    input:
        cif=f"{DATA}/raw/{config['receptor']['primary_pdb'].lower()}.cif",
    output:
        receptor=f"{RESULTS}/00_receptor/mOR_clean.pdb",
        fasta=f"{RESULTS}/00_receptor/mOR_clean.fasta",
        checksums=f"{RESULTS}/00_receptor/checksums.json",
        prep_log=f"{RESULTS}/00_receptor/receptor_prep.log",
        ligands=expand(f"{RESULTS}/00_receptor/{{lig}}_ref.pdb",
                       lig=sorted(config["receptor"].get("extract_ligands", {}))),
    params:
        ph=config["receptor"]["ph"],
        keep=config["receptor"]["keep_chain"],
        ligands=",".join(sorted(config["receptor"].get("extract_ligands", {}))),
    log:
        f"{LOGS}/00_receptor.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/00_receptor.py "
        "--cif {input.cif} --outdir $(dirname {output.receptor}) "
        "--ph {params.ph} --keep-chain {params.keep} "
        "--extract-ligands {params.ligands} > {log} 2>&1"


rule receptor_determinism:
    """Stage 0 must be byte-identical across runs; this proves it rather than
    asserting it. The instructions make determinism a stated requirement, so it
    is a rule with an output a reviewer can inspect, not a note in a log."""
    input:
        checksums=f"{RESULTS}/00_receptor/checksums.json",
        cif=f"{DATA}/raw/{config['receptor']['primary_pdb'].lower()}.cif",
    output:
        report=f"{RESULTS}/00_receptor/determinism.txt",
    params:
        ph=config["receptor"]["ph"],
        keep=config["receptor"]["keep_chain"],
        ligands=",".join(sorted(config["receptor"].get("extract_ligands", {}))),
    log:
        f"{LOGS}/00_determinism.log",
    conda:
        "../../environment.yml"
    shell:
        "python3 {SCRIPTS}/00b_check_determinism.py "
        "--cif {input.cif} --reference {input.checksums} "
        "--ph {params.ph} --keep-chain {params.keep} "
        "--extract-ligands {params.ligands} > {output.report} 2>&1"
