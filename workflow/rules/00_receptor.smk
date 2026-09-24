# Stage 0 - receptor preparation from 6DDF.
#
# STATUS: WORKING. Output verified byte-identical across runs.
#
# Also emits the experimentally resolved DAMGO ligand (chain D, entity 5),
# which Stage 6 uses as its ABFE starting structure rather than a predicted
# pose. See docs/damgo_notes.md.

rule receptor_prep:
    input:
        cif="data/raw/6ddf.cif",
    output:
        receptor="results/00_receptor/mOR_clean.pdb",
        damgo="results/00_receptor/damgo_ref.pdb",
        checksums="results/00_receptor/checksums.json",
        log="results/00_receptor/receptor_prep.log",
    params:
        ph=config["receptor"]["ph"],
    conda:
        "../../environment.yml"
    shell:
        "python3 scripts/00_receptor.py "
        "--cif {input.cif} --outdir results/00_receptor --ph {params.ph}"
