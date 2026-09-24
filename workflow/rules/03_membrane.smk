# Stage 3 - POPC/cholesterol bilayer via packmol-memgen, 150 mM, OPM/PPM orientation.
#
# STATUS: STUBBED. The rule structure, inputs, and outputs are real; the
# script it calls exits non-zero with a STUBBED marker until implemented.

rule membrane_build:
    input:
        "results/02_cofold/scores.tsv"
    output:
        "results/03_membrane/system.parm7"
    log:
        "logs/membrane_build.log"
    conda:
        "../../environment.yml"
    script:
        "../../scripts/03_membrane.py"
