#!/bin/bash
#SBATCH --job-name=mor_envs
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/envs_%j.log
#SBATCH --error=logs/envs_%j.log
#
# Pre-build the workflow's conda environments.
#
# WHY THIS IS A JOB AND NOT A LOGIN-NODE STEP: the libmamba solve for the main
# environment is OOM-killed on the login node. AmberTools makes it a large
# solve and the solver needs several GB. Snakemake would otherwise attempt it
# in-process on the login node and fail with CreateCondaEnvironmentException
# whose output stops at "Solving environment", giving no indication that memory
# was the cause.
#
# --conda-create-envs-only builds the environments and submits no work, so the
# pipeline launch afterwards finds them ready.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"
source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-pilot

echo "=== $(date -Is) building workflow conda environments ==="
free -g | head -2

snakemake -s workflow/Snakefile \
          --profile config/slurm \
          --conda-create-envs-only \
          --cores 8
rc=$?
echo "exit: $rc"

echo "=== environments present ==="
ls -d .snakemake-conda/*/ 2>/dev/null | while read d; do
  printf "  %-60s %s\n" "$d" "$(du -sh "$d" 2>/dev/null | cut -f1)"
done
