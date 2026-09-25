#!/usr/bin/env bash
# Launch the pipeline under the SLURM executor.
#
# Snakemake stays on the login node and only orchestrates: it submits a SLURM
# job per rule instance and waits. That is why it runs detached rather than as
# a job of its own.
#
# PREREQUISITE: build the conda environments first, in a job with enough memory
#
#     sbatch build_envs.sbatch
#
# The libmamba solve for environment.yml is OOM-killed on the login node -
# AmberTools makes it a large solve - and Snakemake reports that as
# CreateCondaEnvironmentException with output truncated at "Solving
# environment", mentioning nothing about memory.
set -uo pipefail

ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"

source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-pilot

# The GROMACS wrappers read this rather than hardcoding a spack path whose
# directory name carries a build hash.
export GMX_PREFIX=$ROOT/spack-install/linux-zen4/gromacs-2025.3-45hxd2evh2mxcjqopqa56cr3atcbwlrv
export PATH="$GMX_PREFIX/bin:$PATH"

mkdir -p logs

echo "=== $(date -Is) launching pipeline ==="
echo "GMX_PREFIX=$GMX_PREFIX"
gmx --version 2>/dev/null | grep -m1 "GROMACS version" || echo "WARNING: gmx not on PATH"

# A killed Snakemake leaves the working directory locked, and the next launch
# then fails with a LockException that reads like a workflow error. Clearing it
# is safe because this script is the only thing that starts the workflow.
snakemake -s workflow/Snakefile --cores 1 --unlock >/dev/null 2>&1 || true

snakemake -s workflow/Snakefile --profile config/slurm --rerun-incomplete --keep-going
rc=$?

echo "=== $(date -Is) snakemake exited $rc ==="
exit $rc
