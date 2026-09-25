#!/usr/bin/env bash
# Launch the full pipeline under the SLURM executor.
#
# Snakemake itself stays on the login node and only orchestrates: it submits a
# SLURM job per rule instance and waits. That is why it runs detached rather
# than as a job of its own.
#
# The first action will be building the two conda environments into
# .snakemake-conda. That is slow - the Chai-1 environment pulls torch - but it
# is the reproducible path: it proves environment.yml and envs/chai.yml
# actually describe a working environment, rather than relying on the two I
# built by hand.
set -uo pipefail

ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"

source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-pilot

# The gmx wrappers read this rather than hardcoding a spack path with a hash.
export GMX_PREFIX=$ROOT/spack-install/linux-zen4/gromacs-2025.3-45hxd2evh2mxcjqopqa56cr3atcbwlrv
export PATH="$GMX_PREFIX/bin:$PATH"

mkdir -p logs

echo "=== $(date -Is) launching pipeline ==="
echo "GMX_PREFIX=$GMX_PREFIX"
gmx --version 2>/dev/null | head -1 || echo "WARNING: gmx not on PATH"

# A killed Snakemake leaves a lock on the working directory. Clearing it here
# is safe because this script is the only thing that launches the workflow, and
# the alternative is a LockException that looks like a workflow error.
snakemake -s workflow/Snakefile --profile config/slurm --unlock >/dev/null 2>&1 || true

snakemake -s workflow/Snakefile \\
          --profile config/slurm \\
          --rerun-incomplete \
          --keep-going \
          2>&1

echo "=== $(date -Is) snakemake exited $? ==="
