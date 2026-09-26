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

# Only one launch may be live at a time.
#
# This script used to run `--unlock` unconditionally, on the reasoning that
# clearing the lock "is safe because this script is the only thing that starts
# the workflow." That reasoning was wrong, and it conflated two different
# claims: that only this script launches the workflow, and that only one launch
# is ever active. Snakemake's directory lock exists precisely to enforce the
# second, so clearing it on every launch meant a new instance silently STOLE
# the lock from a running one instead of being refused by it.
#
# Three instances ended up live against this directory at once. Each of the
# seven benchmark peptides got two concurrent membrane_pack jobs writing the
# same output directory - packmol-memgen cd's into that directory and writes
# fixed filenames - so both copies of every output were interleaved garbage,
# and nothing anywhere reported a problem.
#
# flock is the guard the comment above wrongly assumed. A second launch now
# exits immediately and says which PID holds the lock.
exec 9>"$ROOT/.pipeline.lock"
if ! flock -n 9; then
  echo "ERROR: a pipeline launch is already running (pid $(cat "$ROOT/.pipeline.pid" 2>/dev/null || echo '?'))." >&2
  echo "       Wait for it, or stop it before launching again." >&2
  exit 1
fi
echo $$ > "$ROOT/.pipeline.pid"

# Safe only now: holding the flock proves no other launch is live, so any
# Snakemake lock still present is genuinely stale from a killed run.
snakemake -s workflow/Snakefile --cores 1 --unlock >/dev/null 2>&1 || true

snakemake -s workflow/Snakefile --profile config/slurm --rerun-incomplete --keep-going
rc=$?

echo "=== $(date -Is) snakemake exited $rc ==="
exit $rc
