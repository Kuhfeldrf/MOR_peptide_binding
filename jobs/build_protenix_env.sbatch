#!/usr/bin/env bash
#SBATCH -J ptx_env
#SBATCH -p short
#SBATCH -t 120
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/ptx_env_%j.log
#
# Build the `mor-protenix` environment for the Stage 2 Protenix arm.
#
# WHY A SEPARATE ENVIRONMENT: protenix 2.0.0 pins numpy==2.4.1 and
# torch==2.7.1. `mor-chai` is pinned to numpy<2 (chai-lab 0.6.1, torch 2.6.0).
# The two cannot share an environment, for the same reason mor-chai and
# mor-pilot cannot - see docs/decisions.md D5.
#
# WHY A JOB AND NOT A LOGIN-NODE STEP: see build_envs.sbatch. Large solves are
# OOM-killed on the login node.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"
source "$ROOT/miniforge3/etc/profile.d/conda.sh"

# Keep every cache off the 92 G home quota.
export PIP_CACHE_DIR=$ROOT/.cache/pip
export TORCH_HOME=$ROOT/.cache/torch
mkdir -p "$PIP_CACHE_DIR" "$TORCH_HOME" "$ROOT/logs"

echo "=== $(date -Is) building mor-protenix ==="
free -g | head -2

if [ ! -d "$ROOT/miniforge3/envs/mor-protenix" ]; then
  # protenix requires_python >=3.11. 3.11 rather than 3.12 because several
  # pinned deps (pdbeccdutils, scikit-learn-extra) lag on wheels.
  mamba create -y -n mor-protenix python=3.11 pip || exit 1
fi
conda activate mor-protenix
python -V

# Installed from PyPI with its own pins honoured; no --no-deps, because the
# point of this arm is to run the model as its authors ship it.
python -m pip install --no-input protenix==2.0.0
rc=$?
echo "pip exit: $rc"
[ $rc -ne 0 ] && exit $rc

echo "=== versions ==="
python - <<'PY'
import importlib
for m in ("protenix", "torch", "numpy", "deepspeed", "rdkit", "gemmi", "biotite"):
    try:
        mod = importlib.import_module(m)
        print(f"  {m:12s} {getattr(mod, '__version__', '?')}")
    except Exception as e:
        print(f"  {m:12s} IMPORT FAILED: {type(e).__name__}: {e}")
PY

echo "=== protenix CLI surface (needed to settle the MSA question) ==="
which protenix
protenix --help 2>&1 | head -40
echo "--- pred --help ---"
protenix pred --help 2>&1 | head -80

echo "=== env size ==="
du -sh "$ROOT/miniforge3/envs/mor-protenix" 2>/dev/null
echo "=== $(date -Is) done ==="
