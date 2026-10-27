#!/usr/bin/env bash
#SBATCH -J anc_env
#SBATCH -p short
#SBATCH -t 180
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/MOR_demo_analysis/logs/anchor_env_%j.log
#
# Anchor-set arm B, step 0: rebuild `mor-protenix`.
#
# WHY THIS EXISTS AT ALL: /scratch/kuhfeldr-Kuhfeld_temp was cleared on
# 2026-10-26. The pilot tree, both conda envs and the Protenix/Chai weight
# caches are gone; only a fresh miniforge3 and the separate MOR_pharmacophore
# project remain. Nothing was lost - the repo is intact at
# github.com/Kuhfeldrf/MOR_peptide_binding - but every arm B prerequisite has
# to be rebuilt before an array can run.
#
# LAYOUT CHANGE, deliberate: the old jobs used ROOT=/scratch/kuhfeldr-Kuhfeld_temp
# with the repo AT the scratch root. That root is now shared with
# MOR_pharmacophore, so this repo lives in a subdirectory and REPO != CONDA_ROOT.
# The shared miniforge3 is REUSED rather than reinstalled: adding an env beside
# mor-pharm does not touch mor-pharm.
#
# Pins are not restated here - `protenix==2.0.0` carries its own, and the point
# of this arm is to run the model as its authors ship it. See
# jobs/build_protenix_env.sbatch for why it cannot share mor-chai's environment.
set -uo pipefail
REPO=/scratch/kuhfeldr-Kuhfeld_temp/MOR_demo_analysis
CONDA_ROOT=/scratch/kuhfeldr-Kuhfeld_temp/miniforge3
cd "$REPO"
source "$CONDA_ROOT/etc/profile.d/conda.sh"

# Every cache on scratch, off the 90 G home quota.
export PIP_CACHE_DIR=$REPO/.cache/pip
export TORCH_HOME=$REPO/.cache/torch
mkdir -p "$PIP_CACHE_DIR" "$TORCH_HOME" "$REPO/logs"

echo "=== $(date -Is) building mor-protenix ==="
free -g | head -2
df -h /scratch | tail -1

if [ ! -d "$CONDA_ROOT/envs/mor-protenix" ]; then
  mamba create -y -n mor-protenix python=3.11 pip || exit 1
fi
conda activate mor-protenix || { echo "FATAL: conda activate"; exit 1; }
python -V

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

# gemmi is all Stage 0 needs, so this env also covers receptor prep and the
# gate's RMSD check - no second environment build on the co-folding path.
python -c "import gemmi; print('gemmi OK - Stage 0 can run in this env')" || exit 1

echo "=== env size ==="
du -sh "$CONDA_ROOT/envs/mor-protenix" 2>/dev/null
echo "=== $(date -Is) done ==="
