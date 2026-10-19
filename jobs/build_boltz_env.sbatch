#!/usr/bin/env bash
#SBATCH -J boltz_env
#SBATCH -p short
#SBATCH -t 180
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/boltz_env_%j.log
#
# Build `mor-boltz` for Boltz-2 affinity prediction.
#
# WHY BOLTZ-2. This project showed that co-folding CONFIDENCE metrics (ipTM,
# iLIS, PAE) carry no affinity signal (n=36, 4.4 logs of cited Ki, length
# decorrelated, 88% power; iLIS rho=+0.196 ns). That is a result about
# confidence metrics, NOT about ML in general. Boltz-2 (Passaro et al. 2025,
# bioRxiv 2025.06.14.659707) has a head TRAINED on affinity and reports Pearson
# 0.62 on the held-out FEP+ benchmark - comparable to OpenFE FEP at ~1000x lower
# cost - and won the CASP16 affinity challenge. It is the contemporary method
# and it is the untested one.
#
# WHY IT MATTERS FOR THE HUANG COMPLEX: the machine delivers ~10x for AI
# inference vs ~2x for traditional simulation. An affinity predictor is an AI
# workload; brute-force FEP is not. A pipeline led by Boltz-2 uses the machine's
# actual strength.
#
# SEPARATE ENV REQUIRED: boltz 2.2.1 pins numpy<2.0; mor-protenix has numpy
# 2.4.1. Same incompatibility as decisions.md D5.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"
source /etc/profile.d/z00_lmod.sh
module load cuda/12.9.0-none-none || { echo "FATAL: cuda module"; exit 1; }
source "$ROOT/miniforge3/etc/profile.d/conda.sh"

export PIP_CACHE_DIR=$ROOT/.cache/pip
export TORCH_HOME=$ROOT/.cache/torch
mkdir -p "$PIP_CACHE_DIR" "$TORCH_HOME" "$ROOT/logs"

if [ ! -d "$ROOT/miniforge3/envs/mor-boltz" ]; then
  mamba create -y -n mor-boltz python=3.12 pip || exit 1
fi
conda activate mor-boltz
python -V

python -m pip install --no-input "boltz[cuda]==2.2.1"
rc=$?
echo "pip exit: $rc"
[ $rc -ne 0 ] && exit $rc

python - <<'PY'
import importlib
for m in ("boltz","torch","numpy","rdkit","pytorch_lightning"):
    try:
        mod=importlib.import_module(m)
        print(f"  {m:20s} {getattr(mod,'__version__','?')}")
    except Exception as e:
        print(f"  {m:20s} FAILED: {type(e).__name__}: {e}")
PY
echo "=== boltz CLI ==="
which boltz && boltz --help 2>&1 | head -40
echo "=== predict subcommand options (affinity flags) ==="
boltz predict --help 2>&1 | grep -iE "affinity|diffusion_samples|devices|accelerator|cache|method|sampling" | head -20
echo "=== $(date -Is) done ==="
