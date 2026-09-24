#!/usr/bin/env bash
# Step 3 rebuild: two environments, with verification that fails loudly.
#
# The previous build reported COMPLETED while being unusable. This script
# exits non-zero if any required import fails, so job status is meaningful.
set -eo pipefail
PILOT_ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$PILOT_ROOT"
export PATH="$PILOT_ROOT/miniforge3/bin:$PATH"
source "$PILOT_ROOT/miniforge3/etc/profile.d/conda.sh"

echo "=== $(date -Is) rebuilding split environments ==="

# ---------------------------------------------------------------- main env
echo "--- removing contaminated mor-pilot ---"
conda env remove -n mor-pilot -y 2>/dev/null || true

echo "--- creating mor-pilot (main) ---"
mamba env create -n mor-pilot -f environment.yml

echo "--- creating mor-chai (Stage 2) ---"
conda env remove -n mor-chai -y 2>/dev/null || true
mamba env create -n mor-chai -f envs/chai.yml

# ---------------------------------------------------------------- verify
FAIL=0

echo
echo "=== VERIFY mor-pilot ==="
conda activate mor-pilot
python - <<'PY' || FAIL=1
import importlib, sys
req = ["numpy","scipy","pandas","MDAnalysis","mdtraj","pymbar","alchemlyb",
       "parmed","Bio","gemmi","matplotlib","seaborn","pyarrow"]
bad = []
for m in req:
    try:
        mod = importlib.import_module(m)
        print(f"  OK   {m}: {getattr(mod,'__version__','?')}")
    except Exception as e:
        print(f"  FAIL {m}: {type(e).__name__}: {e}")
        bad.append(m)
import numpy
print(f"  numpy = {numpy.__version__}")
if int(numpy.__version__.split('.')[0]) < 2:
    print("  FAIL numpy must be >=2 in the main env"); bad.append("numpy<2")
sys.exit(1 if bad else 0)
PY

for exe in packmol-memgen pdb2pqr30 snakemake; do
  if command -v "$exe" >/dev/null 2>&1; then echo "  OK   $exe -> $(command -v $exe)"
  else echo "  FAIL $exe not on PATH"; FAIL=1; fi
done
packmol-memgen --help >/dev/null 2>&1 && echo "  OK   packmol-memgen runs" || { echo "  FAIL packmol-memgen does not run"; FAIL=1; }
conda list --explicit > docs/env_solved_main.txt 2>/dev/null || conda list > docs/env_solved_main.txt
conda deactivate

echo
echo "=== VERIFY mor-chai ==="
conda activate mor-chai
python - <<'PY' || FAIL=1
import sys
try:
    import numpy, torch, chai_lab
    print(f"  OK   numpy {numpy.__version__}")
    print(f"  OK   torch {torch.__version__}")
    print(f"  OK   chai_lab {getattr(chai_lab,'__version__','?')}")
    print(f"  cuda available (no GPU on build node is expected): {torch.cuda.is_available()}")
    print(f"  bf16 supported flag: {torch.cuda.is_bf16_supported() if torch.cuda.is_available() else 'n/a (no GPU here)'}")
except Exception as e:
    print(f"  FAIL {type(e).__name__}: {e}"); sys.exit(1)
PY
conda list --explicit > docs/env_solved_chai.txt 2>/dev/null || conda list > docs/env_solved_chai.txt
conda deactivate

echo
if [ "$FAIL" -ne 0 ]; then
  echo "=== $(date -Is) BUILD FAILED VERIFICATION ==="
  exit 1
fi
echo "=== $(date -Is) BOTH ENVIRONMENTS VERIFIED ==="
