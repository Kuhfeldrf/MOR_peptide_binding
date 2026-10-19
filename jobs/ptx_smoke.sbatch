#!/usr/bin/env bash
#SBATCH -J ptx_smoke
#SBATCH -p short
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 90
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/ptx_smoke_%j.log
#
# One (peptide, seed) through Protenix before the 105-job array is written.
#
# Purpose is NOT a result. It is to (a) confirm the model runs MSA-free on this
# receptor at all, (b) record the real output tree and confidence-JSON schema so
# scripts/02_cofold_protenix.py parses what Protenix actually writes rather than
# what its docs imply, and (c) measure runtime per job to size the array.
#
# TWO LESSONS ALREADY PAID FOR, ENCODED HERE:
#
# 1. CUDA_HOME. protenix JIT-compiles a fused LayerNorm CUDA extension at
#    IMPORT time - protenix/model/layer_norm/layer_norm.py does
#    `try: import fast_layer_norm_cuda_v2 / except ImportError: compile(...)`
#    with no environment switch to skip it. Without a CUDA toolkit even
#    `protenix --help` dies with "CUDA_HOME environment variable is not set".
#    The cuda/12.9.0 module sets CUDA_HOME itself. torch here is cu126; nvcc
#    12.9 compiles its extensions fine across 12.x.
#
# 2. EXIT CODES MUST NOT LIE. The first version of this script recorded
#    protenix's rc and carried on, so SLURM reported COMPLETED / 0:0 for a run
#    in which protenix exited 1 and produced no output at all. That is the same
#    false-success shape as the PDBFixer report fixed in b48ce60. Every step
#    below fails the job loudly.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"

source /etc/profile.d/z00_lmod.sh          # `module` is undefined otherwise
module load cuda/12.9.0-none-none || { echo "FATAL: cuda module"; exit 1; }
source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-protenix || { echo "FATAL: conda activate"; exit 1; }

# Weights and the CCD cache are several GB; keep them off the 92 G home quota.
export PROTENIX_ROOT_DIR=$ROOT/.cache/protenix
export TORCH_HOME=$ROOT/.cache/torch
# Pinned rather than autodetected, so the compiled kernel is reproducible and
# does not silently depend on which GPU the build landed on. L40S is Ada, sm_89.
export TORCH_CUDA_ARCH_LIST=8.9
mkdir -p "$PROTENIX_ROOT_DIR" "$TORCH_HOME"

echo "CUDA_HOME=$CUDA_HOME"
nvcc --version | tail -1
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# --- Pre-compile the fused kernel ONCE, here, before any array is launched ---
# build_directory resolves to site-packages when it is writable (it is: the env
# lives in scratch and we own it). That is a SHARED path, so 105 concurrent
# array tasks would each race ninja into the same directory. Compiling it here
# makes every later task hit the plain `import` branch instead.
echo "=== pre-compiling fast_layer_norm_cuda_v2 (first import only) ==="
t0=$SECONDS
python -c "
from protenix.model.layer_norm.layer_norm import FusedLayerNorm
import fast_layer_norm_cuda_v2 as m
print('compiled OK:', m.__file__)
" || { echo "FATAL: fused LayerNorm would not build"; exit 1; }
echo "compile took $((SECONDS-t0))s"

echo; echo "=== protenix pred --help (authoritative flag list) ==="
protenix pred --help 2>&1 | head -100

OUT=$ROOT/results/02_cofold_ptx_smoke
rm -rf "$OUT"; mkdir -p "$OUT"

PEP=pos__beta_endorphin_1_10
RECEPTOR=$(sed -n '2p' results/00_receptor/mOR_clean.fasta)
PEPSEQ=$(awk -F'\t' -v p="$PEP" 'NR>1 && $2==p {print $3}' results/01_library/length_controls.tsv)
echo "receptor ${#RECEPTOR} aa ; peptide $PEP = $PEPSEQ (${#PEPSEQ} aa)"
[ -z "$PEPSEQ" ] && { echo "FATAL: peptide sequence not found"; exit 1; }

# No pairedMsaPath / unpairedMsaPath / templatesPath: MSA-free single-sequence,
# matching the Chai-1 arm (use_msa_server=False). docs/infer_json_format.md
# states all three are optional and inference proceeds without them.
cat > "$OUT/input.json" <<JSON
[
  {
    "name": "${PEP}_seed0",
    "sequences": [
      {"proteinChain": {"sequence": "${RECEPTOR}", "count": 1}},
      {"proteinChain": {"sequence": "${PEPSEQ}", "count": 1}}
    ]
  }
]
JSON

echo; echo "=== run ==="
t0=$SECONDS
protenix pred \
  -i "$OUT/input.json" \
  -o "$OUT/pred" \
  -s 0 \
  -n protenix_base_default_v1.0.0 \
  --use_msa false \
  --use_default_params false \
  --cycle 3 \
  --step 200
rc=$?
echo "protenix exit=$rc  elapsed=$((SECONDS-t0))s"
[ $rc -ne 0 ] && { echo "FATAL: protenix pred failed"; exit $rc; }

echo; echo "=== output tree ==="
find "$OUT/pred" -maxdepth 4 | head -60
echo; echo "=== json files found ==="
find "$OUT/pred" -name '*.json' | head -20
echo; echo "=== cif/pdb files found ==="
find "$OUT/pred" \( -name '*.cif' -o -name '*.pdb' \) | head -20

echo; echo "=== every confidence json, keys and scalar values ==="
python - "$OUT/pred" <<'PY'
import json, pathlib, sys, numpy as np
root = pathlib.Path(sys.argv[1])
for p in sorted(root.rglob("*.json")):
    print(f"\n--- {p.relative_to(root)}")
    try:
        d = json.load(open(p))
    except Exception as e:
        print("   unreadable:", e); continue
    if not isinstance(d, dict):
        print("   type:", type(d).__name__); continue
    for k, v in d.items():
        a = np.asarray(v)
        if a.ndim == 0:
            print(f"   {k:30s} = {v}")
        else:
            print(f"   {k:30s} shape={a.shape} dtype={a.dtype}")
PY
echo; echo "=== $(date -Is) done ==="
