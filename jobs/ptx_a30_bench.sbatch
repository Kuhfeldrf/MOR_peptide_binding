#!/usr/bin/env bash
#SBATCH -J ptx_a30
#SBATCH -p short
#SBATCH --gres=gpu:a30:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 120
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/ptx_a30_%j.log
#
# Two jobs in one, both prerequisites for sizing the Protenix array.
#
# 1. FIX A PORTABILITY BUG I INTRODUCED. jobs/ptx_smoke.sbatch pinned
#    TORCH_CUDA_ARCH_LIST=8.9 for reproducibility. 8.9 is Ada (L40S) ONLY;
#    A30 is Ampere, sm_80. The compiled fast_layer_norm_cuda_v2.so in
#    site-packages therefore has no A30-compatible code, and every A30 task
#    would either fail to load it or silently re-JIT - and re-JIT from 45
#    concurrent tasks races ninja into one shared directory. Recompile ONCE
#    here as a fat binary covering both architectures.
#
# 2. MEASURE THE A30:L40S RATIO. ORCA has 24 L40S but 76 A30, and right now far
#    more A30 are free. Whether to use them depends on the speed ratio, which is
#    NOT obvious from spec sheets: A30 has lower bf16 throughput than L40S but
#    HIGHER memory bandwidth (HBM2 ~933 GB/s vs GDDR6 ~864 GB/s), and at
#    N_token=291 this workload may well be bandwidth-bound. So measure it rather
#    than assume. Same peptide and seed as the L40S smoke run (191366, 1201 s)
#    to make it like-for-like.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"

source /etc/profile.d/z00_lmod.sh
module load cuda/12.9.0-none-none || { echo "FATAL: cuda module"; exit 1; }
source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-protenix || { echo "FATAL: conda activate"; exit 1; }

export PROTENIX_ROOT_DIR=$ROOT/.cache/protenix
export TORCH_HOME=$ROOT/.cache/torch
# Fat binary: 8.0 = A30 (Ampere), 8.9 = L40S (Ada). Covers every GPU on ORCA,
# so the array can schedule onto either without re-JIT.
export TORCH_CUDA_ARCH_LIST="8.0;8.9"

SO=$(python -c "import protenix.model.layer_norm as m, os; print(os.path.dirname(m.__file__))")/fast_layer_norm_cuda_v2.so
echo "existing kernel: $SO"
if [ -f "$SO" ]; then
  echo "architectures currently in it:"
  cuobjdump --list-elf "$SO" 2>/dev/null | head -5 || echo "  (cuobjdump unavailable)"
  mv "$SO" "${SO}.sm89only.bak"
  echo "moved aside, forcing fat recompile"
fi

echo "=== recompiling for 8.0;8.9 ==="
t0=$SECONDS
python -c "
from protenix.model.layer_norm.layer_norm import FusedLayerNorm
import fast_layer_norm_cuda_v2 as m
print('compiled OK:', m.__file__)
" || { echo "FATAL: fat recompile failed"; exit 1; }
echo "compile took $((SECONDS-t0))s"
echo "architectures now present:"
cuobjdump --list-elf "$SO" 2>/dev/null | head -5 || true

nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# Same input as the L40S smoke run, for a like-for-like timing comparison.
OUT=$ROOT/results/02_cofold_ptx_a30bench
rm -rf "$OUT"; mkdir -p "$OUT"
cp "$ROOT/results/02_cofold_ptx_smoke/input.json" "$OUT/input.json" \
  || { echo "FATAL: smoke input.json missing"; exit 1; }

echo; echo "=== run: same peptide+seed as L40S smoke (1201 s baseline) ==="
t0=$SECONDS
protenix pred -i "$OUT/input.json" -o "$OUT/pred" -s 0 \
  -n protenix_base_default_v1.0.0 \
  --use_msa false --use_default_params false --cycle 3 --step 200
rc=$?
W=$((SECONDS-t0))
echo "protenix exit=$rc  A30 wall=${W}s"
[ $rc -ne 0 ] && { echo "FATAL: protenix pred failed on A30"; exit $rc; }

python3 -c "
w=$W; l=1201
print(f'A30/L40S ratio = {w/l:.2f}x  (A30 {w}s vs L40S {l}s)')
print(f'105-task array, all-A30 at 36 concurrent: {3*w/3600:.2f} h wall')
print(f'105-task array, all-L40S at 20 concurrent: {6*l/3600:.2f} h wall')
"
echo "=== $(date -Is) done ==="
