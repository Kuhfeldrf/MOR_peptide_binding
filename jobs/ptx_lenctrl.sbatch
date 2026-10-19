#!/usr/bin/env bash
#SBATCH -J ptxctrl
#SBATCH -p short
#SBATCH --gres=gpu:a30:1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 120
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/ptxctrl/%A_%a.log
#SBATCH --array=1-105%36
#
# Arm B production: the 21 length controls x 5 seeds through Protenix.
# Mirrors controls.sbatch (arm A, Chai-1) task-for-task.
#
# GPU PINNED TO A30 deliberately. Arm A ran on MIXED hardware - its tasks landed
# across orcaga04/05/06 (L40S) and orcaga07-17 (A30) - which leaves GPU
# architecture as an uncontrolled variable in its headline finding, since
# CROSS-SEED AGREEMENT was the only metric that worked there. Pinning one model
# here means arm B's agreement numbers reflect seed, not silicon. A30 is chosen
# over L40S because ORCA has 76 of them against 24 L40S, so the array does not
# queue behind other work. At N_token=291 this model needs a small fraction of
# the A30's 24 GB.
#
# COST, corrected: the smoke run's 1201 s included a one-time data-cache
# download ("Downloading data cache from" in its log). True inference is ~98 s
# per (peptide, seed), so this array is ~2.9 GPU-hours, not the ~35 first
# estimated from the smoke timing.
#
# The fused LayerNorm kernel must already be compiled before this array runs:
# build_directory resolves to site-packages, a SHARED path, so 36 concurrent
# tasks would race ninja building the same extension. It was built by the smoke
# run (191366).
#
# CORRECTION, recorded so the false lesson is not carried forward: an earlier
# version of this file claimed the kernel had to be rebuilt as a "fat binary"
# because TORCH_CUDA_ARCH_LIST=8.9 (Ada/L40S) had been pinned and would exclude
# the A30 (sm_80). That was WRONG. protenix/model/layer_norm/torch_ext_compile.py
# OVERRIDES TORCH_CUDA_ARCH_LIST with its own list derived from the toolkit's
# supported architectures, so the pin was never honoured. The kernel already
# contained sm_70/80/86/89/90/100, and the "rebuilt" .so was byte-identical to
# the original (md5 31df397b...). Setting TORCH_CUDA_ARCH_LIST here is therefore
# a no-op, kept only so the intent is visible.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"

source /etc/profile.d/z00_lmod.sh
module load cuda/12.9.0-none-none || { echo "FATAL: cuda module"; exit 1; }
source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-protenix || { echo "FATAL: conda activate"; exit 1; }

export PROTENIX_ROOT_DIR=$ROOT/.cache/protenix
export TORCH_HOME=$ROOT/.cache/torch
export TORCH_CUDA_ARCH_LIST="8.0;8.9"

# Fail fast if the kernel is not already built, rather than racing 36 JITs.
#
# The check tests for the .so ON DISK. A bare `import fast_layer_norm_cuda_v2`
# does NOT work and is not a valid test: layer_norm.py does
# sys.path.append(current_dir) at import time, so the module is unimportable
# until protenix's own module has been loaded. That false negative failed this
# array's first two tasks while the kernel was present and fine. Importing via
# the protenix path would work but would TRIGGER A COMPILE if absent - exactly
# the race being guarded against - so check the file instead.
python - <<'PYCHK' || { echo "FATAL: fused kernel .so absent; build it in a single serial job first"; exit 1; }
import os, pathlib, sys
import protenix.model.layer_norm as m
so = pathlib.Path(os.path.dirname(m.__file__)) / "fast_layer_norm_cuda_v2.so"
print(f"kernel present: {so}" if so.is_file() else "kernel MISSING")
sys.exit(0 if so.is_file() else 1)
PYCHK

# 21 peptides x 5 seeds, laid out exactly as arm A's controls.sbatch does.
N=$(( (SLURM_ARRAY_TASK_ID - 1) / 5 + 1 ))
SEED=$(( (SLURM_ARRAY_TASK_ID - 1) % 5 ))
PEP=$(sed -n "${N}p" data/templates/controls21.txt)
[ -z "${PEP:-}" ] && { echo "FATAL: no peptide on line $N"; exit 1; }

# Receptor MSA, generated ONCE by jobs/ptx_msa.sbatch (4801 paired / 4178
# unpaired sequences, N_msa 8862 as consumed by the model). Arm B v1 ran without
# it and the receptor never folded - pocket RMSD 11.96-18.34 A vs 1.03 A with it.
MSA=$ROOT/results/02_cofold_ptx_msa/msa/0
for f in "$MSA/pairing.a3m" "$MSA/non_pairing.a3m"; do
  [ -s "$f" ] || { echo "FATAL: missing receptor MSA $f - refusing to run MSA-free"; exit 1; }
done

echo "=== ${PEP} seed ${SEED} ==="
nvidia-smi --query-gpu=name --format=csv,noheader

python scripts/02_cofold_protenix.py \
  --receptor results/00_receptor/mOR_clean.fasta \
  --library results/01_library/length_controls.tsv \
  --peptide "$PEP" --seed "$SEED" \
  --outdir "$ROOT/results/02_cofold_ptx_lenctrl_msa/${PEP}/seed${SEED}" \
  --paired-msa "$MSA/pairing.a3m" --unpaired-msa "$MSA/non_pairing.a3m" \
  --recycles 3 --timesteps 200 --samples 5
rc=$?
echo "exit=$rc"
# Exit codes must not lie - three earlier runs in this experiment reported
# COMPLETED over a dead tool.
[ $rc -ne 0 ] && exit $rc
echo "=== $(date -Is) done ==="
