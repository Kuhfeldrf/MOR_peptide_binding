#!/usr/bin/env bash
#SBATCH -J ptxrnd
#SBATCH -p short
#SBATCH --gres=gpu:a30:1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 120
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/ptxrnd/%A_%a.log
#SBATCH --array=1-25%25
#
# Arm B v2, the FAIR baseline test: 20 real peptides vs 20 length-matched
# scrambles, 5 seeds each.
#
# WHY THIS SET AND NOT THE LENGTH CONTROLS. Arm B v2 scored AUROC 0.815 on the
# length controls, but that set cannot answer whether the method beats the free
# baseline: its positives begin YGGFM and its negatives begin with anything but
# Tyr BY CONSTRUCTION, so "does the sequence start with Tyr?" scores AUROC 1.000
# there. On this set the rule is beatable but strong - 17/20 real peptides start
# with Y, but so do 2/20 scrambles, and 3 real peptides do not, giving the rule
# AUROC 0.875. A method has to clear 0.875 to have earned its GPU time.
#
# The set is EXACTLY length-matched (identical length histograms for both
# classes: 4x4aa, 4x5aa, 6x6aa, 4x7aa, 1x8aa, 1x10aa), which matters because
# arm B's ipTM is severely length-confounded - rho=-0.879 within the length
# controls' negatives. Matched lengths remove that confound from this comparison
# rather than relying on it being small.
#
# IDENTICAL PEPTIDE SET TO ARM C. ADCP ran these same 40 (AUROC 0.562), so this
# gives a direct three-way comparison - ADCP vs Protenix vs the Tyr rule - on
# one peptide set rather than across differently-composed sets.
#
# CAVEAT THAT THIS RUN DOES NOT ADDRESS: the scrambles are assumed non-binders,
# not measured ones. Scrambling YPFP to PYPF keeps the tyrosine, the free amine
# and the identical composition, so a scramble may be a weak binder rather than
# a non-binder. See length_controls_result.md - that ambiguity is exactly why
# the length controls were built. Neither set is sufficient alone.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"

source /etc/profile.d/z00_lmod.sh
module load cuda/12.9.0-none-none || { echo "FATAL: cuda module"; exit 1; }
source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-protenix || { echo "FATAL: conda activate"; exit 1; }

export PROTENIX_ROOT_DIR=$ROOT/.cache/protenix
export TORCH_HOME=$ROOT/.cache/torch

python - <<'PYCHK' || { echo "FATAL: fused kernel .so absent"; exit 1; }
import os, pathlib, sys
import protenix.model.layer_norm as m
so = pathlib.Path(os.path.dirname(m.__file__)) / "fast_layer_norm_cuda_v2.so"
sys.exit(0 if so.is_file() else 1)
PYCHK

# Same receptor MSA as the length-control run: the receptor is identical in every
# task, so it is generated once (jobs/ptx_msa.sbatch) and shared.
MSA=$ROOT/results/02_cofold_ptx_msa/msa/0
for f in "$MSA/pairing.a3m" "$MSA/non_pairing.a3m"; do
  [ -s "$f" ] || { echo "FATAL: missing receptor MSA $f - refusing to run MSA-free, which voided arm B v1"; exit 1; }
done

# 40 peptides x 5 seeds.
N=$(( (SLURM_ARRAY_TASK_ID - 1) / 5 + 1 ))
SEED=$(( (SLURM_ARRAY_TASK_ID - 1) % 5 ))
PEP=$(sed -n "${N}p" data/templates/random5.txt)
[ -z "${PEP:-}" ] && { echo "FATAL: no peptide on line $N of real40.txt"; exit 1; }

echo "=== ${PEP} seed ${SEED} ==="
nvidia-smi --query-gpu=name --format=csv,noheader

python scripts/02_cofold_protenix.py \
  --receptor results/00_receptor/mOR_clean.fasta \
  --library data/templates/random5.tsv \
  --peptide "$PEP" --seed "$SEED" \
  --outdir "$ROOT/results/02_cofold_ptx_rnd5/${PEP}/seed${SEED}" \
  --paired-msa "$MSA/pairing.a3m" --unpaired-msa "$MSA/non_pairing.a3m" \
  --recycles 3 --timesteps 200 --samples 5
rc=$?
echo "exit=$rc"
[ $rc -ne 0 ] && exit $rc
echo "=== $(date -Is) done ==="
