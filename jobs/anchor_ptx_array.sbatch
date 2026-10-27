#!/usr/bin/env bash
#SBATCH -J anc_ptx
#SBATCH -p short
#SBATCH --gres=gpu:a30:1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 120
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/MOR_demo_analysis/logs/anchor_ptx/%A_%a.log
#SBATCH --array=1-95%36
#
# Anchor-set arm B, step 4: Protenix 2.0.0 over the 20-peptide anchor set.
# 19 unique co-fold inputs x 5 seeds = 95 tasks.
#
# THE SET: data/templates/anchor20.tsv, transcribed from
# Known_Opioid_Peptides_v11.xlsx sheet 'Anchor set (20)'. Four blocks, each one
# document and one assay:
#   A  Wang 2023 Cell Table S2    6 peptides, pKi, radioligand, one lab
#   B  Koch 1985 Table 2          8 peptides, KD uM, 2 species x 4 lengths
#   C  Chiba 1989 Table 1         5 peptides, IC50 uM, antagonists
#   D  Tani 1990 Table 1          1 peptide,  IC50 uM
#
# NO CROSS-BLOCK RANKING IS VALID FROM THIS RUN. pKi, KD and IC50 from four
# papers and four assays are not interconvertible, and this is the same pooling
# error already documented for the 36-peptide affinity set. The comparisons this
# array can support are all WITHIN a block: the Wang pKi series (1.86-955 nM
# plus one measured true negative), and Koch's bovine-vs-human casomorphin
# length ladders - including the #9-vs-#12 pair (YPFPGPI 2 uM vs YPFVEPI 6 uM)
# that Aim 2 rests on.
#
# WHY 19 AND NOT 20: anchors #17 and #18 are casoxin-6 free acid and its methyl
# ester - 17x apart in IC50, identical plain sequence (SRYPSY). Co-folding has
# no way to represent the ester, so running both would spend 5 GPU tasks
# reproducing one answer and would invite reading the pair as resolved. #18 is
# recorded in the library and excluded from the run list.
#
# FOUR ANCHORS ARE SHORTER THAN ANY PEPTIDE THIS PIPELINE HAS CO-FOLDED: #19
# (YPY, 3) and the three 4-mers, one of which is endomorphin-1, the
# highest-affinity anchor in the set. peptides_source.yaml sets min_length: 5
# and would have dropped all four; that filter is not applied here, and the
# 3- and 4-mer results must be read knowing they are extrapolation.
#
# GPU: a30 explicitly, matching jobs/ptx_real40.sbatch, so anchor scores stay
# numerically comparable to the real40 arm B run rather than mixing card types.
set -uo pipefail
REPO=/scratch/kuhfeldr-Kuhfeld_temp/MOR_demo_analysis
CONDA_ROOT=/scratch/kuhfeldr-Kuhfeld_temp/miniforge3
cd "$REPO"

source /etc/profile.d/z00_lmod.sh
module load cuda/12.9.0-none-none || { echo "FATAL: cuda module"; exit 1; }
source "$CONDA_ROOT/etc/profile.d/conda.sh"
conda activate mor-protenix || { echo "FATAL: conda activate"; exit 1; }

export PROTENIX_ROOT_DIR=$REPO/.cache/protenix
export TORCH_HOME=$REPO/.cache/torch

python - <<'PYCHK' || { echo "FATAL: fused kernel .so absent"; exit 1; }
import os, pathlib, sys
import protenix.model.layer_norm as m
so = pathlib.Path(os.path.dirname(m.__file__)) / "fast_layer_norm_cuda_v2.so"
sys.exit(0 if so.is_file() else 1)
PYCHK

MSA=$REPO/results/02_cofold_ptx_msa/msa/0
for f in "$MSA/pairing.a3m" "$MSA/non_pairing.a3m"; do
  [ -s "$f" ] || { echo "FATAL: missing receptor MSA $f - refusing to run MSA-free, which voided arm B v1"; exit 1; }
done

N=$(( (SLURM_ARRAY_TASK_ID - 1) / 5 + 1 ))
SEED=$(( (SLURM_ARRAY_TASK_ID - 1) % 5 ))
PEP=$(sed -n "${N}p" data/templates/anchor_run.txt)
[ -z "${PEP:-}" ] && { echo "FATAL: no peptide on line $N of anchor_run.txt"; exit 1; }

echo "=== ${PEP} seed ${SEED} ==="
nvidia-smi --query-gpu=name --format=csv,noheader

python scripts/02_cofold_protenix.py \
  --receptor results/00_receptor/mOR_clean.fasta \
  --library data/templates/anchor20.tsv \
  --peptide "$PEP" --seed "$SEED" \
  --outdir "$REPO/results/02_cofold_ptx_anchor/${PEP}/seed${SEED}" \
  --paired-msa "$MSA/pairing.a3m" --unpaired-msa "$MSA/non_pairing.a3m" \
  --recycles 3 --timesteps 200 --samples 5
rc=$?
echo "exit=$rc"
[ $rc -ne 0 ] && exit $rc
echo "=== $(date -Is) done ==="
