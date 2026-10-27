#!/usr/bin/env bash
#SBATCH -J anc_ptx2
#SBATCH -p short
#SBATCH --gres=gpu:a30:1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 120
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/MOR_demo_analysis/logs/anchor_ptx/%A_%a.log
#SBATCH --array=1-125%36
#
# Anchor set v2: Protenix 2.0.0 over 25 unique co-fold inputs x 5 seeds.
#
# WHAT CHANGED FROM v1 (array 242508, 95/95 COMPLETED). Six free acids were
# added, all from Chiba 1989 Table 1 - the same table, lab and assay that Block
# C already drew on, which the v11 anchor sheet used but did not exhaust:
#
#   YPYY    casoxin B,   h31-34 / b61-64   binding IC50 400 uM   A2 200 uM
#   YPSY    casoxin 4,   b35-38            binding IC50 140 uM   A2 400 uM
#   YPSYG                b35-39            binding IC50  42 uM   A2 500 uM
#   SYPYY                h30-34            binding IC50 700 uM   A2  16 uM
#   PYPYY                b60-64            binding IC50 450 uM   A2  50 uM
#   LPYPYY               b59-64            binding IC50 820 uM   A2 500 uM
#
# Block C now carries ALL TEN free acids from Chiba Table 1 and is the best
# within-document series in the whole anchor set: two N-terminal extension
# ladders (YPSY -> YPSYG -> YPSYGLN, and YPY -> YPYY -> SYPYY/PYPYY -> LPYPYY)
# over 3-10 residues, 20-fold in binding IC50, one assay. It is a better test
# than Block A, whose pKi range is wider but whose six peptides share no ladder.
#
# THE TWO MEASURES ARE NOT INTERCHANGEABLE. Chiba reports a radioreceptor
# binding IC50 and a GPI A2 (the concentration that doubles DAGO's IC50), and
# they disagree in rank: YPSYG is the strongest binder (42 uM) but a weak
# functional antagonist (A2 500 uM), while SYPYY is a weak binder (700 uM) and
# the strongest functional antagonist in the series (A2 16 uM). Co-folding
# confidence is a structural quantity, so the binding IC50 is the column to
# correlate against; the A2 is carried for interpretation, not for ranking.
#
# YPY (#19, 3 residues) IS RETAINED. It was queried on 2026-10-27 on the
# grounds that no anchor should be shorter than 4 residues. Chiba Table 1
# measures it - h31-33 / b61-63, IC50 >1000 uM, A2 >500 uM - and it is the only
# MEASURED near-inactive inside the antagonist series, so it is the one true
# negative Block C has. It is also 3 residues, shorter than anything this
# pipeline has co-folded, and its result must be read as extrapolation.
#
# RESTARTABLE BY DESIGN: 02_cofold_protenix.py treats a run as complete only
# when scores.tsv AND all 5 per-sample confidence JSONs exist, so the 95 tasks
# already computed by array 242508 skip in ~13 s and only the 30 new ones run.
# Submitting the full 125 keeps one job = one complete set.
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
PEP=$(sed -n "${N}p" data/templates/anchor_run_v2.txt)
[ -z "${PEP:-}" ] && { echo "FATAL: no peptide on line $N of anchor_run_v2.txt"; exit 1; }

echo "=== ${PEP} seed ${SEED} ==="
nvidia-smi --query-gpu=name --format=csv,noheader

python scripts/02_cofold_protenix.py \
  --receptor results/00_receptor/mOR_clean.fasta \
  --library data/templates/anchor_library.tsv \
  --peptide "$PEP" --seed "$SEED" \
  --outdir "$REPO/results/02_cofold_ptx_anchor/${PEP}/seed${SEED}" \
  --paired-msa "$MSA/pairing.a3m" --unpaired-msa "$MSA/non_pairing.a3m" \
  --recycles 3 --timesteps 200 --samples 5
rc=$?
echo "exit=$rc"
[ $rc -ne 0 ] && exit $rc
echo "=== $(date -Is) done ==="
