#!/usr/bin/env bash
#SBATCH -J lenctrl
#SBATCH -p short
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 120
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/lenctrl/%A_%a.log
#SBATCH --array=1-105%45
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"
export PATH=$ROOT/.snakemake-conda/f563a67fec0910c6cf8922520a645294_/bin:$PATH
P=$ROOT/.snakemake-conda/f563a67fec0910c6cf8922520a645294_/bin/python3

N=$(( (SLURM_ARRAY_TASK_ID - 1) / 5 + 1 ))
SEED=$(( (SLURM_ARRAY_TASK_ID - 1) % 5 ))
PEP=$(sed -n "${N}p" data/templates/controls21.txt)
echo "=== ${PEP} seed ${SEED} ==="

$P scripts/02_cofold.py \
  --receptor results/00_receptor/mOR_clean.fasta \
  --library results/01_library/length_controls.tsv \
  --peptide "$PEP" --seed "$SEED" \
  --outdir "$ROOT/results/02_cofold_lenctrl/${PEP}/seed${SEED}" \
  --recycles 3 --timesteps 200 \
  --template-m8 data/templates/6ddf_receptor.m8
echo "EXIT=$?"
