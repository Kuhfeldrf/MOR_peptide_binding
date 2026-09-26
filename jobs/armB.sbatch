#!/usr/bin/env bash
#SBATCH -J armB_tmpl
#SBATCH -p short
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 90
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/armB/%A_%a.log
#SBATCH --array=1-185%40
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"
export PATH=$ROOT/.snakemake-conda/f563a67fec0910c6cf8922520a645294_/bin:$PATH
P=$ROOT/.snakemake-conda/f563a67fec0910c6cf8922520a645294_/bin/python3

read -r PEP SEED < <(sed -n "${SLURM_ARRAY_TASK_ID}p" data/templates/tasks.txt)
echo "task=${SLURM_ARRAY_TASK_ID} peptide=${PEP} seed=${SEED}"

$P scripts/02_cofold.py \
  --receptor results/00_receptor/mOR_clean.fasta \
  --library results/01_library/peptides.tsv \
  --peptide "$PEP" --seed "$SEED" \
  --outdir "$ROOT/results/02_cofold_tmpl/${PEP}/seed${SEED}" \
  --recycles 3 --timesteps 200 \
  --template-m8 data/templates/6ddf_receptor.m8
echo "EXIT=$?"
