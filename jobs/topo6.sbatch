#!/usr/bin/env bash
#SBATCH -J topo6
#SBATCH -p normal
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 240
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/topo6/%A_%a.log
#SBATCH --array=1-6
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"
source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-pilot

PEPS=(beta_casomorphin_4_bov beta_casomorphin_5_bov beta_casomorphin_5_hum \
      beta_casomorphin_7_bov beta_casomorphin_7_hum neocasomorphin_6)
PEP=${PEPS[$((SLURM_ARRAY_TASK_ID-1))]}
echo "=== topology: $PEP ==="
python3 scripts/03d_build_topology.py --dir "results/03_membrane/${PEP}"
echo "EXIT=$?"
