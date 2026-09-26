#!/usr/bin/env bash
#SBATCH -J repack200
#SBATCH -p normal
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 1200
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/repack/%A_%a.log
#SBATCH --array=1-7
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"
source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-pilot

PEPS=(met_enkephalin beta_casomorphin_4_bov beta_casomorphin_5_bov \
      beta_casomorphin_5_hum beta_casomorphin_7_bov beta_casomorphin_7_hum \
      neocasomorphin_6)
PEP=${PEPS[$((SLURM_ARRAY_TASK_ID-1))]}
echo "repacking $PEP at nloop_all=200 nloop=40"

bash scripts/03_pack_membrane.sh \
  "results/03_membrane/${PEP}/complex_opm.pdb" \
  "results/03_membrane/${PEP}/membrane_system.pdb" \
  POPC:CHL1 7:3 1.15 18 22.5 0.15 'Na+' 200 40
rc=$?
echo "pack rc=$rc"
echo "=== QC ==="
python3 scripts/03c_membrane_qc.py --dir "results/03_membrane/${PEP}" \
  --expect-chol-frac 0.3 2>&1 | tail -12
