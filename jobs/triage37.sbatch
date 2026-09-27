#!/usr/bin/env bash
#SBATCH -J triage37
#SBATCH -p normal
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 480
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/triage37/%A_%a.log
#SBATCH --array=1-37%37
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"
source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-pilot

PEP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" data/templates/peptides37.txt)
echo "=== ${PEP} (task ${SLURM_ARRAY_TASK_ID}) ==="
CO="results/02_cofold/${PEP}"
MB="results/03_membrane/${PEP}"

# 1. best pose across the 5 seeds
if [ ! -f "${CO}/best_ligand.pdb" ]; then
  echo "--- best_pose ---"
  python3 scripts/02d_best_pose.py --dir "$CO" --out "${CO}/best_pose.pdb" || exit 1
fi

# 2. receptor + peptide in the OPM frame (rigid transfer from the prediction)
if [ ! -f "${MB}/complex_opm.pdb" ]; then
  echo "--- complex ---"
  mkdir -p "$MB"
  python3 scripts/03a_prep_complex.py \
    --receptor results/00_receptor/mOR_clean.pdb \
    --ligand "${CO}/best_ligand.pdb" \
    --pose "${CO}/best_pose.pdb" --pose-chain A --pose-offset 64 \
    --opm data/raw/6ddf_opm.pdb --outdir "$MB" || exit 1
fi

# 3. GAFF2 / AM1-BCC ligand parameters
if [ ! -f "${MB}/params/ligand_unique.pdb" ]; then
  echo "--- ligand params ---"
  python3 scripts/03b_ligand_params.py \
    --ligand "${CO}/best_ligand.pdb" --outdir "${MB}/params" || exit 1
fi

# 4. MM/GBSA triage score
echo "--- mmgbsa ---"
python3 scripts/05_mmgbsa.py \
  --complex "${MB}/complex_opm.pdb" \
  --mol2 "${MB}/params/ligand.mol2" \
  --frcmod "${MB}/params/ligand.frcmod" \
  --ligand-pdb "${MB}/params/ligand_unique.pdb" \
  --outdir "results/05_mmgbsa/${PEP}" \
  --out "results/05_mmgbsa/${PEP}/score.tsv" \
  --peptide "$PEP" --igb 8 --salt 0.15 --minsteps 2000 --disulfide 140,217
echo "EXIT=$?"
