#!/usr/bin/env bash
#SBATCH -J hdk_arr
#SBATCH -p short,normal
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 240
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/haddock_arr/%A_%a.log
#SBATCH --array=1-33%8
#
# Pocket-directed HADDOCK3: 11 peptides x 3 seeds.
#
# Smoke test (192583) confirmed the AIRs bite: cluster_1_model_1 contacted
# 13/14 orthosteric residues, against 0/10 in the prior ranair run's rank-1
# cluster across 20/20 runs. 770 s per run.
#
# 3 seeds per peptide because the prior classification result attributed its
# null to "peptide length and RUN-TO-RUN NOISE" - replicates are required to
# separate those.
set -uo pipefail
R=/scratch/kuhfeldr-Kuhfeld_temp
D=$R/results/11_haddock
cd "$D"
source /home/kuhfeldr/haddock3-env/bin/activate || { echo "FATAL: venv"; exit 1; }
export TMPDIR=$R/tmp; mkdir -p "$TMPDIR"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" runlist.txt)
RID=$(echo "$LINE" | awk '{print $1}')
CFG=$(echo "$LINE" | awk '{print $2}')
[ -z "${RID:-}" ] && { echo "FATAL: no task on line $SLURM_ARRAY_TASK_ID"; exit 1; }

# Skip work already done (the smoke test did one of these).
if [ -d "runs/$RID/6_seletopclusts" ]; then
  echo "$RID already complete, skipping"; exit 0
fi
rm -rf "runs/$RID"
echo "=== $RID ==="
t0=$SECONDS
haddock3 "$CFG"
rc=$?
echo "haddock3 exit=$rc wall=$((SECONDS-t0))s"
[ $rc -ne 0 ] && { echo "FATAL: haddock3 failed for $RID"; exit $rc; }
[ -d "runs/$RID/6_seletopclusts" ] || { echo "FATAL: no seletopclusts output"; exit 1; }
echo "=== $(date -Is) done ==="
