#!/usr/bin/env bash
#SBATCH -J anc_gate
#SBATCH -p short
#SBATCH --gres=gpu:a30:1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 120
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/MOR_demo_analysis/logs/anchor_gate_%j.log
#
# Anchor-set arm B, step 3: THE GATE. One (peptide, seed) with the receptor MSA,
# and a hard pass/fail on receptor accuracy BEFORE the 95-task array is allowed
# to run.
#
# This exists because arm B v1 did the opposite: 105 tasks ran to completion, the
# analysis produced a clean-looking AUROC, and only then was the receptor checked
# - at which point it was 16.2 A from experiment with no orthosteric pocket.
#
# PASS CRITERION, fixed in advance: predicted receptor CA RMSD vs 8F7Q chain R
# <= 4.0 A. Arm A (Chai-1) achieves 2.72 A over 281 CA, so 4.0 A asks only "is
# the 7TM bundle folded", not "is it as good as arm A". v1 scored 16.16 A.
#
# GATE PEPTIDE: beta_endorphin_1_31, anchor #4 - 8F7Q's own ligand. Chosen
# because the gate measures the RECEPTOR only and this makes the comparison to
# 8F7Q as direct as it can be. Its result is written straight into the array's
# own output tree, so array task 16 finds it complete and skips it; the gate is
# task 1 of 95, not an extra GPU-hour.
#
# NOTE ON CIRCULARITY, which applies to the SCORE and not to this gate: anchor
# #4 is 8F7Q's deposited ligand, so its co-folding score is a self-consistency
# ceiling and is never evidence of generalisation. The default checkpoint
# protects the structural comparison - protenix_base_default_v1.0.0 has a
# 2021-09-30 training cutoff and 8F7Q was released 2022-12-14, so it is OUT of
# training. Changing --model to the 20250630 checkpoint would put it IN and
# would void both the gate and anchor #4.
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

MSA=$REPO/results/02_cofold_ptx_msa/msa/0
PAIRED=$MSA/pairing.a3m
UNPAIRED=$MSA/non_pairing.a3m
for f in "$PAIRED" "$UNPAIRED"; do
  [ -s "$f" ] || { echo "FATAL: missing MSA $f"; exit 1; }
  echo "MSA $(basename $f): $(grep -c '^>' "$f") sequences"
done

PEP=beta_endorphin_1_31
OUT=$REPO/results/02_cofold_ptx_anchor/${PEP}/seed0
rm -rf "$OUT"; mkdir -p "$OUT"

echo; echo "=== inference with receptor MSA ==="
t0=$SECONDS
python scripts/02_cofold_protenix.py \
  --receptor results/00_receptor/mOR_clean.fasta \
  --library data/templates/anchor20.tsv \
  --peptide "$PEP" --seed 0 \
  --outdir "$OUT" \
  --paired-msa "$PAIRED" --unpaired-msa "$UNPAIRED" \
  --recycles 3 --timesteps 200 --samples 5
rc=$?
echo "exit=$rc  wall=$((SECONDS-t0))s"
[ $rc -ne 0 ] && { echo "FATAL: inference failed"; exit $rc; }

echo; echo "=== N_msa actually used (must be >> 1) ==="
grep -oE "N_msa [0-9]+" "$REPO/logs/anchor_gate_${SLURM_JOB_ID}.log" | tail -1 || true

echo; echo "=== THE GATE: receptor CA RMSD vs 8F7Q chain R ==="
# The calculation used to live inline here and paired CA atoms BY LIST INDEX.
# Protenix renumbers its output 1-281 while 8F7Q chain R runs 66-352, so that
# pairing was off by one residue and reported 4.47 A on job 242496 - a FAILED
# gate on structures that are 2.59 A once the residue correspondence is right
# (identity at the old offset: 11/281; at the correct one: 276/281).
#
# THE THRESHOLD IS UNCHANGED at 4.00 A. What changed is that the metric now
# compares each residue to itself. scripts/gate_rmsd.py establishes the
# correspondence first and REFUSES TO PRINT AN RMSD if no offset reaches 90%
# residue identity, so a misalignment can never again be read as a folding
# failure. It is not 100% because 6DDF and 8F7Q differ at 5 of 281 positions.
python scripts/gate_rmsd.py \
  --pred-dir "$OUT" \
  --ref-cif data/raw/8f7q.cif --ref-chain R \
  --threshold 4.0 --baseline 2.72
gate=$?
echo "gate exit=$gate"
[ $gate -ne 0 ] && exit $gate
echo "=== $(date -Is) done ==="
