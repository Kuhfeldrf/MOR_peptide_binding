#!/usr/bin/env bash
#SBATCH -J anc_rec
#SBATCH -p short
#SBATCH -t 30
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/MOR_demo_analysis/logs/anchor_stage0_%j.log
#
# Anchor-set arm B, step 1: regenerate Stage 0 from 6DDF.
#
# results/ is gitignored, so mOR_clean.fasta did not survive in the clone and is
# NOT restored from a copy - it is rebuilt by the same script that built it
# before. Determinism is a hard requirement of that script, so a rebuild must
# reproduce the earlier receptor byte-for-byte; the sha256 below is recorded so
# that claim is checkable rather than asserted.
set -uo pipefail
REPO=/scratch/kuhfeldr-Kuhfeld_temp/MOR_demo_analysis
CONDA_ROOT=/scratch/kuhfeldr-Kuhfeld_temp/miniforge3
cd "$REPO"
source "$CONDA_ROOT/etc/profile.d/conda.sh"
conda activate mor-protenix || { echo "FATAL: conda activate"; exit 1; }

# pdb2pqr30 lives in mor-pqr, not here: Stage 0 protonation needs the binary,
# and installing it into this pip-built env would let the conda solver move
# numpy/torch under Protenix. jobs/anchor_pqr_env.sbatch builds that env and
# this one-binary shim; `python` stays mor-protenix's, which is what has gemmi.
export PATH="$REPO/.local/bin:$PATH"
command -v pdb2pqr30 >/dev/null || {
  echo "FATAL: pdb2pqr30 not on PATH - run jobs/anchor_pqr_env.sbatch first."
  echo "       This is what failed job 242438 after 7 s."; exit 1; }
echo "pdb2pqr30: $(command -v pdb2pqr30)"

python scripts/00_receptor.py --cif data/raw/6ddf.cif --outdir results/00_receptor
rc=$?; echo "exit=$rc"; [ $rc -ne 0 ] && exit $rc

echo "=== receptor ==="
head -2 results/00_receptor/mOR_clean.fasta
sha256sum results/00_receptor/mOR_clean.fasta results/00_receptor/*.pdb 2>/dev/null
echo "=== $(date -Is) done ==="
