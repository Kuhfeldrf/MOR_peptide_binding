#!/usr/bin/env bash
#SBATCH -J boltz_sm
#SBATCH -p short
#SBATCH --gres=gpu:a30:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 120
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/boltz_smoke_%j.log
#
# Boltz-2 smoke test: 2 peptides before the full 34. Every prior arm in this
# project that skipped a smoke test produced a void result, so this checks the
# affinity head actually runs on peptide-as-SMILES input and reveals the real
# output schema before any conclusion is drawn.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"
source /etc/profile.d/z00_lmod.sh
module load cuda/12.9.0-none-none || { echo "FATAL: cuda module"; exit 1; }
source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-boltz || { echo "FATAL: conda activate"; exit 1; }
export TORCH_HOME=$ROOT/.cache/torch
export BOLTZ_CACHE=$ROOT/.cache/boltz
mkdir -p "$BOLTZ_CACHE"

MSA=$ROOT/results/02_cofold_ptx_msa/msa/0/non_pairing.a3m
[ -s "$MSA" ] || { echo "FATAL: receptor MSA missing at $MSA"; exit 1; }
echo "receptor MSA: $(grep -c '^>' "$MSA") sequences"

IN=$ROOT/results/10_boltz/smoke_in
rm -rf "$IN" "$ROOT/results/10_boltz/smoke_out"; mkdir -p "$IN"
python scripts/08_boltz_affinity.py \
  --inputs docs/boltz_affinity_inputs.tsv \
  --receptor results/00_receptor/mOR_clean.fasta \
  --msa "$MSA" --outdir "$IN" --limit 2 || exit 1
echo "--- input yaml ---"; head -20 "$IN"/*.yaml | head -24

nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo; echo "=== boltz predict (affinity) ==="
t0=$SECONDS
boltz predict "$IN" \
  --out_dir "$ROOT/results/10_boltz/smoke_out" \
  --cache "$BOLTZ_CACHE" \
  --accelerator gpu --devices 1 \
  --diffusion_samples 1 --sampling_steps 50 \
  --diffusion_samples_affinity 3 \
  --output_format pdb
rc=$?
echo "boltz exit=$rc  wall=$((SECONDS-t0))s"
[ $rc -ne 0 ] && { echo "FATAL: boltz predict failed"; exit $rc; }

echo; echo "=== output tree ==="
find "$ROOT/results/10_boltz/smoke_out" -maxdepth 4 -type f | head -20
echo; echo "=== affinity json, verbatim ==="
for f in $(find "$ROOT/results/10_boltz/smoke_out" -name "affinity*.json" | head -2); do
  echo "--- $f"; cat "$f"; echo
done
echo "=== $(date -Is) done ==="
