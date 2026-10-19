#!/usr/bin/env bash
#SBATCH -J ptx_msa
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 180
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/ptx_msa_%j.log
#
# Arm B v2, step 1: generate ONE receptor MSA and gate on receptor accuracy.
#
# WHY: arm B v1 is void. Protenix run with --use_msa false predicted the
# receptor 16.2 A from 8F7Q chain R over 281 CA - the 7TM bundle was not folded,
# so there was no pocket for any peptide to bind. Arm A only got away with
# single-sequence because Chai-1 supplies use_esm_embeddings=True as an MSA
# surrogate; Protenix has no equivalent. See docs/pose_generation_comparison.md.
#
# ONE MSA, NOT 105. The receptor sequence is identical in every one of the 21x5
# tasks, so the MSA is generated once here and referenced by path from each
# input JSON. Calling the server 105 times for the same sequence would be both
# wasteful and abusive of a free third-party service.
#
# THE PEPTIDES GET NO MSA, deliberately. They are 10-20mer fragments of milk
# proteins and beta-endorphin; an MSA over a 10-mer carries no usable
# co-evolutionary signal, and leaving them single-sequence keeps this arm
# closer to arm A, where neither chain had an MSA.
#
# EXTERNAL DATA TRANSMISSION: --msa_server_mode protenix sends the RECEPTOR
# SEQUENCE ONLY to ByteDance's MSA service. Explicitly authorised by the user
# 2026-09-28. No peptide sequences and no unpublished ligand data leave ORCA.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"

source /etc/profile.d/z00_lmod.sh
module load cuda/12.9.0-none-none || { echo "FATAL: cuda module"; exit 1; }
source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-protenix || { echo "FATAL: conda activate"; exit 1; }
export PROTENIX_ROOT_DIR=$ROOT/.cache/protenix
export TORCH_HOME=$ROOT/.cache/torch

D=$ROOT/results/02_cofold_ptx_msa
rm -rf "$D"; mkdir -p "$D"; cd "$D"

echo "=== protenix msa --help ==="
protenix msa --help 2>&1 | head -30

# Receptor only, one record.
python - "$ROOT" <<'PY'
import pathlib, sys
root = pathlib.Path(sys.argv[1])
txt = (root / "results/00_receptor/mOR_clean.fasta").read_text().splitlines()
seq = "".join(l.strip() for l in txt if not l.startswith(">"))
out = root / "results/02_cofold_ptx_msa/receptor.fasta"
out.write_text(f">mOR_receptor\n{seq}\n")
print(f"receptor: {len(seq)} aa -> {out}")
PY

echo; echo "=== is the MSA service reachable from a compute node? ==="
curl -s -o /dev/null -w "protenix-server.com HTTP %{http_code} in %{time_total}s\n" \
  --max-time 30 https://protenix-server.com || echo "  curl failed (non-fatal probe)"

echo; echo "=== generating receptor MSA ==="
t0=$SECONDS
protenix msa -i receptor.fasta -o "$D/msa" --msa_server_mode protenix
rc=$?
echo "protenix msa exit=$rc  wall=$((SECONDS-t0))s"
if [ $rc -ne 0 ]; then
  echo "FATAL: MSA generation failed. Fall back to a local MSA pipeline"
  echo "       (docs/colabfold_compatible_msa.md) rather than proceeding MSA-free."
  exit $rc
fi

echo; echo "=== what was produced ==="
find "$D/msa" -type f | head -20
echo; echo "=== MSA depth (lines / sequences per a3m) ==="
for f in $(find "$D/msa" -name '*.a3m' | head -5); do
  echo "  $(basename $f): $(grep -c '^>' "$f") sequences"
done
echo "=== $(date -Is) done ==="
