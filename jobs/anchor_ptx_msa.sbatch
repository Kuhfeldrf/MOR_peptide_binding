#!/usr/bin/env bash
#SBATCH -J anc_msa
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 180
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/MOR_demo_analysis/logs/anchor_msa_%j.log
#
# Anchor-set arm B, step 2: ONE receptor MSA, shared by all 95 tasks.
#
# NOT OPTIONAL. Arm B v1 ran Protenix MSA-free and predicted the receptor 16.2 A
# from 8F7Q chain R - the 7TM bundle never folded, so there was no orthosteric
# pocket and every score from that run was void. Chai-1 gets away with
# single-sequence only because it supplies ESM embeddings as an MSA surrogate;
# Protenix has no equivalent. See docs/pose_generation_comparison.md.
#
# ONE CALL, NOT 95: the receptor sequence is identical in every task.
#
# THE 19 ANCHOR PEPTIDES GET NO MSA, deliberately, and for a second reason here
# beyond the usual one. The usual one: a 3-31mer carries no usable
# co-evolutionary signal. The second: these are the lab's ordered synthesis
# candidates, and leaving them single-sequence means no anchor sequence is
# transmitted off ORCA.
#
# EXTERNAL DATA TRANSMISSION: --msa_server_mode protenix sends the RECEPTOR
# SEQUENCE ONLY to ByteDance's MSA service. The receptor is human OPRM1 from the
# deposited structure 6DDF - public. Authorised by the user 2026-09-28 on
# exactly these terms; same receptor, same terms. No peptide leaves the cluster.
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

D=$REPO/results/02_cofold_ptx_msa
rm -rf "$D"; mkdir -p "$D"; cd "$D"

python - "$REPO" <<'PY'
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
  --max-time 30 https://protenix-server.com || echo "  curl probe failed (non-fatal)"

echo; echo "=== generating receptor MSA ==="
t0=$SECONDS
protenix msa -i receptor.fasta -o "$D/msa" --msa_server_mode protenix
rc=$?
echo "protenix msa exit=$rc  wall=$((SECONDS-t0))s"
if [ $rc -ne 0 ]; then
  echo "FATAL: MSA generation failed. Fall back to a local MSA pipeline rather"
  echo "       than proceeding MSA-free, which is what voided arm B v1."
  exit $rc
fi

echo; echo "=== what was produced ==="
find "$D/msa" -type f | head -20
for f in $(find "$D/msa" -name '*.a3m' | head -5); do
  echo "  $(basename $f): $(grep -c '^>' "$f") sequences"
done
echo "=== $(date -Is) done ==="
