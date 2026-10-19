#!/usr/bin/env bash
#SBATCH -J boltz_aff
#SBATCH -p short
#SBATCH --gres=gpu:a30:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 240
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/boltz_full_%A_%a.log
#SBATCH --array=0-1
#
# Boltz-2 affinity on all 34 peptides, run TWICE: with and without
# --affinity_mw_correction.
#
# WHY BOTH. Every score in this project turned out to read peptide SIZE rather
# than chemistry (iLIS rho=-0.834 vs length). Boltz-2 ships an explicit
# molecular-weight correction for its affinity head, so running both ways tests
# directly whether that correction removes the size artefact - the one failure
# mode this project has characterised in depth.
#
# EXPECTATION, recorded before the result: the docs state the affinity head is
# "not recommended ... significantly larger than 56 atoms (limit set during
# training)". Only 17/34 of these peptides are inside that range. So the primary
# analysis is the IN-DISTRIBUTION subset (<=56 heavy atoms, lengths 4-6, still
# 4.00 logs of pKi, length decorrelated rho=-0.094); the full 34 is secondary
# and explicitly off-distribution. A weak result on the off-distribution half
# would be expected behaviour, not a model failure.
#
# Metric is SPEARMAN, not Pearson: Boltz-2 outputs log10(IC50/uM) while the
# labels are pKi, and those differ by an assay-dependent Cheng-Prusoff factor.
# Rank is comparable; absolute value is not.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"
source /etc/profile.d/z00_lmod.sh
module load cuda/12.9.0-none-none || { echo "FATAL: cuda module"; exit 1; }
source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-boltz || { echo "FATAL: conda activate"; exit 1; }
export TORCH_HOME=$ROOT/.cache/torch
export BOLTZ_CACHE=$ROOT/.cache/boltz

MSA=$ROOT/results/02_cofold_ptx_msa/msa/0/non_pairing.a3m
[ -s "$MSA" ] || { echo "FATAL: receptor MSA missing"; exit 1; }

if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
  TAG=nomw; MWFLAG=""
else
  TAG=mw;   MWFLAG="--affinity_mw_correction"
fi
echo "=== variant: $TAG  (mw flag: '${MWFLAG:-none}') ==="

IN=$ROOT/results/10_boltz/full_in
OUT=$ROOT/results/10_boltz/full_out_$TAG
mkdir -p "$IN"; rm -rf "$OUT"
python scripts/08_boltz_affinity.py \
  --inputs docs/boltz_affinity_inputs.tsv \
  --receptor results/00_receptor/mOR_clean.fasta \
  --msa "$MSA" --outdir "$IN" || exit 1
echo "inputs: $(ls "$IN"/*.yaml | wc -l)"

t0=$SECONDS
boltz predict "$IN" \
  --out_dir "$OUT" \
  --cache "$BOLTZ_CACHE" \
  --accelerator gpu --devices 1 \
  --diffusion_samples 1 --sampling_steps 50 \
  --diffusion_samples_affinity 5 \
  $MWFLAG \
  --output_format pdb
rc=$?
echo "boltz exit=$rc  wall=$((SECONDS-t0))s"
[ $rc -ne 0 ] && { echo "FATAL: boltz predict failed"; exit $rc; }

echo "=== collect ==="
python - "$OUT" "$TAG" <<'PY'
import json, pathlib, sys, csv
out=pathlib.Path(sys.argv[1]); tag=sys.argv[2]
rows=[]
for f in sorted(out.rglob("affinity_*.json")):
    d=json.loads(f.read_text())
    rows.append({"name":f.parent.name,
                 "affinity_pred_value":d.get("affinity_pred_value"),
                 "affinity_probability_binary":d.get("affinity_probability_binary")})
p=pathlib.Path(f"/scratch/kuhfeldr-Kuhfeld_temp/docs/boltz_affinity_{tag}.tsv")
with p.open("w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(rows)
print(f"wrote {p} ({len(rows)} peptides)")
PY
echo "=== $(date -Is) done ==="
