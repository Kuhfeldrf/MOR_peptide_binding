#!/usr/bin/env bash
#SBATCH -J adcp_rank
#SBATCH -p normal
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 720
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/adcp_rank/%A_%a.log
#SBATCH --array=1-9
#
# ADCP RANK-STABILITY CHECK at 4x the production budget. Replaces array 191302, whose results are DISCARDED -
# see docs/pose_generation_comparison.md, arm C.
#
# BUDGET: -n 10000000, which is 4x the PRODUCTION budget and 4x ADCP's own
# default of 2500000 (CCSBpckgs/ADCP/runADCP.py:106, argparse default :431).
#
# WHY GO ABOVE THE AUTHORS' DEFAULT. The calibration ladder (191377) showed best
# energy had not plateaued at 2.5M for 10- and 20-mers: the 20mer gained 12.9
# units from 100k->500k and still 11.1 from 500k->2.5M. On a rugged MC landscape
# the absolute best may never plateau, so "has the energy converged" is the
# wrong question. The right one is whether the ORDERING of length-matched
# peptides is stable, because that ordering is all a discrimination statistic
# consumes. If Spearman rho between the 2.5M and 10M rankings of these 9
# peptides is high, the production budget is sufficient for ranking even though
# its absolute energies are still drifting. If it is low, ADCP cannot rank at
# this length at any affordable budget - which is itself the finding.
#
# -N 50 replicas, also ADCP's default. NOTE the calibration ladder used -N 30,
# so calibration best-energies are a CONSERVATIVE LOWER BOUND on what this run
# will find (more replicas can only improve a best-of). N is held fixed across
# all 61 peptides, so every within-experiment comparison here is valid.
#
# EFFORT IS FIXED ACROSS LENGTHS, deliberately (decision D-C1). Length is the
# variable under test; scaling effort with it would let the sampling schedule
# manufacture a length trend. Calibration showed that is exactly what happened
# at 100k.
#
# WARNING carried forward from calibration: raw best-target-energy is NOT
# comparable across peptide lengths. At converged sampling it orders by size
# (20aa -35.33 < 10aa -29.16 < 4aa -22.72). Only per-length comparisons -
# positive vs its 6 length-matched negatives, and real vs length-matched
# scramble - are interpretable. Do not rank the library on raw energy.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
export PATH=$ROOT/tools/adfr/bin:$PATH

# The 10-aa group only: 1 positive + 6 length-matched negatives. Rank stability
# WITHIN a single length is the convergence criterion that matters for
# discrimination - absolute best energy never plateaus on a rugged MC landscape,
# but the ORDERING of length-matched peptides either is stable between budgets
# or the assay cannot rank at all. Compared against the same 7 peptides at
# 2.5M from the production array via Spearman rho.
NAME=$(awk -F'\t' 'NR>1 && $4==10 {print $1}' "$ROOT/data/templates/dock_set.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -z "${NAME:-}" ] && { echo "FATAL: no 10-aa peptide at index $SLURM_ARRAY_TASK_ID"; exit 1; }
IFS=$'\t' read -r _ SEQ CLS LEN < <(awk -F'\t' -v p="$NAME" '$1==p {print}' "$ROOT/data/templates/dock_set.tsv")

# Lower case = coil. ADCP reads case as a secondary-structure hint; all-coil is
# the unbiased choice and avoids imposing a helix on peptides with no reason to
# be helical.
SEQL=$(echo "$SEQ" | tr '[:upper:]' '[:lower:]')

D="$ROOT/results/08_dock/runs_rank10M/${NAME}"
rm -rf "$D"; mkdir -p "$D"; cd "$D"

# Relative -t is mandatory: ADCP does
# os.path.join('./tmp_<job>/', target[:-4], 'translationPoints.npy') and
# os.path.join discards everything before an absolute path. An absolute -t
# killed all 61 tasks of array 191241 while SLURM reported COMPLETED.
ln -sf "$ROOT/results/08_dock/mor_ortho.trg" mor_ortho.trg

echo "=== ${NAME} (${CLS}, ${LEN} aa) seq=${SEQL}  n=10000000 N=50 ==="
t0=$SECONDS
adcp -t mor_ortho.trg -s "$SEQL" -N 50 -n 10000000 -c "$SLURM_CPUS_PER_TASK" \
     -S 42 -o "$NAME" -O 2>&1 | tail -5
rc=${PIPESTATUS[0]}
echo "adcp exit=$rc  wall=$((SECONDS-t0))s"
# Exit codes must not lie: three earlier runs in this experiment reported
# COMPLETED over a dead tool. Fail loudly.
[ $rc -ne 0 ] && { echo "FATAL: adcp failed for $NAME"; exit $rc; }

echo "--- per-replica best target energy ---"
for f in ${NAME}_*.out; do
  [ -e "$f" ] || continue
  printf "%s\t%s\n" "$(basename "$f" .out)" \
    "$(grep -m1 'best target energy' "$f" | awk '{print $NF}')"
done | tee energies.tsv

n_out=$(wc -l < energies.tsv)
[ "$n_out" -eq 0 ] && { echo "FATAL: no energies parsed for $NAME"; exit 1; }

# Sentinel rejection is explicit, never averaged in. ADCP surfaces 99999/2e6
# failure sentinels as totalE when a replica never escapes clashes.
python3 - <<'PY'
vals = []
for line in open("energies.tsv"):
    p = line.rstrip("\n").split("\t")
    if len(p) == 2:
        try:
            vals.append(float(p[1]))
        except ValueError:
            pass
valid = [v for v in vals if v < 0]
if not valid:
    print(f"SUMMARY n={len(vals)} valid=0 ALL REPLICAS SENTINEL")
else:
    print(f"SUMMARY n={len(vals)} valid={len(valid)} "
          f"best={min(valid):+.3f} valid_spread={max(valid)-min(valid):.3f} "
          f"sentinel={len(vals)-len(valid)}")
PY
echo "=== $(date -Is) done ==="
