#!/usr/bin/env bash
#SBATCH -J adcp_cal
#SBATCH -p normal
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 480
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/adcp_cal/%A_%a.log
#SBATCH --array=1-9
#
# ADCP SAMPLING CALIBRATION - run BEFORE any production docking.
#
# WHY THIS EXISTS. The first production array (191302) used `-n 100000`.
# ADCP's own default, in CCSBpckgs/ADCP/runADCP.py:106 and its argparse
# default at :431, is numSteps = 2500000 - so that run sampled at 1/25th of
# the value the authors ship, with -N 30 against a default of -N 50.
#
# The consequences were visible in the output and are all one defect:
#   - replica best energies spread by ~300-430 units for 15-20mers, vs ~8 for
#     4mers (rho=+0.937 with length); 34/61 peptides had a spread larger than
#     their own best energy
#   - 268/1830 replicas returned POSITIVE energies, 157 of them above +50,
#     max +417, which are ADCP's 99999/2e6 failure sentinels surfacing as
#     totalE because the search never escaped clashes
#   - replicas "finished" in a median of 5 s
# Those numbers measure sampling failure, not binding. Any discrimination
# statistic computed on them is meaningless and must not be reported.
#
# NOT A DEFECT, recorded so it is not "fixed" later by mistake: every run logs
# "Constraining amino acids: 1". ADCP writes the `constrains` file itself,
# unconditionally, at runADCP.py:182-184, and wires it into its own default
# option string `-p Bias=NULL,external=5,constrains,1.0,1.0` at :226. It is
# normal ADCP behaviour, not a stray file from the dry run.
#
# WHAT THIS JOB DECIDES. Rather than assume 2.5M is enough (or that it is not
# wastefully more than enough), measure convergence directly: 3 peptides
# spanning the length range x 3 step budgets, 30 replicas each. Production
# effort is then set from the LONGEST peptide's requirement and held fixed
# across all lengths - keeping the original, correct decision not to scale
# effort with length, since length is the variable under test.
#
# Acceptance criterion, fixed in advance: a budget is converged for a peptide
# when the replica spread is a small fraction of |best energy| AND no replica
# returns a non-negative energy. Production uses the smallest budget meeting
# that at 20 aa.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
export PATH=$ROOT/tools/adfr/bin:$PATH

# 3 peptides x 3 budgets. Peptides span the range: a 4mer real peptide, and the
# 10mer and 20mer positives (the 20mer is the hardest case in the whole set).
PEPS=(alpha_lactorphin pos__beta_endorphin_1_10 pos__beta_endorphin_1_20)
STEPS=(100000 500000 2500000)

i=$(( (SLURM_ARRAY_TASK_ID - 1) / 3 ))
j=$(( (SLURM_ARRAY_TASK_ID - 1) % 3 ))
PEP=${PEPS[$i]}
N=${STEPS[$j]}

SEQ=$(awk -F'\t' -v p="$PEP" 'NR>1 && $1==p {print $2}' "$ROOT/data/templates/dock_set.tsv")
[ -z "$SEQ" ] && { echo "FATAL: no sequence for $PEP"; exit 1; }
SEQL=$(echo "$SEQ" | tr '[:upper:]' '[:lower:]')

D="$ROOT/results/08_dock/calib/${PEP}_n${N}"
rm -rf "$D"; mkdir -p "$D"; cd "$D"

# ADCP must be given a RELATIVE target path: it does
# os.path.join('./tmp_<job>/', targetFile[:-4], 'translationPoints.npy'), and
# os.path.join discards everything before an absolute path, so an absolute -t
# collapses the unzip directory away. This is what killed array 191241.
ln -sf "$ROOT/results/08_dock/mor_ortho.trg" mor_ortho.trg

echo "=== ${PEP} (${#SEQL} aa) seq=${SEQL}  steps=${N}  replicas=30 ==="
t0=$SECONDS
adcp -t mor_ortho.trg -s "$SEQL" -N 30 -n "$N" -c "$SLURM_CPUS_PER_TASK" -S 42 \
     -o "$PEP" -O 2>&1 | tail -5
rc=${PIPESTATUS[0]}
echo "adcp exit=$rc  wall=$((SECONDS-t0))s"
[ $rc -ne 0 ] && { echo "FATAL: adcp failed"; exit $rc; }

# Per-replica energies, so convergence is measured rather than assumed.
echo "--- per-replica best target energy ---"
for f in ${PEP}_*.out; do
  [ -e "$f" ] || continue
  printf "%s\t%s\n" "$(basename "$f" .out)" \
    "$(grep -m1 'best target energy' "$f" | awk '{print $NF}')"
done | tee energies.tsv

python3 - <<'PY'
import sys
vals = []
for line in open("energies.tsv"):
    parts = line.split("\t")
    if len(parts) == 2:
        try:
            vals.append(float(parts[1]))
        except ValueError:
            pass
if not vals:
    print("SUMMARY no energies parsed"); sys.exit(0)
best, worst = min(vals), max(vals)
bad = sum(1 for v in vals if v >= 0)
print(f"SUMMARY n={len(vals)} best={best:+.3f} worst={worst:+.3f} "
      f"spread={worst-best:.3f} nonneg={bad} "
      f"spread_over_absbest={(worst-best)/abs(best) if best else float('nan'):.2f}")
PY
echo "=== $(date -Is) done ==="
