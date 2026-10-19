#!/usr/bin/env bash
#SBATCH -J ptx_gate
#SBATCH -p short
#SBATCH --gres=gpu:a30:1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 90
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/ptx_gate_%j.log
#
# Arm B v2, step 2: THE GATE. One (peptide, seed) with the receptor MSA, and a
# hard pass/fail on receptor accuracy BEFORE any array is allowed to run.
#
# This exists because arm B v1 did the opposite: 105 tasks ran to completion,
# the analysis produced a clean-looking AUROC of 0.278, and only then was the
# receptor checked - at which point it turned out to be 16.2 A from experiment
# with no orthosteric pocket at all. Every number from that run was void.
#
# PASS CRITERION, fixed in advance: predicted receptor CA RMSD vs 8F7Q chain R
# must be <= 4.0 A. Arm A (Chai-1, single-sequence + ESM) achieves 2.72 A over
# 281 CA, so 4.0 A is a deliberately loose bar - it asks only "is the 7TM bundle
# actually folded", not "is it as good as arm A". v1 scored 16.2 A.
#
# If the gate fails, NO discrimination statistic is computed and the array is
# not launched.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"

source /etc/profile.d/z00_lmod.sh
module load cuda/12.9.0-none-none || { echo "FATAL: cuda module"; exit 1; }
source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-protenix || { echo "FATAL: conda activate"; exit 1; }
export PROTENIX_ROOT_DIR=$ROOT/.cache/protenix
export TORCH_HOME=$ROOT/.cache/torch

MSA=$ROOT/results/02_cofold_ptx_msa/msa/0
PAIRED=$MSA/pairing.a3m
UNPAIRED=$MSA/non_pairing.a3m
for f in "$PAIRED" "$UNPAIRED"; do
  [ -s "$f" ] || { echo "FATAL: missing MSA $f"; exit 1; }
  echo "MSA $(basename $f): $(grep -c '^>' "$f") sequences"
done

PEP=pos__beta_endorphin_1_10
OUT=$ROOT/results/02_cofold_ptx_gate/${PEP}/seed0
rm -rf "$OUT"; mkdir -p "$OUT"

echo; echo "=== inference with receptor MSA ==="
t0=$SECONDS
python scripts/02_cofold_protenix.py \
  --receptor results/00_receptor/mOR_clean.fasta \
  --library results/01_library/length_controls.tsv \
  --peptide "$PEP" --seed 0 \
  --outdir "$OUT" \
  --paired-msa "$PAIRED" --unpaired-msa "$UNPAIRED" \
  --recycles 3 --timesteps 200 --samples 5
rc=$?
echo "exit=$rc  wall=$((SECONDS-t0))s"
[ $rc -ne 0 ] && { echo "FATAL: inference failed"; exit $rc; }

echo; echo "=== N_msa actually used (must be >> 1) ==="
grep -oE "N_msa [0-9]+" "$ROOT/logs/ptx_gate_${SLURM_JOB_ID}.log" | tail -1 || true

echo; echo "=== THE GATE: receptor CA RMSD vs 8F7Q chain R ==="
python - "$ROOT" "$OUT" <<'PY'
import gemmi, numpy as np, pathlib, sys
ROOT, OUT = sys.argv[1], sys.argv[2]

def ca(st, idx):
    pts = []
    for res in st[0][idx]:
        a = res.find_atom("CA", "*")
        if a:
            pts.append([a.pos.x, a.pos.y, a.pos.z])
    return np.asarray(pts)

def kabsch_rmsd(P, Q):
    n = min(len(P), len(Q)); P, Q = P[:n], Q[:n]
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Rot = Vt.T @ np.diag([1, 1, d]) @ U.T
    Pa = (Rot @ (P - pc).T).T + qc
    return float(np.sqrt(np.mean(np.sum((Pa - Q) ** 2, axis=1)))), n

ex = gemmi.read_structure(f"{ROOT}/data/raw/8f7q.cif")
ex.setup_entities(); ex.remove_ligands_and_waters()
exR = ca(ex, [c.name for c in ex[0]].index('R'))

cifs = sorted(pathlib.Path(OUT).glob("pred/*/seed_*/predictions/*_sample_*.cif"))
if not cifs:
    sys.exit("FATAL: no predicted CIF found")
worst = 0.0
for cif in cifs:
    st = gemmi.read_structure(str(cif)); st.setup_entities()
    r, n = kabsch_rmsd(ca(st, 0), exR)
    print(f"  {cif.name}: {r:.2f} A over {n} CA")
    worst = max(worst, r)
best = min(kabsch_rmsd(ca(gemmi.read_structure(str(c)), 0), exR)[0] for c in cifs)
print(f"\n  best {best:.2f} A | worst {worst:.2f} A")
print(f"  arm A baseline 2.72 A | arm B v1 was 16.16 A | gate threshold 4.00 A")
if best <= 4.0:
    print("\n  GATE PASSED - the 7TM bundle is folded; the array may run.")
else:
    print("\n  GATE FAILED - receptor still not folded. Do NOT launch the array")
    print("  and do NOT compute discrimination statistics from this condition.")
    sys.exit(1)
PY
gate=$?
echo "gate exit=$gate"
[ $gate -ne 0 ] && exit $gate
echo "=== $(date -Is) done ==="
