#!/usr/bin/env bash
#SBATCH -J hdk_sm
#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 240
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/haddock_smoke_%j.log
#
# ONE pocket-directed HADDOCK3 run before the 33-run array. Every arm in this
# project that skipped a smoke test produced a void result.
#
# What this must confirm: (a) CNS runs, (b) the AIRs are ACCEPTED and actually
# restrain to the pocket (not silently ignored), (c) the ensemble input is read
# as 3 models, (d) runtime, so the array can be sized.
set -uo pipefail
R=/scratch/kuhfeldr-Kuhfeld_temp
D=$R/results/11_haddock
cd "$D"
source /home/kuhfeldr/haddock3-env/bin/activate || { echo "FATAL: venv"; exit 1; }
export TMPDIR=$R/tmp; mkdir -p "$TMPDIR"

RID=op_met_enkephalin__seed917
rm -rf "runs/$RID"
echo "=== $RID  (pocket-directed, ambig.tbl) ==="
grep -c '^assign' ambig.tbl | xargs echo "AIR restraints:"
t0=$SECONDS
haddock3 "configs/${RID}.cfg"
rc=$?
echo "haddock3 exit=$rc  wall=$((SECONDS-t0))s"
[ $rc -ne 0 ] && { echo "FATAL: haddock3 failed"; exit $rc; }

echo; echo "=== run tree ==="
ls "runs/$RID" 2>/dev/null
echo; echo "=== did the AIRs actually bite? contacts with pocket residues in top model ==="
TOP=$(find "runs/$RID" -path "*seletopclusts*" -name "*.pdb" | sort | head -1)
echo "top model: $TOP"
python3 - "$TOP" <<'PY'
import sys, math
POCKET={146,149,150,153,219,220,221,223,238,298,299,302,320,324}
try: pdb=sys.argv[1]
except IndexError: sys.exit("no model")
A=[];B=[]
for l in open(pdb):
    if not l.startswith("ATOM"): continue
    ch=l[21]; r=int(l[22:26]); xyz=(float(l[30:38]),float(l[38:46]),float(l[46:54]))
    (A if ch=="A" else B).append((r,xyz))
if not B: sys.exit("no chain B in top model")
hit=set()
for r,a in A:
    if r not in POCKET: continue
    for _,b in B:
        if sum((a[i]-b[i])**2 for i in range(3)) < 25.0:  # 5 A
            hit.add(r); break
print(f"  peptide atoms: {len(B)}, receptor atoms: {len(A)}")
print(f"  pocket residues contacted (<5 A): {len(hit)}/14  -> {sorted(hit)}")
print("  VERDICT:", "AIRs WORKED - peptide is in the orthosteric pocket"
      if len(hit)>=5 else "AIRs did NOT place it in the pocket - investigate")
PY
echo "=== $(date -Is) done ==="
