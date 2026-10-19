#!/usr/bin/env bash
#SBATCH -J fpd_smoke
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 120
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/logs/fpd_smoke_%j.log
#
# Arm D, step 1: prepare the 8F7Q self-docking input and validate the
# FlexPepDock pipeline end to end at small nstruct before scaling.
#
# SCORING: ref2015, standard FlexPepDock. NO membrane score function - decided
# 2026-09-28. Rationale in docs/pose_generation_comparison.md (arm D): the only
# membrane term in franklin2019 is fa_water_to_bilayer, a per-residue transfer
# energy by bilayer depth. The receptor is identical across peptides, and the
# peptide binds in the extracellular vestibule at/above the hydrophobic slab
# boundary, so the term is a near-constant offset that CANCELS in ranking -
# and ranking is all a discrimination statistic consumes. Its real benefit is
# preventing TM-helix distortion during minimization, which is obtained more
# cheaply and more controllably by holding the receptor backbone fixed.
# franklin is deferred to a controlled second arm, run only if this one passes,
# so the comparison tests one variable rather than two.
#
# TARGET: 8F7Q, not 6DDF. 8F7Q is beta-endorphin-bound MOR, ALL CANONICAL
# residues, and beta-endorphin is the parent of this pilot's length-control
# positives. DAMGO (6DDF) carries DAL / MEA / ETA - a D-amino acid, an
# N-methylated residue and a C-terminal alcohol - and MEA and ETA have no valid
# canonical Rosetta params (see docs/damgo_notes.md). 8F7Q is also HELD OUT of
# the Protenix checkpoint's 2021-09-30 training cutoff (released 2022-12-14),
# which makes it the common validation target for arms B and D.
# FLAG NAMES, established from the binary and source rather than guessed:
# the prepack flag is -flexPepDocking:flexpep_prepack. An earlier attempt used
# -flexPepDocking:prepack_only, which does NOT exist - that option belongs to
# carbohydrates:glycan_dock - and Rosetta exited 255 before doing any work.
# -flexPepDocking:min_receptor_bb false is set EXPLICITLY, not left to its
# default, because a fixed receptor backbone is the premise on which membrane
# scoring was judged unnecessary for this arm.
set -uo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
ROS=$ROOT/tools/rosetta
BIN=$ROS/source/bin/FlexPepDocking.linuxgccrelease
D=$ROOT/results/09_fpd
mkdir -p "$D"; cd "$D"

source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-protenix || { echo "FATAL: conda"; exit 1; }

# --- Step 1: extract the native complex -------------------------------------
# Chain R = MOR receptor (287 aa), chain P = beta-endorphin (21 aa,
# YGGFMTSEKSQTPLVTLFKNA). Chains A/B/C are Gi, E/F antibody, M/Q a second copy
# of the same complex in the asymmetric unit - all discarded.
#
# Cholesterol (CLR) on chain R IS DISCARDED for this run. It is lipid-facing,
# not in the orthosteric pocket, and would need custom Rosetta params. Recorded
# as a deviation because decisions.md notes GPCR function depends on
# cholesterol - it matters biologically, just not for peptide pose recovery
# against a fixed receptor backbone.
echo "=== extracting native complex from 8f7q.cif ==="
python - <<'PY'
import gemmi
st = gemmi.read_structure('/scratch/kuhfeldr-Kuhfeld_temp/data/raw/8f7q.cif')
st.setup_entities()
st.remove_alternative_conformations()
st.remove_hydrogens()
st.remove_waters()
st.remove_ligands_and_waters()   # drops CLR too; see note above
m = st[0]
keep = {'R', 'P'}
for ch in [c.name for c in m]:
    if ch not in keep:
        m.remove_chain(ch)
st.setup_entities()
# Rosetta wants the peptide as its own chain; R then P is the order on disk.
st.write_pdb('/scratch/kuhfeldr-Kuhfeld_temp/results/09_fpd/native.pdb')
for ch in m:
    print(f"  kept chain {ch.name}: {len(ch)} residues")
PY
[ -s native.pdb ] || { echo "FATAL: native.pdb not written"; exit 1; }
grep -c '^ATOM' native.pdb | xargs echo "ATOM records:"
echo "chains present:"; awk '/^ATOM/ {print substr($0,22,1)}' native.pdb | sort -u | tr '\n' ' '; echo

# --- Step 2: prepack --------------------------------------------------------
# Stage 0 of the standard protocol: repack both partners separately so the
# refinement does not spend its budget fixing input sidechain clashes, which
# would otherwise be scored as if they were interface energy.
echo; echo "=== prepack ==="
t0=$SECONDS
$BIN -s native.pdb \
     -native native.pdb \
     -flexPepDocking:receptor_chain R \
     -flexPepDocking:peptide_chain P \
     -flexPepDocking:flexpep_prepack \
     -ex1 -ex2aro -use_input_sc \
     -out:prefix ppk_ \
     -overwrite \
     > prepack.log 2>&1
rc=$?
echo "prepack exit=$rc  ${SECONDS-t0}s"
if [ $rc -ne 0 ]; then
  echo "FATAL: prepack failed. Rosetta said:"; tail -25 prepack.log; exit $rc
fi
PPK=$(ls -1 ppk_*.pdb 2>/dev/null | head -1)
[ -z "$PPK" ] && { echo "FATAL: no prepacked PDB produced"; tail -25 prepack.log; exit 1; }
echo "prepacked: $PPK"

# --- Step 3: refinement smoke test -----------------------------------------
# nstruct 10 only. This validates the pipeline and reveals the score-term names
# and RMSD columns in the silent file; it is NOT a result. Receptor backbone is
# NOT flexible: pep_refine optimises peptide backbone, rigid-body orientation
# and interface sidechains, leaving the receptor backbone fixed - which is the
# decision that makes membrane scoring unnecessary here.
echo; echo "=== refinement smoke, nstruct=10, from the native pose ==="
t0=$SECONDS
$BIN -s "$PPK" \
     -native native.pdb \
     -flexPepDocking:receptor_chain R \
     -flexPepDocking:peptide_chain P \
     -flexPepDocking:pep_refine \
     -flexPepDocking:min_receptor_bb false \
     -nstruct 10 \
     -ex1 -ex2aro -use_input_sc \
     -out:file:silent smoke_refine.silent \
     -out:file:scorefile smoke_refine.sc \
     -overwrite \
     > refine.log 2>&1
rc=$?
echo "refine exit=$rc  wall=$((SECONDS-t0))s"
if [ $rc -ne 0 ]; then
  echo "FATAL: refinement failed. Rosetta said:"; tail -30 refine.log; exit $rc
fi

echo; echo "=== scorefile header (column names we will analyse) ==="
head -2 smoke_refine.sc
echo; echo "=== all 10 decoys: score vs RMSD ==="
cat smoke_refine.sc
echo; echo "=== per-decoy wall cost ==="
python3 -c "
import sys
w=$SECONDS
print(f'10 decoys in {w}s => ~{w/10:.1f}s per decoy')
print(f'a 500-decoy production run on 1 core: {500*w/10/3600:.2f} h')
"
echo "=== $(date -Is) done ==="
