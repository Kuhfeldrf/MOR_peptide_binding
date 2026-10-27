#!/usr/bin/env bash
#SBATCH -J anc_pqr
#SBATCH -p short
#SBATCH -t 60
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -o /scratch/kuhfeldr-Kuhfeld_temp/MOR_demo_analysis/logs/anchor_pqr_env_%j.log
#
# Anchor-set arm B, step 0b: supply `pdb2pqr30`, which step 0 did not.
#
# WHY THIS EXISTS: jobs/anchor_stage0.sbatch (job 242438) failed in 7 s with
#   FileNotFoundError: [Errno 2] No such file or directory: 'pdb2pqr30'
# after extracting chain R and DAMGO correctly. jobs/anchor_ptx_env.sbatch
# asserts "gemmi is all Stage 0 needs"; that is wrong. Stage 0 also SHELLS OUT
# to pdb2pqr30 for the pH 7.4 / PROPKA protonation of D15, and `pip install
# protenix` has no reason to carry it. The claim was never tested, and the
# chain's three dependent jobs (242439 msa, 242440 gate, 242441 array) sat in
# DependencyNeverSatisfied behind it.
#
# WHY A SEPARATE ENV AND NOT `mamba install -n mor-protenix pdb2pqr`:
# mor-protenix is a pip-built environment whose whole point is Protenix 2.0.0
# with the pins its authors ship. pdb2pqr pulls numpy, docutils and mmcif-pdbx
# through the conda solver, which is free to move numpy or torch underneath
# Protenix to satisfy them. Stage 0 needs pdb2pqr30 as an EXECUTABLE, not as an
# importable module in Protenix's interpreter, so there is no reason to accept
# that risk. mor-pqr is built beside it and mor-protenix is not touched.
#
# VERSIONS ARE PINNED TO THE PILOT, not left to the solver:
#   pdb2pqr 3.6.1  (docs/env_solved_main.txt:382)
#   propka  3.5.1  (docs/env_solved_main.txt:237)
#   python  3.12   (environment.yml, the pilot's interpreter)
# Protonation state assignment is PROPKA's output, and a different PROPKA can
# legitimately flip a borderline residue. 6DDF H297 is exactly the kind of
# residue that can flip. Reusing the pilot's versions is what makes the
# rebuilt receptor comparable to the pilot's rather than merely similar.
#
# NOTE ON THE DETERMINISM CLAIM IN STEP 0: anchor_stage0.sbatch says a rebuild
# must reproduce the earlier receptor byte-for-byte and prints sha256 so the
# claim is "checkable rather than asserted". There is NOTHING TO CHECK IT
# AGAINST - no receptor sha256 was recorded in docs/provenance.md before
# /scratch was cleared, and results/ is gitignored. So the sha256 this rebuild
# prints establishes the reference for future rebuilds; it does not confirm a
# match with the pilot. Said here rather than left implied.
set -uo pipefail
REPO=/scratch/kuhfeldr-Kuhfeld_temp/MOR_demo_analysis
CONDA_ROOT=/scratch/kuhfeldr-Kuhfeld_temp/miniforge3
cd "$REPO"
source "$CONDA_ROOT/etc/profile.d/conda.sh"

export PIP_CACHE_DIR=$REPO/.cache/pip
mkdir -p "$PIP_CACHE_DIR" "$REPO/logs"

echo "=== $(date -Is) building mor-pqr ==="
if [ ! -d "$CONDA_ROOT/envs/mor-pqr" ]; then
  mamba create -y -n mor-pqr -c conda-forge \
    python=3.12 pdb2pqr=3.6.1 propka=3.5.1 || exit 1
fi
conda activate mor-pqr || { echo "FATAL: conda activate"; exit 1; }

echo; echo "=== versions ==="
python -V
pdb2pqr30 --version 2>&1 | head -3
python -c "import propka; print('propka', propka.__version__)" 2>&1

# A ONE-BINARY shim directory, not the env's whole bin/ on PATH. Prepending
# $CONDA_ROOT/envs/mor-pqr/bin would also shadow `python`, and Stage 0 runs
# `python scripts/00_receptor.py` - it would execute in mor-pqr, which has no
# gemmi. Only pdb2pqr30 is exposed; the symlink's shebang sends it to mor-pqr's
# own interpreter, so it carries its own propka regardless of the caller's env.
mkdir -p "$REPO/.local/bin"
ln -sf "$CONDA_ROOT/envs/mor-pqr/bin/pdb2pqr30" "$REPO/.local/bin/pdb2pqr30"

echo; echo "=== shim resolves from mor-protenix ==="
conda activate mor-protenix || { echo "FATAL: conda activate mor-protenix"; exit 1; }
export PATH="$REPO/.local/bin:$PATH"
echo "which python:     $(which python)"
echo "which pdb2pqr30:  $(which pdb2pqr30)"
python -c "import gemmi; print('gemmi still OK in mor-protenix:', gemmi.__version__)" || exit 1
pdb2pqr30 --version >/dev/null 2>&1 || { echo "FATAL: shim does not run"; exit 1; }
echo "shim OK - Stage 0 can protonate"

echo; echo "=== $(date -Is) done ==="
