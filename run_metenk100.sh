#!/usr/bin/env bash
# Launch the met-enkephalin 100 ns pose-stability campaign.
#
# Five systems, one peptide, four independent starting poses plus a scramble:
#
#   metenk_exp   experimental 8F7Q pose      positive reference - sets the
#                                            scale for what "stable" means
#   metenk_hdk   HADDOCK3 blind, emref_67
#   metenk_chai  Chai-1 seed0 model1
#   metenk_ptx   Protenix+MSA seed1 sample0
#   metenk_scr   scrambled met-enkephalin    negative control - tests whether
#                                            staying bound discriminates at all
#
# The two controls are not optional. Every method in this pilot already places
# every peptide in the orthosteric pocket, and MM/GBSA over these same
# complexes scored decoys BETTER than real peptides (AUROC 0.483). Without
# metenk_scr, "the peptide stayed bound" cannot be read as evidence about
# binding; without metenk_exp, the RMSD numbers have no yardstick.
#
# This targets the five prod.xtc files explicitly rather than the `all` rule,
# so the DAG stays at 51 jobs and does not reach into the triage tier.
#
# Same flock discipline as run_pipeline.sh: a second launch is refused rather
# than stealing the Snakemake lock from a live one.
set -uo pipefail

ROOT=/scratch/kuhfeldr-Kuhfeld_temp
cd "$ROOT"

source "$ROOT/miniforge3/etc/profile.d/conda.sh"
conda activate mor-pilot

export GMX_PREFIX=$ROOT/spack-install/linux-zen4/gromacs-2025.3-45hxd2evh2mxcjqopqa56cr3atcbwlrv
export PATH="$GMX_PREFIX/bin:$PATH"

mkdir -p logs

TARGETS=(
  results/04_md/metenk_exp/prod.xtc
  results/04_md/metenk_hdk/prod.xtc
  results/04_md/metenk_chai/prod.xtc
  results/04_md/metenk_ptx/prod.xtc
  results/04_md/metenk_scr/prod.xtc
)

echo "=== $(date -Is) launching met-enkephalin 100 ns campaign ==="
gmx --version 2>/dev/null | grep -m1 "GROMACS version" || echo "WARNING: gmx not on PATH"

exec 9>"$ROOT/.pipeline.lock"
if ! flock -n 9; then
  echo "ERROR: a pipeline launch is already running (pid $(cat "$ROOT/.pipeline.pid" 2>/dev/null || echo '?'))." >&2
  exit 1
fi
echo $$ > "$ROOT/.pipeline.pid"

snakemake -s workflow/Snakefile --cores 1 --unlock >/dev/null 2>&1 || true

snakemake -s workflow/Snakefile --profile config/slurm \
          --rerun-incomplete --keep-going "${TARGETS[@]}"
rc=$?
echo "=== $(date -Is) snakemake exited $rc ==="
exit $rc
