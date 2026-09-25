#!/usr/bin/env bash
# Staged restraint release: run NPT once per force constant, each starting from
# the last, then report the diagnostics that would reveal a failed equilibration.
#
# Usage:
#   04_gmx_equilibrate.sh <outdir> <top> <ndx> <mdp> <fc1> <fc2> ...
#
# Releasing restraints in one step lets the solute relax into a lipid and water
# arrangement that has not yet adapted to it, which is the usual way a membrane
# protein ends up distorted during equilibration.
set -euo pipefail

OUTDIR=$1; TOP=$2; NDX=$3; MDP=$4; shift 4
STAGES=("$@")

: "${GMX:=$(command -v gmx || true)}"
if [ -z "$GMX" ] && [ -n "${GMX_PREFIX:-}" ]; then GMX="$GMX_PREFIX/bin/gmx"; fi
if [ -z "$GMX" ]; then echo "FATAL: gmx not found" >&2; exit 1; fi

TOP=$(realpath "$TOP"); NDX=$(realpath "$NDX"); MDP=$(realpath "$MDP")
cd "$OUTDIR"

prev=nvt
for FC in "${STAGES[@]}"; do
  echo "########## NPT, POSRES_FC = ${FC} ##########"
  # grompp reads preprocessor defines from the mdp, not the command line.
  grep -v '^define' "$MDP" > "npt_${FC}.mdp"
  echo "define = -DPOSRES -DPOSRES_FC=${FC}.0" >> "npt_${FC}.mdp"

  "$GMX" grompp -f "npt_${FC}.mdp" -c "${prev}.gro" -r "${prev}.gro" \
         -t "${prev}.cpt" -p "$TOP" -n "$NDX" -o "npt_${FC}.tpr" -maxwarn 5
  "$GMX" mdrun -deffnm "npt_${FC}" -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" \
         -nb gpu -pme gpu
  prev="npt_${FC}"
done

cp "${prev}.gro" npt_final.gro
cp "${prev}.cpt" npt_final.cpt

# Diagnostics. Per-group temperatures matter: a single coupling group can hide
# a solute that never equilibrated behind correctly heated solvent.
{
  echo "Staged restraint release: ${STAGES[*]} kJ/mol/nm^2"
  for FC in "${STAGES[@]}"; do
    echo "--- POSRES_FC = ${FC} (last 250 ps) ---"
    printf "Temperature\nPressure\nBox-X\nBox-Y\nBox-Z\nT-Protein_DAM\nT-MEMB\nT-SOLV\n\n" | \
      "$GMX" energy -f "npt_${FC}.edr" -b 250 -o "ener_${FC}.xvg" 2>/dev/null | \
      grep -E '^(Temperature|Pressure|Box-|T-)' || true
  done
} > equilibration.txt

cat equilibration.txt
