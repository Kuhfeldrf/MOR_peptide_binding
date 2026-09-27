#!/usr/bin/env bash
# One GROMACS stage: grompp then mdrun.
#
# Usage:
#   04_gmx_run.sh <name> <mdp> <conf> <ref> <top> <ndx> <outdir> [cpt]
#
# The GROMACS prefix comes from config rather than being hardcoded, because it
# is a spack build whose path carries a spec hash (docs/provenance.md).
#
# -update gpu is deliberately NOT used: OPC water carries virtual sites, which
# that code path does not support. Nonbonded and PME still run on the GPU.
set -euo pipefail

NAME=$1; MDP=$2; CONF=$3; REF=$4; TOP=$5; NDX=$6; OUTDIR=$7; CPT=${8:-}

: "${GMX:=$(command -v gmx || true)}"
if [ -z "$GMX" ] && [ -n "${GMX_PREFIX:-}" ]; then
  GMX="$GMX_PREFIX/bin/gmx"
fi
if [ -z "$GMX" ]; then
  echo "FATAL: gmx not found. Set GMX or GMX_PREFIX." >&2
  exit 1
fi

mkdir -p "$OUTDIR"
CONF=$(realpath "$CONF"); REF=$(realpath "$REF")
TOP=$(realpath "$TOP");   NDX=$(realpath "$NDX"); MDP=$(realpath "$MDP")
[ -n "$CPT" ] && CPT=$(realpath "$CPT")
cd "$OUTDIR"

echo "=== grompp ${NAME} ==="
if [ -n "$CPT" ]; then
  "$GMX" grompp -f "$MDP" -c "$CONF" -r "$REF" -t "$CPT" \
         -p "$TOP" -n "$NDX" -o "${NAME}.tpr" -maxwarn 5
else
  "$GMX" grompp -f "$MDP" -c "$CONF" -r "$REF" \
         -p "$TOP" -n "$NDX" -o "${NAME}.tpr" -maxwarn 5
fi

# PME cannot run on the GPU for a NON-DYNAMICAL integrator, which is what
# energy minimisation uses (steepest descents). GROMACS refuses outright:
#
#   Cannot compute PME interactions on a GPU, because:
#     PME GPU does not support: Non-dynamical integrator (use md, sd, etc).
#
# The flags were unconditional, so every minimisation failed before doing any
# work. It went unnoticed because the DAMGO system was minimised by hand
# outside this wrapper - the same DAMGO-only blind spot as the DAM
# moleculetype and the ligand naming.
#
# Short-range non-bonded still goes to the GPU; only PME falls back to CPU, and
# minimisation is short enough that this costs nothing worth optimising.
PME_TARGET="gpu"
if [ "$NAME" = "min" ]; then
  PME_TARGET="cpu"
  echo "note: minimisation is non-dynamical, so PME runs on CPU"
fi

echo "=== mdrun ${NAME} ==="
"$GMX" mdrun -deffnm "${NAME}" -ntmpi 1 -ntomp "${OMP_NUM_THREADS:-8}" \
       -nb gpu -pme "$PME_TARGET"

echo "=== ${NAME} summary ==="
grep -E "Steepest Descents|Potential Energy|Maximum force|Performance" \
     "${NAME}.log" 2>/dev/null | tail -4 || true
