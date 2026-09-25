#!/usr/bin/env bash
# Pack one membrane system with packmol-memgen.
#
# packmol-memgen writes into the current directory, so this cd's into the
# peptide's output directory. Doing that inside the rule's shell command meant
# the log redirect - a relative path - resolved against the new directory after
# the cd and the shell failed before packmol started, leaving an empty log and
# an exit code with no explanation.
#
# --pbc is mandatory. Without it the lipid patch overhangs the box and atoms
# overlap their own periodic images, which produces a system that looks clean
# in Cartesian space while both GROMACS and sander report infinite forces
# (decision D29).
#
# --apl_offset leaves room for the receptor footprint. packmol-memgen's default
# lipid count fills the whole box cross-section, over-packing the bilayer by
# roughly 20 lipids per leaflet (D27). Under-packing is recoverable because
# semi-isotropic NPT relaxes the area per lipid; over-packing is not.
set -euo pipefail

COMPLEX=$1; OUTPDB=$2; LIPIDS=$3; RATIO=$4; APL=$5
DIST=$6; WAT=$7; SALTCON=$8; CATION=$9

OUTDIR=$(dirname "$OUTPDB")
mkdir -p "$OUTDIR"
COMPLEX=$(realpath "$COMPLEX")
cd "$OUTDIR"

echo "=== packing $(basename "$OUTDIR") ==="
echo "lipids ${LIPIDS} ratio ${RATIO} apl_offset ${APL} cation ${CATION}"

packmol-memgen \
  --pdb "$COMPLEX" \
  --lipids "$LIPIDS" --ratio "$RATIO" \
  --apl_offset "$APL" \
  --pbc \
  --preoriented \
  --salt --saltcon "$SALTCON" --salt_c "$CATION" \
  --dist "$DIST" --dist_wat "$WAT" \
  --notprotonate --nottrim --keepligs \
  --output "$(basename "$OUTPDB")" \
  --log packmol_memgen.log
rc=$?
echo "packmol-memgen exit: $rc"

# packmol reports its own convergence, and that verdict is more reliable than
# any geometric proxy computed downstream. It is surfaced here so the QC rule
# does not have to be the first thing to notice a failed pack.
if grep -q "Success" packmol.log 2>/dev/null; then
  echo "PACKMOL: converged"
else
  echo "PACKMOL: did NOT reach tolerance (wrote best solution found)"
  grep -iE "STOP|best solution" packmol.log 2>/dev/null | tail -3 || true
fi
echo "atoms: $(grep -c '^ATOM\|^HETATM' "$(basename "$OUTPDB")" 2>/dev/null || echo 0)"
exit $rc
