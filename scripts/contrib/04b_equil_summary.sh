#!/usr/bin/env bash
# Stage 4b summary: did equilibration actually equilibrate?
#
# Checks the quantities that would reveal a failure:
#   per-group temperature - a cold solute behind hot solvent is the classic
#                           way a single-group thermostat hides a problem
#   pressure              - should sit near 1 bar with large fluctuations
#   box x and z           - should stop drifting
#   area per lipid        - the membrane observable; POPC with 30% cholesterol
#                           is expected near 48-52 A^2
set -e
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pilot
GMX=/scratch/kuhfeldr-Kuhfeld_temp/spack-install/linux-zen4/gromacs-2025.3-45hxd2evh2mxcjqopqa56cr3atcbwlrv/bin/gmx
cd /scratch/kuhfeldr-Kuhfeld_temp/results/03_membrane_prod

# 156 POPC + 66 CHL = 222 lipids, 111 per leaflet
NLEAF=111

for FC in 1000 500 100 10; do
  echo "=========== POSRES_FC = ${FC} (last 250 ps) ==========="
  printf "17\n19\n21\n22\n23\n47\n48\n49\n\n" | \
    "$GMX" energy -f npt_${FC}.edr -b 250 -o ener_${FC}.xvg 2>/dev/null | \
    awk '/^(Temperature|Pressure|Box-X|Box-Y|Box-Z|T-)/ {
           printf "  %-16s %12.3f  +- %8.3f\n", $1, $2, $3 }'
  BX=$(printf "21\n\n" | "$GMX" energy -f npt_${FC}.edr -b 250 -o /dev/null 2>/dev/null | awk '/^Box-X/{print $2}')
  BY=$(printf "22\n\n" | "$GMX" energy -f npt_${FC}.edr -b 250 -o /dev/null 2>/dev/null | awk '/^Box-Y/{print $2}')
  if [ -n "$BX" ] && [ -n "$BY" ]; then
    awk -v x="$BX" -v y="$BY" -v n="$NLEAF" \
      'BEGIN{ printf "  area per lipid   %12.2f  A^2  (box %.2f x %.2f nm)\n", x*y*100/n, x, y }'
  fi
done
