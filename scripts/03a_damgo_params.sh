#!/usr/bin/env bash
# Stage 3a - GAFF2/AM1-BCC parameters for DAMGO.
#
# DAMGO is treated as a single molecule rather than as residues (decision D10),
# so no N-methyl amino-acid residue library is needed. Charges are AM1-BCC
# (D8); escalation to RESP is documented if Stage 6 convergence is poor.
#
# Net charge is +1: the N-terminal amine (pKa ~9.5) is protonated at pH 7.4,
# and that ammonium forms the canonical salt bridge with D147, which Stage 0
# confirmed is deprotonated. Getting this wrong would delete the defining
# interaction of opioid binding.
set -eo pipefail

D=/scratch/kuhfeldr-Kuhfeld_temp/results/03_membrane/damgo_param
cd "$D"

echo "=== input ==="
python3 -c "
from openbabel import pybel
m = next(pybel.readfile('pdb','damgo_ph74.pdb'))
print('formula      :', m.formula)
print('total charge :', m.charge)
print('atoms        :', len(m.atoms))
assert m.charge == 1, f'expected net charge +1, got {m.charge}'
"

echo
echo "=== antechamber: GAFF2 atom types, AM1-BCC charges, net charge +1 ==="
antechamber -i damgo_ph74.pdb -fi pdb \
            -o damgo.mol2 -fo mol2 \
            -c bcc -nc 1 -at gaff2 -rn DAM -s 2 -pf y
echo "antechamber exit: $?"

echo
echo "=== parmchk2: find parameters GAFF2 does not cover ==="
parmchk2 -i damgo.mol2 -f mol2 -o damgo.frcmod -s gaff2

echo
echo "=== verification ==="
NATOM=$(grep -c "^ *[0-9]" damgo.mol2 2>/dev/null || true)
python3 - <<'PY'
import re, sys, pathlib

mol2 = pathlib.Path("damgo.mol2").read_text()
block = mol2.split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0].strip().splitlines()
charges = [float(l.split()[-1]) for l in block]
total = sum(charges)
print(f"atoms in mol2      : {len(block)}")
print(f"sum of AM1-BCC q   : {total:+.4f}")
if abs(total - 1.0) > 0.01:
    sys.exit(f"FATAL: charges sum to {total:+.4f}, expected +1.000")
print("charge check       : OK (+1)")

frc = pathlib.Path("damgo.frcmod").read_text()
attn = [l for l in frc.splitlines() if "ATTN" in l]
print(f"frcmod ATTN lines  : {len(attn)}")
for l in attn[:15]:
    print("   ", l.strip())
if attn:
    print("\nNOTE: ATTN lines mark parameters parmchk2 had to guess by analogy.")
    print("They are not fatal, but they are the parameters least supported by")
    print("GAFF2 and are the first place to look if Stage 6 misbehaves.")
else:
    print("No ATTN lines: GAFF2 covers every parameter directly.")
PY

echo
echo "=== done ==="
ls -la damgo.mol2 damgo.frcmod
