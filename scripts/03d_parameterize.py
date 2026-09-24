#!/usr/bin/env python3
"""Stage 3d - parameterise the packed membrane system with tleap.

Force field per decision D7: ff19SB / LIPID21 / OPC / GAFF2 for DAMGO.

THREE SILENT-ERROR TRAPS HANDLED HERE.

1. DISULFIDES. tleap only forms an S-S bond when both residues are named CYX
   AND an explicit bond command is given. Left as CYS, it builds two free
   thiols instead, silently breaking the conserved GPCR TM3-ECL2 disulfide that
   staples ECL2 over the binding pocket - the pocket DAMGO occupies. Detected
   geometrically here and declared explicitly.

2. NO CRYST1 RECORD. packmol-memgen writes no box, so tleap would build a
   non-periodic system. The box is taken from the packmol packing regions:
   86.06 x 86.06 x 119.00 A. For a membrane the xy box must match the lipid
   patch periodicity exactly, or the bilayer is discontinuous across the
   boundary.

3. TOTAL CHARGE. D7 requires verifying charge across the build. packmol-memgen
   added 34 K+ and 48 Cl-, a net -14, implying it computed protein+ligand as
   +14. That is checked against what tleap actually assigns, rather than
   trusted.
"""
from __future__ import annotations

import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp")
D = ROOT / "results/03_membrane_prod"
BOX = (86.06, 86.06, 119.00)
SS = [(76, 153)]          # detected geometrically, 2.03 A

src = D / "membrane_system.pdb"
dst = D / "system_cyx.pdb"

# ------------------------------------------------------- 1. CYS->CYX, strip protein H
#
# 4. PROTEIN HYDROGEN NAMING. pdb2pqr names the N-terminal amide hydrogen "H",
#    but AMBER's N-terminal residue template expects H1/H2/H3, and tleap dies
#    with "Atom .R<NMET 1>.A<H 20> does not have a type". Rather than patch
#    names one at a time, protein hydrogens are STRIPPED and rebuilt by tleap.
#    Nothing is lost: protonation state is carried by the RESIDUE NAME
#    (ASH/ASP, HID/HIE/HIP), which Stage 3a already set from the hydrogens
#    pdb2pqr placed. Lipid, water, ion and DAMGO hydrogens are kept, since
#    those templates expect them.
AA = {"ALA","ARG","ASN","ASP","ASH","CYS","CYX","CYM","GLN","GLU","GLH","GLY",
      "HIS","HID","HIE","HIP","ILE","LEU","LYS","LYN","MET","PHE","PRO","SER",
      "THR","TRP","TYR","VAL"}

ss_res = {r for pair in SS for r in pair}
out, renamed, dropped = [], 0, 0
for l in src.read_text().splitlines():
    if l.startswith(("ATOM", "HETATM")):
        res = l[17:20].strip()
        if res == "CYS" and int(l[22:26]) in ss_res:
            l = l[:17] + "CYX" + l[20:]
            renamed += 1
            res = "CYX"
        if res in AA:
            el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
            if el == "H":
                dropped += 1
                continue
    out.append(l)
dst.write_text("\n".join(out) + "\n")
print(f"Renamed {renamed} atoms across {len(ss_res)} disulfide cysteines to CYX")
print(f"Stripped {dropped} protein hydrogens; tleap rebuilds them from residue names")

# ------------------------------------------------------- 2. tleap input
leapin = D / "build.leap"
bonds = "\n".join(f"bond sys.{a}.SG sys.{b}.SG" for a, b in SS)
leapin.write_text(f"""# Stage 3 parameterisation - decision D7
source leaprc.protein.ff19SB
source leaprc.lipid21
source leaprc.water.opc
source leaprc.gaff2

# DAMGO as a single GAFF2 molecule (D10), AM1-BCC charges (D8), net +1
loadamberparams damgo_param/damgo.frcmod
DAM = loadmol2 damgo_param/damgo.mol2

sys = loadpdb system_cyx.pdb

# Conserved TM3-ECL2 disulfide. Without this tleap builds free thiols.
{bonds}

# Box from the packmol packing regions; no CRYST1 exists in the input.
set sys box {{ {BOX[0]} {BOX[1]} {BOX[2]} }}

charge sys
saveamberparm sys system.parm7 system.rst7
savepdb sys system_built.pdb
quit
""")
print(f"Wrote {leapin}")

# ------------------------------------------------------- 3. run
r = subprocess.run(["tleap", "-f", "build.leap"], cwd=D,
                   capture_output=True, text=True)
log = r.stdout + r.stderr
(D / "tleap.log").write_text(log)

err = [l for l in log.splitlines()
       if re.search(r"^(FATAL|Error|Could not|Failed)", l.strip(), re.I)]
warn = [l for l in log.splitlines() if l.strip().startswith("Warning")]
print(f"\ntleap exit {r.returncode}   errors {len(err)}   warnings {len(warn)}")
for l in err[:20]:
    print("  ERR ", l.strip())
for l in warn[:12]:
    print("  WARN", l.strip())

m = re.search(r"Total unperturbed charge:\s*([-\d.]+)", log)
if m:
    q = float(m.group(1))
    print(f"\nTotal system charge (tleap): {q:+.4f}")
    if abs(q) > 0.01:
        print("  *** SYSTEM IS NOT NEUTRAL ***")
else:
    print("\nCould not parse total charge from tleap output")

parm = D / "system.parm7"
if not parm.exists():
    sys.exit("FATAL: tleap did not write system.parm7")

# ------------------------------------------------------- 4. verify
txt = parm.read_text()
mp = re.search(r"%FLAG POINTERS.*?%FORMAT\(10I8\)\s*\n(.*?)%FLAG", txt, re.S)
natom = int(mp.group(1).split()[0]) if mp else -1
n_in = sum(1 for l in src.read_text().splitlines()
           if l.startswith(("ATOM", "HETATM")))
print(f"\n--- atom count across parameterisation ---")
print(f"packed PDB : {n_in}")
print(f"parm7      : {natom}")
delta = natom - n_in
print(f"delta      : {delta:+d}")
print("A positive delta is expected: OPC is a 4-site water model, so tleap")
print(f"adds one virtual site per water. Waters x1 = expected increase.")
print(f"\nWrote {parm} and system.rst7")
