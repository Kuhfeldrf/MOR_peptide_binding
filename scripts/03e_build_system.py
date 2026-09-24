#!/usr/bin/env python3
"""Stage 3e - build the parameterised system, with DAMGO loaded separately.

WHY NOT loadpdb ALONE: loading a PDB that contains DAMGO made tleap instantiate
the DAM unit EIGHT times - 8 x 73 atoms, total charge +7.976 instead of +1.
The membrane is therefore loaded WITHOUT the ligand, DAMGO is loaded from its
own file, and the two are combined. That makes the ligand count structural
rather than something to hope for, and it is asserted below.

NEUTRALISATION: packmol-memgen neutralised the protein (+14) but not DAMGO
(+1), so the system carries a net +1. One water molecule in the bulk, far from
both protein and bilayer, is converted to a chloride. Done explicitly rather
than with addIons so the placement is deterministic and inspectable.
"""
from __future__ import annotations

import pathlib
import re
import subprocess
import sys

import numpy as np

D = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp/results/03_membrane_prod")
BOX = (86.06, 86.06, 119.00)
SS = [(76, 153)]

src = (D / "system_cyx.pdb").read_text().splitlines()
atoms = [l for l in src if l.startswith(("ATOM", "HETATM"))]

# TER RECORDS MUST BE PRESERVED. An earlier version filtered to ATOM/HETATM
# only, dropping all 18,214 TER records; tleap then chained the last protein
# residue straight into the first lipid fragment and died in chirality.c with
# "Atom named C11 from PC did not match". Molecule boundaries are structural
# information, not formatting.
dam = [l for l in atoms if l[17:20].strip() == "DAM"]
rest = [l for l in src
        if (l.startswith("TER")
            or (l.startswith(("ATOM", "HETATM")) and l[17:20].strip() != "DAM"))]
print(f"DAMGO atoms {len(dam)}   everything else {len(rest)}")
if len(dam) != 73:
    sys.exit(f"FATAL: expected 73 DAMGO atoms, found {len(dam)}")

# ---------------------------------------------------- neutralising chloride
prot_xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in rest if l.startswith(("ATOM", "HETATM"))
                     if l[17:20].strip() not in
                     ("WAT", "HOH", "PC", "PA", "OL", "CHL", "K+", "Cl-")])
wat = {}
for l in rest:
    if not l.startswith(("ATOM", "HETATM")):
        continue
    if l[17:20].strip() in ("WAT", "HOH") and l[12:16].strip() in ("O", "OW"):
        key = (l[21], l[22:27])
        wat[key] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])

best, best_score = None, -1e9
for key, p in wat.items():
    if abs(p[2]) < 30:                      # keep clear of the bilayer
        continue
    dmin = float(np.min(np.linalg.norm(prot_xyz - p, axis=1)))
    if dmin > best_score:
        best_score, best = dmin, key
if best is None:
    sys.exit("FATAL: no suitable bulk water found for the counterion")
print(f"Converting water {best} to Cl-; {best_score:.1f} A from nearest protein atom")

out, removed = [], 0
for l in rest:
    if not l.startswith(("ATOM", "HETATM")):
        out.append(l)
        continue
    if (l[21], l[22:27]) == best and l[17:20].strip() in ("WAT", "HOH"):
        if l[12:16].strip() in ("O", "OW"):
            # Columns matter: atom name occupies [12:16], residue name [17:20].
            out.append(l[:12] + " Cl- Cl-" + l[20:])
        removed += 1
        continue
    out.append(l)
print(f"removed {removed} water atoms, added 1 chloride")

(D / "system_nodam.pdb").write_text("\n".join(out) + "\nEND\n")
(D / "damgo_frame.pdb").write_text("\n".join(dam) + "\nTER\nEND\n")
print(f"TER records preserved: {sum(1 for l in out if l.startswith('TER'))}")

# ---------------------------------------------------- leap
bonds = "\n".join(f"bond memb.{a}.SG memb.{b}.SG" for a, b in SS)
(D / "build.leap").write_text(f"""# Stage 3 build - decision D7
source leaprc.protein.ff19SB
source leaprc.lipid21
source leaprc.water.opc
source leaprc.gaff2

loadamberparams damgo_param/damgo.frcmod
DAM = loadmol2 damgo_param/damgo.mol2

memb = loadpdb system_nodam.pdb
{bonds}

lig = loadpdb damgo_frame.pdb

sys = combine {{ memb lig }}
set sys box {{ {BOX[0]} {BOX[1]} {BOX[2]} }}
charge sys
saveamberparm sys system.parm7 system.rst7
savepdb sys system_built.pdb
quit
""")

r = subprocess.run(["tleap", "-f", "build.leap"], cwd=D,
                   capture_output=True, text=True)
log = r.stdout + r.stderr
(D / "tleap.log").write_text(log)
errs = [l for l in log.splitlines() if re.match(r"\s*(FATAL|Error)", l, re.I)]
print(f"\ntleap exit {r.returncode}   errors {len(errs)}")
for l in errs[:10]:
    print("  ", l.strip())

m = re.search(r"Total unperturbed charge:\s*([-\d.]+)", log)
q = float(m.group(1)) if m else None
print(f"Total system charge: {q:+.4f}" if q is not None else "charge not parsed")

if not (D / "system.parm7").exists():
    sys.exit("FATAL: no system.parm7 written")

# ---------------------------------------------------- verify
sys.path.insert(0, "")
import parmed
p = parmed.load_file(str(D / "system.parm7"))
ndam = sum(1 for res in p.residues if res.name == "DAM")
qtot = sum(a.charge for a in p.atoms)
print("\n--- verification ---")
print(f"atoms in topology : {len(p.atoms)}")
print(f"DAM residues      : {ndam}")
print(f"total charge      : {qtot:+.4f}")
ok = True
if ndam != 1:
    print(f"  FAIL: expected exactly 1 DAMGO, found {ndam}"); ok = False
if abs(qtot) > 0.01:
    print(f"  FAIL: system is not neutral ({qtot:+.4f})"); ok = False
print("PASS: one DAMGO, system neutral" if ok else "BUILD REJECTED")
sys.exit(0 if ok else 1)
