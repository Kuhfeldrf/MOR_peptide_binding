#!/usr/bin/env python3
"""Give DAMGO unique atom names, then regenerate its GAFF2 parameters.

WHY: DAMGO is treated as ONE residue (D10), but its atoms came from five
source residues (TYR, DAL, GLY, MEA, ETA), each contributing its own N, CA, C
and O. Inside a single residue those names collide. tleap states the
consequence explicitly - "same-name atoms are handled by using the first
occurrence and by ignoring the rest" - so the molecule's connectivity silently
collapsed and tleap died with:

    Atom .R<DAM 1>.A<CZ 31> has force field coordination 4
    but only 3 bonded neighbors.

That error is the lucky case. A subtler name collision could have produced a
topology that built without complaint and was simply wrong.

Atom ORDER is preserved through obabel, antechamber and packmol, so names are
reassigned by index: element symbol plus a per-element counter. The same
mapping is then applied to the DAM atoms already packed into the membrane
system, so no re-packing is needed.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys

ROOT = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp")
P = ROOT / "results/03_membrane_prod/damgo_param"
SYS = ROOT / "results/03_membrane_prod/membrane_system.pdb"

src = P / "damgo_ph74.pdb"
lines = [l for l in src.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
print(f"DAMGO atoms: {len(lines)}")

counts: dict[str, int] = {}
names: list[str] = []
for l in lines:
    el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
    counts[el] = counts.get(el, 0) + 1
    names.append(f"{el}{counts[el]}")

if len(set(names)) != len(names):
    sys.exit("FATAL: generated names are not unique")
print(f"unique names generated: {len(set(names))}  e.g. {names[:6]} ... {names[-3:]}")
print(f"per-element: {counts}")

# ---------------------------------------------------- rewrite the ligand PDB
# Also collapse to a SINGLE residue. The source PDB carried five residue
# numbers from TYR/DAL/GLY/MEA/ETA. antechamber -rn DAM renames them all to
# DAM but leaves them as five DISTINCT residues, and tleap then instantiates
# the ligand once per residue: the built system contained 8 DAM residues
# carrying +7.976 instead of one carrying +1. D10 treats DAMGO as one
# molecule, so it must be one residue.
out = []
for l, n in zip(lines, names):
    l = l[:12] + f"{n:<4}" + l[16:]          # unique atom name
    l = l[:17] + "DAM" + " " + "L" + "   1" + l[26:]   # one residue, DAM 1
    out.append(l)
uniq = P / "damgo_unique.pdb"
uniq.write_text("\n".join(out) + "\nEND\n")
print(f"\nWrote {uniq}")
resids = {l[22:26] for l in out}
print(f"residue ids in ligand PDB: {sorted(resids)}  (must be exactly one)")
if len(resids) != 1:
    sys.exit("FATAL: ligand must be a single residue")

# ---------------------------------------------------- regenerate parameters
print("\n=== antechamber (GAFF2, AM1-BCC, net +1) ===")
r = subprocess.run(
    ["antechamber", "-i", "damgo_unique.pdb", "-fi", "pdb",
     "-o", "damgo.mol2", "-fo", "mol2",
     "-c", "bcc", "-nc", "1", "-at", "gaff2", "-rn", "DAM", "-s", "2", "-pf", "y"],
    cwd=P, capture_output=True, text=True)
print(f"exit {r.returncode}")
if r.returncode != 0:
    print(r.stdout[-2500:]); print(r.stderr[-2500:])
    sys.exit("FATAL: antechamber failed")

r = subprocess.run(["parmchk2", "-i", "damgo.mol2", "-f", "mol2",
                    "-o", "damgo.frcmod", "-s", "gaff2"],
                   cwd=P, capture_output=True, text=True)
print(f"parmchk2 exit {r.returncode}")

mol2 = (P / "damgo.mol2").read_text()
rids = set()
for l in mol2.split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0].strip().splitlines():
    rids.add(l.split()[6])
print(f"mol2 residue ids: {sorted(rids)}")
if len(rids) != 1:
    sys.exit(f"FATAL: mol2 has {len(rids)} residues; tleap would build that "
             f"many copies of the ligand")
block = mol2.split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0].strip().splitlines()
m2names = [l.split()[1] for l in block]
q = sum(float(l.split()[-1]) for l in block)
print(f"\nmol2 atoms {len(block)}   unique names {len(set(m2names))}   charge {q:+.4f}")
if len(set(m2names)) != len(m2names):
    sys.exit("FATAL: mol2 still contains duplicate atom names")
if abs(q - 1.0) > 0.01:
    sys.exit(f"FATAL: charge {q:+.4f}, expected +1")
attn = [l for l in (P / "damgo.frcmod").read_text().splitlines() if "ATTN" in l]
print(f"frcmod ATTN lines: {len(attn)}")

# ---------------------------------------------------- patch the packed system
txt = SYS.read_text().splitlines()
dam_idx = [i for i, l in enumerate(txt)
           if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "DAM"]
print(f"\nDAM atoms in packed system: {len(dam_idx)}")
if len(dam_idx) != len(m2names):
    sys.exit(f"FATAL: {len(dam_idx)} DAM atoms in system vs {len(m2names)} in mol2")

for i, n in zip(dam_idx, m2names):
    txt[i] = txt[i][:12] + f"{n:<4}" + txt[i][16:]
SYS.write_text("\n".join(txt) + "\n")
print(f"Renamed DAM atoms in {SYS.name} to match the mol2, by index")
print("\nDone. Re-run scripts/03d_parameterize.py")
