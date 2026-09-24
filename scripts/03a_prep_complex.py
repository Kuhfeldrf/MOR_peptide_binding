#!/usr/bin/env python3
"""Stage 3 prep - build the receptor+DAMGO complex with AMBER-correct naming.

TWO SILENT-ERROR TRAPS ARE HANDLED HERE.

1. tleap infers histidine protonation FROM THE RESIDUE NAME (HID/HIE/HIP).
   pdb2pqr's --pdb-output writes the generic name HIS regardless of the state
   it actually assigned. Feeding that to tleap would silently default every
   histidine to HIE and discard the pH 7.4 assignment Stage 0 made - including
   H297, which Stage 0 determined is HID. The same applies to ASP/ASH and
   GLU/GLH. Residues are therefore RENAMED here from the hydrogens present.

2. DAMGO must carry the residue name the GAFF2 parameters were built for (DAM)
   and keep its experimental coordinates.

Emits the complex plus a report of every rename, so the change is auditable.
"""
from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp")
OUT = ROOT / "results/03_membrane"
OUT.mkdir(parents=True, exist_ok=True)


def read_atoms(p: pathlib.Path) -> list[str]:
    return [l for l in p.read_text().splitlines()
            if l.startswith(("ATOM", "HETATM"))]


def residues(lines: list[str]) -> dict[tuple, list[str]]:
    out: dict[tuple, list[str]] = {}
    for l in lines:
        key = (l[21], int(l[22:26]))
        out.setdefault(key, []).append(l)
    return out


# ---------------------------------------------------------------- receptor
rec_lines = read_atoms(ROOT / "results/00_receptor/mOR_clean.pdb")
res = residues(rec_lines)
renames: list[str] = []
fixed: list[str] = []

for key in sorted(res, key=lambda k: k[1]):
    block = res[key]
    name = block[0][17:20].strip()
    anames = {l[12:16].strip() for l in block}
    new = name

    if name in ("HIS", "HID", "HIE", "HIP"):
        hd1, he2 = "HD1" in anames, "HE2" in anames
        if hd1 and he2:
            new = "HIP"
        elif hd1:
            new = "HID"
        elif he2:
            new = "HIE"
        else:
            sys.exit(f"FATAL: His{key[1]} has neither HD1 nor HE2")
    elif name in ("ASP", "ASH"):
        new = "ASH" if "HD2" in anames else "ASP"
    elif name in ("GLU", "GLH"):
        new = "GLH" if "HE2" in anames else "GLU"
    elif name in ("CYS", "CYX", "CYM"):
        # Leave cysteine alone: disulfide detection is tleap's job and
        # guessing here would be worse than letting it decide.
        new = name

    if new != name:
        renames.append(f"  {name}{key[1]} -> {new}")
    for l in block:
        fixed.append(l[:17] + f"{new:>3}" + l[20:])

print(f"Receptor: {len(res)} residues, {len(rec_lines)} atoms")
print(f"Renamed for AMBER: {len(renames)}")
for r in renames:
    print(r)

# The two residues the instructions name explicitly.
for target in (147, 297):
    hit = [l for l in fixed if int(l[22:26]) == target]
    if hit:
        print(f"  CHECK residue {target}: {hit[0][17:20].strip()}")

# ---------------------------------------------------------------- DAMGO
dam = read_atoms(OUT / "damgo_param/damgo_ph74.pdb")
dam_fixed = []
for i, l in enumerate(dam, start=1):
    l = l[:17] + "DAM" + l[20:]          # residue name the GAFF2 params use
    l = l[:21] + "L" + f"{1:>4}" + l[26:]  # single residue, chain L
    dam_fixed.append(l)
print(f"\nDAMGO: {len(dam_fixed)} atoms, renamed to DAM, chain L")

# ---------------------------------------------------------------- write
def renumber(lines: list[str]) -> list[str]:
    out = []
    for i, l in enumerate(lines, start=1):
        out.append(l[:6] + f"{i:>5}" + l[11:])
    return out


complex_pdb = OUT / "complex_damgo.pdb"
body = renumber(fixed + ["TER"] + dam_fixed if False else fixed) + ["TER"] + \
       renumber(dam_fixed)
complex_pdb.write_text("\n".join(body) + "\nTER\nEND\n")

(OUT / "prep_renames.log").write_text(
    "Residues renamed for AMBER (state read from placed hydrogens):\n"
    + ("\n".join(renames) if renames else "  none") + "\n")

print(f"\nWrote {complex_pdb}")
print(f"  receptor atoms: {len(fixed)}")
print(f"  DAMGO atoms   : {len(dam_fixed)}")
print(f"  total         : {len(fixed) + len(dam_fixed)}")
