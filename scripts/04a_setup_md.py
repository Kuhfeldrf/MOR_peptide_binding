#!/usr/bin/env python3
"""Stage 4a - index groups and position restraints for the MD protocol.

Produces:
  index.ndx    Protein_DAM / MEMB / SOLV for thermostat coupling, and Reduced
               for the trajectory output policy.
  posre_*.itp  position restraints on the solute, released in stages.

THERMOSTAT GROUPS: solute, membrane and solvent are coupled separately so that
a cold solute cannot be hidden by hot solvent in the reported temperature - a
single coupling group is a classic way to miss a solute that never equilibrated.

REDUCED GROUP: receptor + DAMGO + lipids with any atom within a cutoff of the
solute. This is what the build instructions' trajectory policy writes at
analysis frequency; the full system goes out sparsely. Uncontrolled trajectory
writing is named there as the single most likely way to fill the filesystem.
"""
from __future__ import annotations

import pathlib
import re
import sys

import numpy as np
import parmed

D = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp/results/03_membrane_prod")
PROXIMAL = 12.0       # A, lipid atoms this close to the solute are "proximal"

p = parmed.load_file(str(D / "system.parm7"), xyz=str(D / "min.gro"))
xyz = np.array(p.coordinates)
n = len(p.atoms)
print(f"atoms {n}")

LIPID = {"PC", "PA", "OL", "CHL"}
WATER = {"WAT", "HOH"}
ION = {"K+", "Cl-", "NA+"}

names = np.array([a.residue.name for a in p.atoms])
solute = np.array([nm not in LIPID | WATER and nm.strip() not in ION
                   for nm in names])
memb = np.isin(names, list(LIPID))
solv = np.array([nm in WATER or nm.strip() in ION for nm in names])
print(f"  solute {solute.sum()}   membrane {memb.sum()}   solvent {solv.sum()}")
if solute.sum() + memb.sum() + solv.sum() != n:
    sys.exit("FATAL: groups do not partition the system")

# proximal lipids: whole lipid molecules with any atom near the solute
from scipy.spatial import cKDTree
tree = cKDTree(xyz[solute])
lip_idx = np.nonzero(memb)[0]
near = tree.query_ball_point(xyz[lip_idx], PROXIMAL)
hit = {lip_idx[i] for i, nb in enumerate(near) if nb}
# expand to whole residues
res_of = {}
for a in p.atoms:
    res_of.setdefault(a.residue.idx, []).append(a.idx)
hit_res = {p.atoms[i].residue.idx for i in hit}
proximal = sorted(i for r in hit_res for i in res_of[r])
print(f"  proximal lipid atoms within {PROXIMAL} A: {len(proximal)} "
      f"({len(hit_res)} lipid fragments)")

reduced = sorted(set(np.nonzero(solute)[0].tolist()) | set(proximal))
print(f"  Reduced group: {len(reduced)} atoms "
      f"({100*len(reduced)/n:.1f} % of the system)")


def write_group(fh, name, idx):
    fh.write(f"[ {name} ]\n")
    for k in range(0, len(idx), 15):
        fh.write(" ".join(f"{i+1:6d}" for i in idx[k:k+15]) + "\n")
    fh.write("\n")


with (D / "index.ndx").open("w") as fh:
    write_group(fh, "System", list(range(n)))
    write_group(fh, "Protein_DAM", np.nonzero(solute)[0].tolist())
    write_group(fh, "MEMB", np.nonzero(memb)[0].tolist())
    write_group(fh, "SOLV", np.nonzero(solv)[0].tolist())
    write_group(fh, "Reduced", reduced)
print(f"Wrote {D/'index.ndx'}")

# ---------------------------------------------------------------- restraints
# Restraints go inside the moleculetype they belong to, with indices local to
# it. system1 is the receptor (1 molecule); DAM is the ligand.
top = (D / "system.top").read_text()

def moleculetype_atom_count(text: str, name: str) -> int:
    m = re.search(r"\[ moleculetype \]\s*\n;[^\n]*\n\s*" + re.escape(name)
                  + r"\s+\d+\s*\n(.*?)(?=\[ moleculetype \]|\Z)", text, re.S)
    if not m:
        sys.exit(f"FATAL: moleculetype {name} not found")
    block = m.group(1)
    am = re.search(r"\[ atoms \]\s*\n;[^\n]*\n(.*?)(?=\n\s*\[|\Z)", block, re.S)
    if not am:
        sys.exit(f"FATAL: no atoms block in {name}")
    return sum(1 for l in am.group(1).splitlines()
               if l.strip() and not l.strip().startswith(";"))

for mt, fc in (("system1", "POSRES_FC"), ("DAM", "POSRES_FC")):
    cnt = moleculetype_atom_count(top, mt)
    # heavy atoms only: restrain positions of atoms with mass > 2
    heavy = []
    # walk the topology's atom lines for masses
    m = re.search(r"\[ moleculetype \]\s*\n;[^\n]*\n\s*" + re.escape(mt)
                  + r"\s+\d+\s*\n.*?\[ atoms \]\s*\n;[^\n]*\n(.*?)(?=\n\s*\[|\Z)",
                  top, re.S)
    # The counter must advance only on real atom lines: blanks and comments do
    # not occupy an atom index, and counting them pushed restraint indices past
    # the end of the molecule (grompp: "Atom index (4613) out of bounds (1-4610)").
    i = 0
    for l in m.group(1).splitlines():
        if not l.strip() or l.strip().startswith(";"):
            continue
        parts = l.split()
        try:
            mass = float(parts[7])
        except (IndexError, ValueError):
            continue
        i += 1
        if mass > 2.0:
            heavy.append(i)
    out = D / f"posre_{mt}.itp"
    with out.open("w") as fh:
        fh.write(f"; position restraints for {mt}, heavy atoms only\n")
        fh.write("[ position_restraints ]\n")
        fh.write(";  i funct       fcx        fcy        fcz\n")
        for i in heavy:
            fh.write(f"{i:6d}     1  POSRES_FC  POSRES_FC  POSRES_FC\n")
    print(f"Wrote {out.name}: {len(heavy)} of {cnt} atoms restrained")

# insert #include guards at the end of each moleculetype
new = top
for mt in ("system1", "DAM"):
    inc = (f'\n#ifdef POSRES\n#include "posre_{mt}.itp"\n#endif\n')
    if f'posre_{mt}.itp' in new:
        continue
    # place immediately before the NEXT [ moleculetype ] or [ system ]
    pat = re.compile(r"(\[ moleculetype \]\s*\n;[^\n]*\n\s*" + re.escape(mt)
                     + r"\s+\d+\s*\n.*?)(?=\[ moleculetype \]|\[ system \])", re.S)
    m = pat.search(new)
    if not m:
        sys.exit(f"FATAL: could not locate end of moleculetype {mt}")
    new = new[:m.end(1)] + inc + new[m.end(1):]
(D / "system.top").write_text(new)
print("Patched system.top with POSRES includes")
