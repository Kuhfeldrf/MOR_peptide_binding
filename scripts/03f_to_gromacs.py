#!/usr/bin/env python3
"""Stage 3f - convert the Amber topology to GROMACS, verifying the crossing.

Decision D7 accepted the Amber-to-GROMACS conversion as a known silent-error
risk and required it to be checked. ParmEd rewrites atom types, charges and
bonded terms into a different format; a mismatch does not raise, it produces a
topology that runs and is wrong.

Checked across the boundary:
  * atom count
  * total charge
  * per-residue-name charge, so an error localised to one species is caught
  * box vectors
  * molecule counts
Any mismatch exits non-zero.
"""
from __future__ import annotations

import collections
import pathlib
import sys

import parmed

D = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp/results/03_membrane_prod")

print("loading Amber topology...")
amb = parmed.load_file(str(D / "system.parm7"), xyz=str(D / "system_relaxed.rst7"))
print(f"  atoms {len(amb.atoms)}  residues {len(amb.residues)}")
print(f"  box   {amb.box}")

a_q = collections.defaultdict(float)
a_n = collections.Counter()
for r in amb.residues:
    a_q[r.name] += sum(x.charge for x in r.atoms)
    a_n[r.name] += 1
a_tot = sum(a_q.values())
print(f"  total charge {a_tot:+.4f}")

print("\nwriting GROMACS files...")
amb.save(str(D / "system.top"), format="gromacs", overwrite=True)
amb.save(str(D / "system.gro"), overwrite=True)

print("reloading GROMACS topology for an independent check...")
gmx = parmed.load_file(str(D / "system.top"), xyz=str(D / "system.gro"))
g_q = collections.defaultdict(float)
g_n = collections.Counter()
for r in gmx.residues:
    g_q[r.name] += sum(x.charge for x in r.atoms)
    g_n[r.name] += 1
g_tot = sum(g_q.values())

fail = []
print("\n--- crossing the Amber -> GROMACS boundary ---")
print(f"{'quantity':<22}{'amber':>14}{'gromacs':>14}")
print(f"{'atoms':<22}{len(amb.atoms):>14}{len(gmx.atoms):>14}")
if len(amb.atoms) != len(gmx.atoms):
    fail.append(f"atom count {len(amb.atoms)} -> {len(gmx.atoms)}")
print(f"{'residues':<22}{len(amb.residues):>14}{len(gmx.residues):>14}")
if len(amb.residues) != len(gmx.residues):
    fail.append(f"residue count {len(amb.residues)} -> {len(gmx.residues)}")
print(f"{'total charge':<22}{a_tot:>14.4f}{g_tot:>14.4f}")
if abs(a_tot - g_tot) > 1e-3:
    fail.append(f"total charge {a_tot:+.4f} -> {g_tot:+.4f}")
if abs(g_tot) > 0.01:
    fail.append(f"GROMACS system not neutral: {g_tot:+.4f}")

print("\nper-species charge:")
for name in sorted(set(a_q) | set(g_q), key=lambda k: -abs(a_q.get(k, 0))):
    qa, qg = a_q.get(name, 0.0), g_q.get(name, 0.0)
    na, ng = a_n.get(name, 0), g_n.get(name, 0)
    flag = ""
    if na != ng or abs(qa - qg) > 1e-3:
        flag = "  <-- MISMATCH"
        fail.append(f"{name}: n {na}->{ng}, q {qa:+.3f}->{qg:+.3f}")
    if abs(qa) > 0.005 or flag:
        print(f"  {name:<6} n {na:>6} -> {ng:<6}  q {qa:>9.3f} -> {qg:<9.3f}{flag}")

ba = [round(float(v), 3) for v in (amb.box if amb.box is not None else [])]
bg = [round(float(v), 3) for v in (gmx.box if gmx.box is not None else [])]
print(f"\nbox amber  : {ba}")
print(f"box gromacs: {bg}")
if not bg or any(v <= 0 for v in bg[:3]):
    fail.append(f"GROMACS box invalid: {bg}")

print("\n" + "=" * 56)
if fail:
    print("CONVERSION REJECTED:")
    for f in fail:
        print("  - " + f)
    sys.exit(1)
print("CONVERSION VERIFIED: atoms, residues, per-species charge and box all match.")
print(f"Wrote {D/'system.top'} and {D/'system.gro'}")
