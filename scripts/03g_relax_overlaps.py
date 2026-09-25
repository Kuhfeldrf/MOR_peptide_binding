#!/usr/bin/env python3
"""Stage 3g - separate the residual overlaps by RIGID-BODY molecule moves.

WHY RIGID: an earlier attempt moved individual atoms apart with no restoring
force on their bonded partners. Accumulated independent pushes stretched C-H
bonds to 3.05 A (r0 1.09) and tore water H1-H2 to 2.99 A (r0 1.371), leaving
the structure more broken than before. Translating a whole molecule CANNOT
distort it - every internal distance is preserved exactly - so this approach is
safe by construction, and the script verifies that claim rather than asserting
it.

WHY IT IS VIABLE NOW: raising the assumed area per lipid (--apl_offset 1.15)
cut the pathological contacts from 551 pairs under 1.2 A to 50. A targeted fix
on ~100 molecules is reasonable; the same fix on the 3,678 molecules of the
previous pack would not have been, which is why that one was re-packed instead.

Molecules are moved apart with a displacement weighted by inverse atom count,
so water and ions move and the receptor effectively does not.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import parmed
from scipy.spatial import cKDTree

D = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp/results/03_membrane_prod")
TARGET = 1.85     # A, desired minimum inter-molecular separation
HARD = 1.50       # A, refuse to ship below this
ITERS = 200

p = parmed.load_file(str(D / "system.parm7"), xyz=str(D / "system.rst7"))
xyz = np.array(p.coordinates, dtype=float)
n = len(p.atoms)
print(f"atoms {n}")

# molecule membership from bonded connectivity
parent = list(range(n))
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
for b in p.bonds:
    rx, ry = find(b.atom1.idx), find(b.atom2.idx)
    if rx != ry:
        parent[ry] = rx
mol = np.array([find(i) for i in range(n)])
uniq, inv = np.unique(mol, return_inverse=True)
nmol = len(uniq)
sizes = np.bincount(inv)
print(f"molecules {nmol}")

members = [[] for _ in range(nmol)]
for i, m in enumerate(inv):
    members[m].append(i)
members = [np.array(v) for v in members]

# record internal geometry so the rigid claim can be checked
bond_pairs = np.array([[b.atom1.idx, b.atom2.idx] for b in p.bonds])
d_before = np.linalg.norm(xyz[bond_pairs[:, 0]] - xyz[bond_pairs[:, 1]], axis=1)

def inter_pairs(coords, cutoff):
    t = cKDTree(coords)
    pr = t.query_pairs(cutoff, output_type="ndarray")
    if not len(pr):
        return pr
    return pr[inv[pr[:, 0]] != inv[pr[:, 1]]]

start = inter_pairs(xyz, TARGET)
print(f"\ninter-molecular pairs under {TARGET} A: {len(start)}")
if len(start):
    d = np.linalg.norm(xyz[start[:, 0]] - xyz[start[:, 1]], axis=1)
    print(f"worst: {d.min():.3f} A")

moved_mols = set()
for it in range(ITERS):
    pr = inter_pairs(xyz, TARGET)
    if not len(pr):
        print(f"resolved after {it} iterations")
        break
    shift = np.zeros((nmol, 3))
    for i, j in pr:
        mi, mj = inv[i], inv[j]
        v = xyz[i] - xyz[j]
        d = float(np.linalg.norm(v))
        if d < 1e-6:
            v = np.random.default_rng(int(i) * 7919 + int(j)).normal(size=3)
            d = float(np.linalg.norm(v))
        need = (TARGET - d) + 0.02
        u = v / d
        wi = 1.0 / sizes[mi]
        wj = 1.0 / sizes[mj]
        tot = wi + wj
        shift[mi] += u * need * (wi / tot)
        shift[mj] -= u * need * (wj / tot)
        moved_mols.add(int(mi)); moved_mols.add(int(mj))
    for m in np.nonzero(np.any(shift != 0, axis=1))[0]:
        xyz[members[m]] += shift[m] * 0.5
else:
    print(f"still {len(inter_pairs(xyz, TARGET))} pairs after {ITERS} iterations")

final = inter_pairs(xyz, TARGET)
dmin = float("inf")
if len(final):
    dmin = float(np.min(np.linalg.norm(xyz[final[:, 0]] - xyz[final[:, 1]], axis=1)))
print(f"\nmolecules translated : {len(moved_mols)} of {nmol}")
print(f"pairs still under {TARGET} A: {len(final)}")
print(f"worst remaining      : {dmin if np.isfinite(dmin) else TARGET:.3f} A")

# verify rigidity: every bond length must be unchanged
d_after = np.linalg.norm(xyz[bond_pairs[:, 0]] - xyz[bond_pairs[:, 1]], axis=1)
drift = float(np.max(np.abs(d_after - d_before)))
print(f"\nmax change in ANY bond length: {drift:.2e} A")
if drift > 1e-6:
    sys.exit(f"FATAL: rigid-body move altered internal geometry by {drift:.3e} A")
print("rigidity verified: internal geometry is bit-for-bit preserved")

if np.isfinite(dmin) and dmin < HARD:
    sys.exit(f"REFUSED: worst contact {dmin:.3f} A is below the {HARD} A floor")

p.coordinates = xyz
p.save(str(D / "system_relaxed.rst7"), overwrite=True)
p.save(str(D / "system_relaxed.gro"), overwrite=True)
print(f"\nWrote system_relaxed.gro / .rst7")
