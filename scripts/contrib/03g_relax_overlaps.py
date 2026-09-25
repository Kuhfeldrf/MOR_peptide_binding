#!/usr/bin/env python3
"""Stage 3g - rigid-body overlap relaxation, under MINIMUM IMAGE.

Two properties, both verified rather than assumed:

RIGID: whole molecules are translated, so internal geometry is preserved
exactly. An earlier version moved individual atoms and stretched C-H bonds to
3.05 A and tore water molecules apart (D28).

PERIODIC: separations are measured under minimum image. The earlier version
used raw Cartesian distances and was blind to atoms overlapping their own
periodic images, which was the governing defect in Stage 3 (D29).
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import parmed
from scipy.spatial import cKDTree

D = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp/results/03_membrane_prod")
TARGET = 1.85
HARD = 1.00      # below this the LJ term overflows; at ~1 A it is large but
                 # FINITE, which is what steepest descent with a small step is
                 # for. An earlier 1.5 A floor refused a 1.058 A structure that
                 # minimisation can relieve.
ITERS = 600

p = parmed.load_file(str(D / "system.parm7"), xyz=str(D / "system.rst7"))
xyz = np.array(p.coordinates, dtype=float)
box = np.array(p.box[:3], dtype=float)
n = len(p.atoms)
print(f"atoms {n}   box {box}")

par = list(range(n))
def find(x):
    while par[x] != x:
        par[x] = par[par[x]]; x = par[x]
    return x
for b in p.bonds:
    rx, ry = find(b.atom1.idx), find(b.atom2.idx)
    if rx != ry: par[ry] = rx
mol = np.array([find(i) for i in range(n)])
uniq, inv = np.unique(mol, return_inverse=True)
nmol = len(uniq)
sizes = np.bincount(inv)
members = [[] for _ in range(nmol)]
for i, m in enumerate(inv):
    members[m].append(i)
members = [np.array(v) for v in members]
print(f"molecules {nmol}")

bond_pairs = np.array([[b.atom1.idx, b.atom2.idx] for b in p.bonds])
d_before = np.linalg.norm(xyz[bond_pairs[:, 0]] - xyz[bond_pairs[:, 1]], axis=1)

origin = xyz.min(0)

def wrapped(c):
    w = c - origin
    return w - np.floor(w / box) * box

def mic(v):
    """minimum-image displacement"""
    return v - np.round(v / box) * box

def close(c, cutoff):
    t = cKDTree(wrapped(c), boxsize=box)
    pr = t.query_pairs(cutoff, output_type="ndarray")
    if not len(pr):
        return pr
    return pr[inv[pr[:, 0]] != inv[pr[:, 1]]]

start = close(xyz, TARGET)
print(f"\nminimum-image inter-molecular pairs under {TARGET} A: {len(start)}")
if len(start):
    d = np.linalg.norm(mic(xyz[start[:, 0]] - xyz[start[:, 1]]), axis=1)
    print(f"worst: {d.min():.3f} A")

moved = set()
for it in range(ITERS):
    pr = close(xyz, TARGET)
    if not len(pr):
        print(f"resolved after {it} iterations")
        break
    shift = np.zeros((nmol, 3))
    for i, j in pr:
        mi, mj = inv[i], inv[j]
        v = mic(xyz[i] - xyz[j])
        d = float(np.linalg.norm(v))
        if d < 1e-6:
            v = np.random.default_rng(int(i) * 7919 + int(j)).normal(size=3)
            d = float(np.linalg.norm(v))
        u = v / d
        need = (TARGET - d) + 0.02
        wi, wj = 1.0 / sizes[mi], 1.0 / sizes[mj]
        tot = wi + wj
        shift[mi] += u * need * (wi / tot)
        shift[mj] -= u * need * (wj / tot)
        moved.add(int(mi)); moved.add(int(mj))
    for m in np.nonzero(np.any(shift != 0, axis=1))[0]:
        xyz[members[m]] += shift[m] * 0.5
else:
    print(f"still {len(close(xyz, TARGET))} pairs after {ITERS} iterations")

final = close(xyz, TARGET)
dmin = float("inf")
if len(final):
    dmin = float(np.min(np.linalg.norm(mic(xyz[final[:, 0]] - xyz[final[:, 1]]), axis=1)))
print(f"\nmolecules translated : {len(moved)} of {nmol}")
print(f"pairs still under {TARGET} A (min image): {len(final)}")
print(f"worst remaining      : {dmin if np.isfinite(dmin) else TARGET:.3f} A")

d_after = np.linalg.norm(xyz[bond_pairs[:, 0]] - xyz[bond_pairs[:, 1]], axis=1)
drift = float(np.max(np.abs(d_after - d_before)))
print(f"max change in ANY bond length: {drift:.2e} A")
if drift > 1e-6:
    sys.exit(f"FATAL: rigid move altered geometry by {drift:.3e} A")
print("rigidity verified")

if np.isfinite(dmin) and dmin < HARD:
    sys.exit(f"REFUSED: worst minimum-image contact {dmin:.3f} A below {HARD} A")

p.coordinates = xyz
p.save(str(D / "system_relaxed.rst7"), overwrite=True)
p.save(str(D / "system_relaxed.gro"), overwrite=True)
print("Wrote system_relaxed.rst7 / .gro")
