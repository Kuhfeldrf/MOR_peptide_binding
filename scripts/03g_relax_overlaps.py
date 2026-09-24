#!/usr/bin/env python3
"""Stage 3g - relax the residual packing overlaps before minimisation.

WHY THIS EXISTS: packmol did not reach its 2.0 A tolerance. Its log ends with
"packing problem with the desired distance tolerance ... contains the best
solution found" and STOP 173 - it wrote the best configuration it managed
rather than a converged one. The result carries 113 NON-BONDED pairs closer
than 0.8 A, all lipid-lipid. At 0.5 A the Lennard-Jones term overflows, which
is why GROMACS reported an infinite force and steepest descent quit after 16
steps at 5.2e17 kJ/mol.

WHAT THIS DOES: separates non-bonded atom pairs closer than a floor by pushing
each atom half the deficit along their separation vector, iterating because a
push can create a new contact. Bond geometry is perturbed slightly; the
following minimisation restores it, which is exactly the work minimisation is
good at once the forces are finite.

This is a documented remedy for an unconverged pack, not a silent repair: the
number of pairs touched and the worst remaining contact are reported, and the
script refuses to write a structure that still contains an overflow-level
contact.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import parmed
from scipy.spatial import cKDTree

D = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp/results/03_membrane_prod")
FLOOR = 1.80          # A, target minimum non-bonded separation
HARD = 1.20           # A, refuse to ship anything below this
ITERS = 60

print("loading topology and coordinates...")
p = parmed.load_file(str(D / "system.parm7"), xyz=str(D / "system.rst7"))
xyz = np.array(p.coordinates, dtype=float)
print(f"atoms {len(p.atoms)}")

# Exclude 1-2 and 1-3 neighbours: those are meant to be close.
excl = set()
for b in p.bonds:
    excl.add(frozenset((b.atom1.idx, b.atom2.idx)))
for a in p.angles:
    excl.add(frozenset((a.atom1.idx, a.atom3.idx)))
print(f"excluded bonded/angle pairs: {len(excl)}")

def close_pairs(coords, cutoff):
    t = cKDTree(coords)
    out = []
    for i, j in t.query_pairs(cutoff, output_type="ndarray"):
        if frozenset((int(i), int(j))) not in excl:
            out.append((int(i), int(j)))
    return out

start = close_pairs(xyz, FLOOR)
if not start:
    print("no non-bonded contacts below the floor; nothing to do")
    sys.exit(0)
d0 = [float(np.linalg.norm(xyz[i] - xyz[j])) for i, j in start]
print(f"\nnon-bonded pairs under {FLOOR} A : {len(start)}")
print(f"worst contact                 : {min(d0):.3f} A")

touched = set()
for it in range(ITERS):
    pairs = close_pairs(xyz, FLOOR)
    if not pairs:
        print(f"resolved after {it} iterations")
        break
    shift = np.zeros_like(xyz)
    for i, j in pairs:
        v = xyz[i] - xyz[j]
        d = float(np.linalg.norm(v))
        if d < 1e-6:
            v = np.random.default_rng(i * 7919 + j).normal(size=3)
            d = float(np.linalg.norm(v))
        need = (FLOOR - d) / 2.0 + 1e-3
        u = v / d
        shift[i] += u * need
        shift[j] -= u * need
        touched.add(i); touched.add(j)
    # damp so overlapping corrections do not overshoot
    xyz += shift * 0.6
else:
    print(f"still {len(close_pairs(xyz, FLOOR))} pairs after {ITERS} iterations")

final = close_pairs(xyz, FLOOR)
dmin = min((float(np.linalg.norm(xyz[i] - xyz[j])) for i, j in final),
           default=FLOOR)
print(f"\natoms moved        : {len(touched)} of {len(p.atoms)} "
      f"({100*len(touched)/len(p.atoms):.3f} %)")
print(f"pairs still under {FLOOR} A: {len(final)}")
print(f"worst remaining    : {dmin:.3f} A")

moved = np.linalg.norm(xyz - np.array(p.coordinates, dtype=float), axis=1)
print(f"largest single displacement: {moved.max():.3f} A")
print(f"mean displacement of moved atoms: {moved[list(touched)].mean():.3f} A"
      if touched else "")

if dmin < HARD:
    sys.exit(f"REFUSED: a contact at {dmin:.3f} A remains, below the {HARD} A "
             f"floor; minimisation would overflow again")

p.coordinates = xyz
p.save(str(D / "system_relaxed.rst7"), overwrite=True)
p.save(str(D / "system_relaxed.gro"), overwrite=True)
print(f"\nWrote system_relaxed.gro / .rst7")
