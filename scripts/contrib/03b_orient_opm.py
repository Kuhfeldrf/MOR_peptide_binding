#!/usr/bin/env python3
"""Stage 3b - put the complex into the OPM membrane frame.

The build instructions require OPM/PPM orientation. packmol-memgen's default is
MEMEMBED, a different method, so using it would be a silent substitution.
Instead the OPM-deposited 6DDF is downloaded (it carries DUM pseudo-atoms
marking the bilayer boundaries and a header stating half-thickness 15.7 A), and
our protonated complex is superposed onto it.

Because the OPM entry is the SAME structure as data/raw/6ddf.cif, merely
rotated and translated, this is an exact rigid-body transform: no coordinates
are distorted and Stage 0's protonation assignment is preserved intact. The
membrane normal then lies along z, as packmol-memgen --preoriented expects.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp")
OUT = ROOT / "results/03_membrane"


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def ca_map(lines, chain=None):
    out = {}
    for l in lines:
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if l[12:16].strip() != "CA":
            continue
        if chain and l[21] != chain:
            continue
        out[int(l[22:26])] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return out


opm = (ROOT / "data/raw/6ddf_opm.pdb").read_text().splitlines()
cpx = (OUT / "complex_damgo.pdb").read_text().splitlines()

opm_ca = ca_map(opm, chain="R")
cpx_ca = ca_map([l for l in cpx if l[21:22] != "L"])
common = sorted(set(opm_ca) & set(cpx_ca))
print(f"OPM chain R CA: {len(opm_ca)}   complex receptor CA: {len(cpx_ca)}")
print(f"matched by residue number: {len(common)}")
if len(common) < 200:
    sys.exit(f"FATAL: only {len(common)} CA matched; numbering may differ")

P = np.array([cpx_ca[i] for i in common])
Q = np.array([opm_ca[i] for i in common])
R, t = kabsch(P, Q)
fit = float(np.sqrt(np.mean(np.sum(((R @ P.T).T + t - Q) ** 2, axis=1))))
print(f"superposition RMSD: {fit:.3f} A")
if fit > 1.0:
    sys.exit(f"FATAL: RMSD {fit:.3f} A - these should be the same structure")

# membrane extent from the OPM DUM pseudo-atoms
dum_z = [float(l[46:54]) for l in opm
         if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "DUM"]
print(f"OPM bilayer markers: {len(dum_z)} DUM atoms, "
      f"z from {min(dum_z):.1f} to {max(dum_z):.1f} A "
      f"(half-thickness {(max(dum_z)-min(dum_z))/2:.1f} A)")

out_lines = []
for l in cpx:
    if l.startswith(("ATOM", "HETATM")):
        v = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        v = R @ v + t
        l = l[:30] + f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}" + l[54:]
    out_lines.append(l)

dest = OUT / "complex_damgo_opm.pdb"
dest.write_text("\n".join(out_lines) + "\n")

z = [float(l[46:54]) for l in out_lines
     if l.startswith(("ATOM", "HETATM"))]
print(f"\nComplex now spans z {min(z):.1f} to {max(z):.1f} A "
      f"(membrane normal along z)")
print(f"Wrote {dest}")
