#!/usr/bin/env python3
"""Stage 3a - receptor + peptide complex, AMBER-named, in the OPM frame.

Combines what were three hardcoded scripts. Takes every path as an argument, so
the workflow can call it per peptide.

TWO SILENT-ERROR TRAPS HANDLED HERE.

1. tleap infers histidine protonation FROM THE RESIDUE NAME, but pdb2pqr writes
   the generic HIS regardless of the state it assigned. Feeding that to tleap
   silently defaults every histidine to HIE and discards the pH assignment.
   Residues are renamed from the hydrogens actually present.

2. The build instructions require OPM/PPM orientation; packmol-memgen's default
   is MEMEMBED, a different method. The complex is superposed onto the
   OPM-deposited structure instead - an exact rigid-body transform, so the
   protonation assignment is preserved.
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

AA = {"ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYX", "CYM", "GLN", "GLU",
      "GLH", "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN",
      "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


def atoms(path: pathlib.Path) -> list[str]:
    return [l for l in path.read_text().splitlines()
            if l.startswith(("ATOM", "HETATM"))]


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
        if l[12:16].strip() != "CA":
            continue
        if chain and l[21] != chain:
            continue
        out[int(l[22:26])] = np.array(
            [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True, type=pathlib.Path)
    ap.add_argument("--ligand", required=True, type=pathlib.Path)
    ap.add_argument("--opm", required=True, type=pathlib.Path)
    ap.add_argument("--outdir", required=True, type=pathlib.Path)
    ap.add_argument("--ligand-resname", default="LIG")
    ap.add_argument("--opm-chain", default="R")
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------- 1. AMBER residue naming
    rec = atoms(a.receptor)
    byres: dict[tuple, list[str]] = {}
    for l in rec:
        byres.setdefault((l[21], int(l[22:26])), []).append(l)

    renames, fixed = [], []
    for key in sorted(byres, key=lambda k: k[1]):
        block = byres[key]
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
        if new != name:
            renames.append(f"  {name}{key[1]} -> {new}")
        for l in block:
            fixed.append(l[:17] + f"{new:>3}" + l[20:])

    print(f"Receptor: {len(byres)} residues, {len(rec)} atoms")
    print(f"Renamed for AMBER: {len(renames)}")
    for r in renames:
        print(r)

    # ---------------------------------------------- 2. ligand
    lig = atoms(a.ligand)
    lig_fixed = [l[:17] + f"{a.ligand_resname:>3}" + " L" + f"{1:>4}" + l[26:]
                 for l in lig]
    print(f"Ligand: {len(lig_fixed)} atoms as {a.ligand_resname}")

    body = fixed + ["TER"] + lig_fixed
    body = [l if l == "TER" else l[:6] + f"{i+1:>5}" + l[11:]
            for i, l in enumerate(body)]

    # ---------------------------------------------- 3. OPM frame
    opm = atoms(a.opm)
    opm_ca = ca_map(opm, chain=a.opm_chain)
    cpx_ca = ca_map([l for l in body if l != "TER" and l[21:22] != "L"])
    common = sorted(set(opm_ca) & set(cpx_ca))
    print(f"\nOPM chain {a.opm_chain}: {len(opm_ca)} CA; matched {len(common)}")
    if len(common) < 100:
        sys.exit(f"FATAL: only {len(common)} CA matched the OPM structure")

    P = np.array([cpx_ca[i] for i in common])
    Q = np.array([opm_ca[i] for i in common])
    R, t = kabsch(P, Q)
    fit = float(np.sqrt(np.mean(np.sum(((R @ P.T).T + t - Q) ** 2, axis=1))))
    print(f"superposition RMSD: {fit:.3f} A")
    if fit > 1.0:
        sys.exit(f"FATAL: RMSD {fit:.3f} A - these should be the same structure")

    out = []
    for l in body:
        if l == "TER":
            out.append(l)
            continue
        v = R @ np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) + t
        out.append(l[:30] + f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}" + l[54:])

    dest = a.outdir / "complex_opm.pdb"
    dest.write_text("\n".join(out) + "\nTER\nEND\n")
    (a.outdir / "prep_renames.log").write_text(
        "Residues renamed for AMBER (state read from placed hydrogens):\n"
        + ("\n".join(renames) if renames else "  none") + "\n")

    dum = [float(l[46:54]) for l in opm if l[17:20].strip() == "DUM"]
    if dum:
        print(f"OPM bilayer: z {min(dum):.1f} to {max(dum):.1f} A "
              f"(half-thickness {(max(dum)-min(dum))/2:.1f} A)")
    print(f"Wrote {dest}")


if __name__ == "__main__":
    main()
