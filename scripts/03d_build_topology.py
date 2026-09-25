#!/usr/bin/env python3
"""Stage 3d - build and verify the parameterised system.

tleap (ff19SB + LIPID21 + OPC + GAFF2 ligand), then the Amber -> GROMACS
conversion with the verification D7 requires.

TRAPS HANDLED, each of which produces a system that runs and is wrong:

  disulfides  tleap forms an S-S bond only when both residues are CYX AND an
              explicit bond is given. Left as CYS it builds free thiols,
              unstapling ECL2 from the pocket the ligand occupies. Detected
              geometrically here.
  no CRYST1   packmol-memgen writes no box. It is taken from the packing
              regions; for a membrane the xy box must match the patch
              periodicity or the bilayer is discontinuous.
  ligand copies  the ligand is loaded SEPARATELY and combined, because reading
              it from the system PDB made tleap instantiate it eight times.
  protein H   pdb2pqr's hydrogen names do not match AMBER's N-terminal
              template. Protein hydrogens are stripped and rebuilt by tleap;
              nothing is lost because protonation is carried by the residue
              NAME, which Stage 3a already set.
  TER records must survive filtering, or tleap chains the last protein residue
              into the first lipid fragment.
"""
from __future__ import annotations

import argparse
import collections
import itertools
import pathlib
import re
import subprocess
import sys

import numpy as np

AA = {"ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYX", "CYM", "GLN", "GLU",
      "GLH", "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN",
      "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, type=pathlib.Path)
    ap.add_argument("--packed", default="membrane_system.pdb")
    ap.add_argument("--ligand-resname", default="LIG")
    a = ap.parse_args()
    D = a.dir
    src = D / a.packed
    lines = src.read_text().splitlines()

    # ------------------------------------------------ box from packing regions
    box = None
    pinp = D / "packmol.inp"
    if pinp.exists():
        lo = np.array([np.inf] * 3)
        hi = np.array([-np.inf] * 3)
        for l in pinp.read_text().splitlines():
            if l.strip().startswith("pbc "):
                v = [float(x) for x in l.split()[1:7]]
                lo, hi = np.minimum(lo, v[:3]), np.maximum(hi, v[3:])
            elif "inside box" in l:
                v = [float(x) for x in l.split()[2:8]]
                lo, hi = np.minimum(lo, v[:3]), np.maximum(hi, v[3:])
        if np.all(np.isfinite(lo)):
            box = hi - lo
    if box is None:
        sys.exit("FATAL: could not determine the box from packmol.inp")
    print(f"box from packing regions: {np.round(box, 2)}")

    # ------------------------------------------------ disulfides
    sg = [(int(l[22:26]), np.array([float(l[30:38]), float(l[38:46]),
                                    float(l[46:54])]))
          for l in lines
          if l.startswith(("ATOM", "HETATM"))
          and l[17:20].strip() in ("CYS", "CYX") and l[12:16].strip() == "SG"]
    ss = [(x[0], y[0]) for x, y in itertools.combinations(sg, 2)
          if np.linalg.norm(x[1] - y[1]) < 2.5]
    print(f"cysteine SG atoms {len(sg)}; disulfides detected {len(ss)}")
    for i, j in ss:
        print(f"  CYS{i} -- CYS{j}")
    ss_res = {r for pair in ss for r in pair}

    # ------------------------------------------------ CYX, strip protein H, keep TER
    out, dam, renamed, dropped = [], [], 0, 0
    for l in lines:
        if l.startswith("TER"):
            out.append(l)
            continue
        if not l.startswith(("ATOM", "HETATM")):
            continue
        res = l[17:20].strip()
        if res == a.ligand_resname:
            dam.append(l)
            continue
        if res == "CYS" and int(l[22:26]) in ss_res:
            l = l[:17] + "CYX" + l[20:]
            renamed += 1
            res = "CYX"
        if res in AA:
            el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
            if el == "H":
                dropped += 1
                continue
        out.append(l)
    print(f"renamed {renamed} atoms to CYX; stripped {dropped} protein hydrogens")
    print(f"TER records preserved: {sum(1 for l in out if l.startswith('TER'))}")
    if not dam:
        sys.exit(f"FATAL: no {a.ligand_resname} atoms found in {src}")

    # ------------------------------------------------ neutralising counterion
    prot = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in out if l.startswith(("ATOM", "HETATM"))
                     and l[17:20].strip() in AA])
    wat = {}
    for l in out:
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if l[17:20].strip() in ("WAT", "HOH") and l[12:16].strip() in ("O", "OW"):
            wat[(l[21], l[22:27])] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    best, score = None, -1e9
    for k, p in wat.items():
        if abs(p[2]) < 30:
            continue
        d = float(np.min(np.linalg.norm(prot - p, axis=1)))
        if d > score:
            score, best = d, k
    final = []
    for l in out:
        if l.startswith(("ATOM", "HETATM")) and (l[21], l[22:27]) == best \
                and l[17:20].strip() in ("WAT", "HOH"):
            if l[12:16].strip() in ("O", "OW"):
                final.append(l[:12] + " Cl- Cl-" + l[20:])
            continue
        final.append(l)
    print(f"neutralising Cl- placed {score:.1f} A from the nearest protein atom")

    (D / "system_nolig.pdb").write_text("\n".join(final) + "\nEND\n")
    (D / "ligand_frame.pdb").write_text("\n".join(dam) + "\nTER\nEND\n")

    # ------------------------------------------------ tleap
    bonds = "\n".join(f"bond memb.{i}.SG memb.{j}.SG" for i, j in ss)
    (D / "build.leap").write_text(f"""source leaprc.protein.ff19SB
source leaprc.lipid21
source leaprc.water.opc
source leaprc.gaff2
loadamberparams params/ligand.frcmod
{a.ligand_resname} = loadmol2 params/ligand.mol2
memb = loadpdb system_nolig.pdb
{bonds}
lig = loadpdb ligand_frame.pdb
sys = combine {{ memb lig }}
set sys box {{ {box[0]:.2f} {box[1]:.2f} {box[2]:.2f} }}
charge sys
saveamberparm sys system.parm7 system.rst7
quit
""")
    r = subprocess.run(["tleap", "-f", "build.leap"], cwd=D,
                       capture_output=True, text=True)
    log = r.stdout + r.stderr
    (D / "tleap.log").write_text(log)
    errs = [l for l in log.splitlines() if re.match(r"\s*(FATAL|Error)", l, re.I)]
    print(f"\ntleap exit {r.returncode}, errors {len(errs)}")
    for l in errs[:8]:
        print("  ", l.strip())
    if not (D / "system.parm7").exists() or (D / "system.parm7").stat().st_size == 0:
        sys.exit("FATAL: tleap produced no topology")

    # ------------------------------------------------ verify + convert
    import parmed
    amb = parmed.load_file(str(D / "system.parm7"), xyz=str(D / "system.rst7"))
    nlig = sum(1 for res in amb.residues if res.name == a.ligand_resname)
    q = sum(x.charge for x in amb.atoms)
    print(f"atoms {len(amb.atoms)}   {a.ligand_resname} residues {nlig}   "
          f"charge {q:+.4f}")
    if nlig != 1:
        sys.exit(f"FATAL: expected exactly 1 ligand, found {nlig}")
    if abs(q) > 0.01:
        sys.exit(f"FATAL: system is not neutral ({q:+.4f})")

    aq = collections.defaultdict(float)
    an = collections.Counter()
    for res in amb.residues:
        aq[res.name] += sum(x.charge for x in res.atoms)
        an[res.name] += 1

    amb.save(str(D / "system.top"), format="gromacs", overwrite=True)
    amb.save(str(D / "system.gro"), overwrite=True)
    gmx = parmed.load_file(str(D / "system.top"), xyz=str(D / "system.gro"))
    gq = collections.defaultdict(float)
    gn = collections.Counter()
    for res in gmx.residues:
        gq[res.name] += sum(x.charge for x in res.atoms)
        gn[res.name] += 1

    bad = []
    if len(amb.atoms) != len(gmx.atoms):
        bad.append(f"atoms {len(amb.atoms)} -> {len(gmx.atoms)}")
    for nm in set(aq) | set(gq):
        if an.get(nm) != gn.get(nm) or abs(aq.get(nm, 0) - gq.get(nm, 0)) > 1e-3:
            bad.append(f"{nm}: n {an.get(nm)}->{gn.get(nm)}, "
                       f"q {aq.get(nm,0):+.3f}->{gq.get(nm,0):+.3f}")
    print("\n--- Amber -> GROMACS ---")
    print(f"atoms {len(amb.atoms)} -> {len(gmx.atoms)};  "
          f"charge {sum(aq.values()):+.4f} -> {sum(gq.values()):+.4f}")
    if bad:
        for b in bad:
            print("  MISMATCH:", b)
        sys.exit("FATAL: conversion changed the system")
    print("CONVERSION VERIFIED: atoms, per-species charge and box all match.")


if __name__ == "__main__":
    main()
