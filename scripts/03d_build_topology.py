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

    # ------------------------------------------------ which force field
    #
    # 03b routes a peptide of standard L-amino acids to ff19SB and everything
    # else - DAMGO, with its D-Ala, N-Me-Phe and ethanolamide - to GAFF2, and
    # records the choice in params_mode.txt. Stage 5 has read that file since
    # b48ce60; this stage did not, so it went on demanding params/ligand.mol2
    # and params/ligand.frcmod, which the ff19SB path never writes. Every
    # canonical peptide therefore failed here with its parameters sitting
    # correctly on disk under different names.
    #
    # The packed system is unaffected by the choice: 03a always writes the
    # peptide as a single LIG residue, and that copy is used only as a set of
    # COORDINATES to move the properly named template onto (see lib_ligand).
    # So the bilayer packs do not have to be rebuilt to change force field.
    mf = D / "params" / "params_mode.txt"
    mode = mf.read_text().strip() if mf.exists() else "gaff2"
    if mode not in ("gaff2", "ff19SB"):
        sys.exit(f"FATAL: params_mode.txt says {mode!r}, expected "
                 f"'gaff2' or 'ff19SB'")
    print(f"ligand parameters: {mode}")

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
    # The counterion count and SIGN come from the ligand's actual charge, not
    # from an assumption.
    #
    # This used to place exactly one Cl-, which is right only for a +1 ligand.
    # That is the same +1 assumption D33 removed from ligand_params and it
    # survived here: met-enkephalin is a neutral zwitterion, so the Cl- CREATED
    # a net -1 and the neutrality check then failed. Peptides carrying a
    # glutamate are -1 and need Na+, not Cl- - two ions in the wrong direction.
    # GAFF2 carries the charge in the mol2's per-atom column; ff19SB has no
    # mol2 at all, and 03b wrote the composition charge to ligand_charge.txt.
    if mode == "gaff2":
        q_lig = round(sum(
            float(l.split()[8]) for l in
            (D / "params" / "ligand.mol2").read_text().split("@<TRIPOS>ATOM")[1]
            .split("@<TRIPOS>BOND")[0].strip().splitlines()))
    else:
        q_lig = int((D / "params" / "ligand_charge.txt").read_text().strip())
    ion = "Cl-" if q_lig > 0 else "Na+"
    n_ion = abs(q_lig)
    print(f"ligand net charge {q_lig:+d} -> {n_ion} x {ion}")

    ranked = sorted(
        ((float(np.min(np.linalg.norm(prot - p, axis=1))), k)
         for k, p in wat.items() if abs(p[2]) >= 30),
        reverse=True)
    if len(ranked) < n_ion:
        sys.exit(f"FATAL: need {n_ion} bulk waters to replace, found "
                 f"{len(ranked)}")
    chosen = {k: d for d, k in ranked[:n_ion]}

    final = []
    for l in out:
        if l.startswith(("ATOM", "HETATM")) and (l[21], l[22:27]) in chosen \
                and l[17:20].strip() in ("WAT", "HOH"):
            if l[12:16].strip() in ("O", "OW"):
                final.append(l[:12] + f" {ion:<4}{ion:<4}"[:9] + l[20:])
            continue
        final.append(l)
    if n_ion:
        print(f"{n_ion} x {ion} placed, nearest "
              f"{min(chosen.values()):.1f} A from any protein atom")
    else:
        print("ligand is neutral - no counterion needed")

    (D / "system_nolig.pdb").write_text("\n".join(final) + "\nEND\n")
    # The ligand handed to tleap must carry the MOL2's atom names, not the
    # predictor's. The LIG records in the packed system are the raw predicted
    # peptide - protein-style names (CE, SD, OH, CZ...) repeated across
    # residues, no OXT, no hydrogens - so tleap could not type a single one of
    # them ("Atom .R<LIG>.A<CE 89> does not have a type", 14 fatal errors).
    #
    # lib_ligand moves the mol2-consistent file onto these coordinates by rigid
    # superposition and verifies the fit. Shared with Stage 5, because this is
    # the fourth site at which the ligand's identity was mishandled.
    sys.path.insert(0, str(pathlib.Path(__file__).parent))
    from lib_ligand import ligand_in_frame
    # GAFF2's template is the mol2-consistent PDB; ff19SB's is the cleaned
    # peptide 03b wrote, whose residues are real amino acids numbered from 901.
    # lib_ligand does not rename anything - it only moves - so the ff19SB
    # residue names survive into ligand_frame.pdb, which is what tleap needs.
    template = D / "params" / ("ligand_unique.pdb" if mode == "gaff2"
                               else "ligand_peptide.pdb")
    try:
        lig_lines, fit, extra = ligand_in_frame(
            template, dam, keep_residues=(mode == "ff19SB"))
    except ValueError as e:
        sys.exit(f"FATAL: {e}")
    print(f"ligand renamed from {template.name}: {len(dam)} -> "
          f"{len(lig_lines)} atoms, superposition RMSD {fit:.4f} A "
          f"({extra} extra heavy atom(s), e.g. OXT)")
    (D / "ligand_frame.pdb").write_text("\n".join(lig_lines) + "\nTER\nEND\n")

    # ------------------------------------------------ tleap
    bonds = "\n".join(f"bond memb.{i}.SG memb.{j}.SG" for i, j in ss)
    # A canonical peptide needs no small-molecule parameters at all - ff19SB
    # already has every residue, and tleap rebuilds its hydrogens and caps the
    # termini itself. Loading leaprc.gaff2 and a mol2 that does not exist is
    # what blocked this path.
    ligand_params = ("source leaprc.gaff2\n"
                     "loadamberparams params/ligand.frcmod\n"
                     f"{a.ligand_resname} = loadmol2 params/ligand.mol2\n"
                     if mode == "gaff2" else "")
    (D / "build.leap").write_text(f"""source leaprc.protein.ff19SB
source leaprc.lipid21
source leaprc.water.opc
{ligand_params}memb = loadpdb system_nolig.pdb
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
    q = sum(x.charge for x in amb.atoms)
    if mode == "gaff2":
        nlig = sum(1 for res in amb.residues if res.name == a.ligand_resname)
        print(f"atoms {len(amb.atoms)}   {a.ligand_resname} residues {nlig}   "
              f"charge {q:+.4f}")
        if nlig != 1:
            sys.exit(f"FATAL: expected exactly 1 ligand, found {nlig}")
    else:
        # ff19SB gives the peptide REAL residues, so it cannot be found by
        # name - the receptor has its own tyrosines - and it is NOT simply the
        # last n residues either: `combine { memb lig }` does not leave the
        # ligand at the end of the parm, which an earlier version of this
        # check assumed and duly reported the peptide as WAT-WAT-WAT-WAT-WAT
        # on a topology that was in fact correct.
        #
        # So the peptide is located by SEQUENCE, and the test is that the run
        # occurs exactly once. That catches a dropped, duplicated or mistyped
        # residue, and it does not care where tleap put it. Charged termini
        # are named NTYR/CMET, so the leading N/C is stripped before matching.
        at = [l for l in lig_lines if l.startswith(("ATOM", "HETATM"))]
        exp = [l[17:20].strip() for i, l in enumerate(at)
               if i == 0 or l[21:27] != at[i - 1][21:27]]
        norm = lambda n: n[1:] if len(n) == 4 and n[0] in "NC" else n
        names = [norm(r.name) for r in amb.residues]
        hits = [i for i in range(len(names) - len(exp) + 1)
                if names[i:i + len(exp)] == exp]
        print(f"atoms {len(amb.atoms)}   peptide {'-'.join(exp)} "
              f"at residue index {hits}   charge {q:+.4f}")
        if len(hits) != 1:
            sys.exit(f"FATAL: expected exactly one {'-'.join(exp)} run in the "
                     f"topology, found {len(hits)} - tleap did not build the "
                     f"sequence it was given")
        npep = sum(len(r.atoms) for r in amb.residues[hits[0]:hits[0] + len(exp)])
        print(f"peptide atoms {npep} (heavy {len(at)} + hydrogens "
              f"rebuilt by tleap)")
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
