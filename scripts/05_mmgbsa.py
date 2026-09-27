#!/usr/bin/env python3
"""Stage 5 - MM/GBSA triage scoring.

PURPOSE IS RANK ORDERING, NOT AFFINITY. The number this writes is an
interaction energy, not a binding free energy: there is no configurational
entropy, no explicit solvent, and no membrane. It is reported in kcal/mol
because that is its unit, not because it is comparable to an experimental dG.

WHY THIS RUNS WITHOUT MD
------------------------
The obvious implementation scores an ensemble from the production trajectory.
That is the better calculation and the wrong tool: a filter that first requires
50 ns of MD per peptide is not a filter, it is post-hoc analysis of the
expensive step it was supposed to gate.

This scores a MINIMISED COMPLEX instead, which costs minutes rather than a day,
so every peptide in the library can be ranked - including the scrambled decoys,
which is the only way to know whether the ranking means anything. The
trajectory-ensemble variant belongs downstream, on survivors.

WHAT IS APPROXIMATED, EXPLICITLY
--------------------------------
1. NO MEMBRANE. Implicit solvent treats everything as aqueous. The orthosteric
   pocket is largely occluded from lipid and every peptide is scored at the
   SAME site, so the membrane term is roughly constant across the comparison
   and cancels in a RANKING. It would not cancel in an absolute number.
2. NO ENTROPY. Normal-mode or quasi-harmonic entropy is omitted, as is
   conventional for MM/GBSA ranking. Peptides differ in flexibility, so this is
   a real and unequal error, not a constant offset.
3. SINGLE STRUCTURE. One minimised conformation, not an ensemble. No
   conformational averaging.
4. The pose comes from co-folding, which Stage 7 showed cannot distinguish
   binders from decoys. A physics score on a wrong pose is a precise number
   for the wrong structure. This is the main reason the decoy control matters.
"""
from __future__ import annotations

import argparse
import pathlib
import re
import shutil
import subprocess
import sys

import numpy as np


def atoms(path):
    return [l for l in path.read_text().splitlines()
            if l.startswith(("ATOM", "HETATM"))]


def xyz(lines):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in lines])


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def run(cmd: list[str] | str, cwd: pathlib.Path, log: str) -> subprocess.CompletedProcess:
    shell = isinstance(cmd, str)
    p = subprocess.run(cmd, cwd=cwd, shell=shell, capture_output=True, text=True)
    (cwd / log).write_text((p.stdout or "") + (p.stderr or ""))
    return p


LEAP = """source leaprc.protein.ff19SB
source leaprc.gaff2
set default PBRadii mbondi3
loadamberparams {frcmod}
LIG = loadmol2 {mol2}
rec = loadpdb receptor.pdb
{bonds}
lig = loadpdb ligand_frame.pdb
cpx = combine {{ rec lig }}
saveamberparm cpx complex.parm7 complex.rst7
quit
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--complex", required=True, type=pathlib.Path,
                    help="complex_opm.pdb - receptor chain R + LIG chain L")
    ap.add_argument("--mol2", required=True, type=pathlib.Path)
    ap.add_argument("--frcmod", required=True, type=pathlib.Path)
    ap.add_argument("--ligand-pdb", required=True, type=pathlib.Path,
                    help="params/ligand_unique.pdb - the ligand with the SAME "
                         "atom names and count as the mol2.")
    ap.add_argument("--outdir", required=True, type=pathlib.Path)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    ap.add_argument("--peptide", required=True)
    ap.add_argument("--igb", type=int, default=8,
                    help="GB model. 8 (GBn2) with mbondi3 radii is the "
                         "AmberTools-recommended pairing.")
    ap.add_argument("--salt", type=float, default=0.15)
    ap.add_argument("--minsteps", type=int, default=2000)
    ap.add_argument("--disulfide", default="140,217",
                    help="TM3-ECL2 CYX pair, receptor numbering (D-notes).")
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)
    for f in (a.mol2, a.frcmod, a.complex):
        shutil.copy(f, a.outdir / f.name)

    # -------------------------------------------- 1. put the ligand in frame
    #
    # complex_opm.pdb's LIG cannot be used directly. It holds the 39 raw heavy
    # atoms the co-folding model emitted, flattened into ONE residue - so its
    # atom names repeat (N, CA, C, O once per original residue) and it lacks
    # the C-terminal OXT and every hydrogen. The mol2 that carries the GAFF
    # parameters has 75 uniquely named atoms. tleap matched them by name,
    # decided 78 atoms were "missing", built them from its template, and
    # emitted a Fatal Error while still writing a file.
    #
    # The ligand geometry is recovered instead from params/ligand_unique.pdb,
    # which HAS the mol2's naming but sits in the prediction's own frame. Its
    # heavy atoms are in the same order as the raw ligand, so the rigid
    # transform between the two is exact and recoverable - and asserted below,
    # because a silent mis-superposition would place a correctly-named ligand
    # in the wrong place.
    cpx = atoms(a.complex)
    lig_raw = [l for l in cpx if l[17:20].strip() == "LIG"]
    rec = [l for l in cpx if l[17:20].strip() != "LIG"]
    if not lig_raw:
        sys.exit(f"FATAL: no LIG residue in {a.complex}")

    lu = atoms(a.ligand_pdb)
    lu_heavy = [l for l in lu
                if (l[76:78].strip() or l[12:16].strip()[0]).upper() != "H"]
    n = len(lig_raw)
    if len(lu_heavy) < n:
        sys.exit(f"FATAL: {a.ligand_pdb} has {len(lu_heavy)} heavy atoms, "
                 f"fewer than the {n} in the complex - not the same molecule.")

    R, tvec = kabsch(xyz(lu_heavy[:n]), xyz(lig_raw))
    fit = float(np.sqrt(np.mean(np.sum(
        ((R @ xyz(lu_heavy[:n]).T).T + tvec - xyz(lig_raw)) ** 2, axis=1))))
    print(f"ligand frame transfer: {n} heavy atoms, RMSD {fit:.4f} A "
          f"({len(lu_heavy) - n} extra heavy atom(s) carried along, e.g. OXT)")
    if fit > 0.1:
        sys.exit(f"FATAL: {fit:.3f} A - these should be the SAME conformation "
                 f"in two frames. Atom order does not correspond.")

    moved = []
    for l in lu:
        v = R @ np.array([float(l[30:38]), float(l[38:46]),
                          float(l[46:54])]) + tvec
        moved.append(l[:17] + "LIG L   1" + l[26:30]
                     + f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}" + l[54:])
    (a.outdir / "ligand_frame.pdb").write_text("\n".join(moved) + "\nTER\nEND\n")

    # The receptor needs two edits before tleap will accept it, both of which
    # Stage 3 also performs - and both of which fail LOUDLY here rather than
    # silently, which is the only reason they were found.
    #
    # 1. STRIP PROTEIN HYDROGENS. These were placed by PROPKA/pdb2pqr, whose
    #    naming does not match ff19SB's templates - tleap stopped with
    #    "Atom .R<NMET 65>.A<H 20> does not have a type". tleap rebuilds them
    #    correctly from its own templates, so the placed ones are discarded.
    #    Protonation STATE is preserved regardless, because that is carried by
    #    the residue NAME (HID/HIE/ASH), not by which hydrogens are present.
    #
    # 2. RENAME THE DISULFIDE CYSTEINES TO CYX. The file says CYS, and a CYS
    #    loaded as CYS keeps its HG and cannot form the bond. Naming them CYX
    #    and declaring the bond is what actually creates the TM3-ECL2
    #    disulfide, which the 50 ns DAMGO run showed contacting the ligand in
    #    94% of frames.
    ss = [x.strip() for x in a.disulfide.split(",")] if a.disulfide else []
    kept, dropped, renamed = [], 0, 0
    for l in rec:
        el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        if el == "H":
            dropped += 1
            continue
        if l[22:26].strip() in ss and l[17:20].strip() == "CYS":
            l = l[:17] + "CYX" + l[20:]
            renamed += 1
        kept.append(l)
    print(f"receptor: stripped {dropped} hydrogens, renamed {renamed} atoms "
          f"to CYX at {ss}")
    if ss and renamed == 0:
        sys.exit(f"FATAL: no CYS found at {ss} - the disulfide cannot be "
                 f"formed and the residue numbering is not what was assumed.")
    (a.outdir / "receptor.pdb").write_text("\n".join(kept) + "\nTER\nEND\n")

    # ---------------------------------------------------------- 2. topology
    # The disulfide must be declared or tleap silently leaves two free
    # cysteines - the same failure class as the HID naming in Stage 3.
    bonds = ""
    if a.disulfide:
        i, j = a.disulfide.split(",")
        bonds = f"bond rec.{i}.SG rec.{j}.SG"
    (a.outdir / "leap.in").write_text(LEAP.format(
        frcmod=a.frcmod.name, mol2=a.mol2.name, bonds=bonds))
    run(["tleap", "-f", "leap.in"], a.outdir, "leap.log")

    # tleap writes a topology file even when it reports a Fatal Error, so
    # existence is not success. This check was added after exactly that: the
    # file was there, the run had failed, and the next tool reported an
    # unrelated "could not determine file type".
    leaplog = (a.outdir / "leap.log").read_text()
    m = re.search(r"Errors\s*=\s*(\d+)", leaplog)
    nerr = int(m.group(1)) if m else -1
    parm = a.outdir / "complex.parm7"
    if nerr != 0 or not parm.exists() or parm.stat().st_size == 0:
        for line in [l for l in leaplog.splitlines()
                     if re.search(r"error|fatal|missing", l, re.I)][:10]:
            print("   ", line.strip())
        sys.exit(f"FATAL: tleap reported {nerr} error(s). "
                 f"See {a.outdir}/leap.log")
    print(f"tleap: 0 errors, topology written")

    # ---------------------------------------------------- 2. split topologies
    # Single-trajectory MM/GBSA: receptor and ligand topologies are DERIVED
    # from the complex, so all three share identical parameters and the
    # interaction energy is a difference of like with like.
    p = run(["ante-MMPBSA.py", "-p", "complex.parm7",
             "-c", "com.parm7", "-r", "rec.parm7", "-l", "lig.parm7",
             "-s", ":WAT,Na+,Cl-", "-n", ":LIG", "--radii=mbondi3"],
            a.outdir, "ante.log")
    for f in ("com.parm7", "rec.parm7", "lig.parm7"):
        if not (a.outdir / f).exists():
            sys.exit(f"FATAL: ante-MMPBSA.py did not write {f}. "
                     f"See {a.outdir}/ante.log")

    # ------------------------------------------------------- 3. minimisation
    # Co-folded poses arrive with real clashes - closest contacts around 1.0-1.3
    # A - and an unminimised MM/GBSA is dominated by those, producing large
    # positive van der Waals terms that say more about the predictor than the
    # peptide. Minimisation is in the SAME implicit solvent the scoring uses,
    # so the structure is relaxed on the surface it is scored on.
    (a.outdir / "min.in").write_text(
        "minimise complex in implicit solvent\n"
        " &cntrl\n"
        f"  imin=1, maxcyc={a.minsteps}, ncyc={a.minsteps // 2},\n"
        f"  igb={a.igb}, saltcon={a.salt}, gbsa=1,\n"
        "  cut=999.0, ntb=0, ntpr=200,\n"
        " /\n")
    p = run(["sander", "-O", "-i", "min.in", "-p", "complex.parm7",
             "-c", "complex.rst7", "-r", "min.rst7", "-o", "min.out"],
            a.outdir, "sander.log")
    if not (a.outdir / "min.rst7").exists():
        sys.exit(f"FATAL: minimisation failed. See {a.outdir}/min.out")

    # cpptraj converts the restart into a one-frame trajectory, which is the
    # input form MMPBSA.py expects.
    (a.outdir / "traj.in").write_text(
        "parm complex.parm7\ntrajin min.rst7\ntrajout min.mdcrd\nrun\nquit\n")
    run(["cpptraj", "-i", "traj.in"], a.outdir, "cpptraj.log")

    # --------------------------------------------------------- 4. MM/GBSA
    (a.outdir / "mmgbsa.in").write_text(
        "MM/GBSA, single trajectory, one minimised frame\n"
        "&general\n  startframe=1, endframe=1, interval=1, verbose=2,\n/\n"
        f"&gb\n  igb={a.igb}, saltcon={a.salt},\n/\n")
    p = run(["MMPBSA.py", "-O", "-i", "mmgbsa.in", "-o", "mmgbsa.dat",
             "-sp", "complex.parm7", "-cp", "com.parm7",
             "-rp", "rec.parm7", "-lp", "lig.parm7", "-y", "min.mdcrd"],
            a.outdir, "mmpbsa.log")

    dat = a.outdir / "mmgbsa.dat"
    if not dat.exists():
        sys.exit(f"FATAL: MMPBSA.py wrote no output. See {a.outdir}/mmpbsa.log")

    # ----------------------------------------------------------- 5. parse
    text = dat.read_text()
    m = re.search(r"DELTA TOTAL\s+(-?\d+\.\d+)", text)
    if not m:
        sys.exit(f"FATAL: no 'DELTA TOTAL' in {dat} - the calculation ran but "
                 f"produced no interaction energy.")
    total = float(m.group(1))

    # Parse components ONLY from the "Differences" block.
    #
    # mmgbsa.dat lists Complex, Receptor, Ligand and then Differences, each
    # with the same component names. A plain search finds the COMPLEX's
    # absolute energies first - which is what the first version of this
    # reported, giving a breakdown summing to about -24,000 kcal/mol beside a
    # total of -40.7 and looking superficially like a component table.
    #
    # The sum check below is what makes that non-recurring: the components of
    # a difference must add up to the difference.
    diff = text.split("Differences (Complex - Receptor - Ligand)")
    if len(diff) < 2:
        sys.exit(f"FATAL: no Differences block in {dat}")
    block = diff[1]

    terms = {}
    for key in ("VDWAALS", "EEL", "EGB", "ESURF"):
        mm = re.search(rf"^{key}\s+(-?\d+\.\d+)", block, re.M)
        if mm:
            terms[key] = float(mm.group(1))

    ssum = sum(terms.values())
    if abs(ssum - total) > 0.5:
        sys.exit(f"FATAL: components sum to {ssum:.2f} but DELTA TOTAL is "
                 f"{total:.2f}. The wrong section was parsed.")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w") as fh:
        fh.write("peptide\tmetric\tvalue_kcal_per_mol\n")
        fh.write(f"{a.peptide}\tdelta_total\t{total:.3f}\n")
        for k, v in terms.items():
            fh.write(f"{a.peptide}\t{k.lower()}\t{v:.3f}\n")

    print(f"{a.peptide}: DELTA TOTAL = {total:.2f} kcal/mol")
    for k, v in terms.items():
        print(f"    {k:<9} {v:>10.2f}")
    print("\nRANKING SCORE ONLY - not a binding free energy (no entropy, no "
          "membrane, single structure).")


if __name__ == "__main__":
    main()
