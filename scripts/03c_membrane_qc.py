#!/usr/bin/env python3
"""Stage 3 QC v2 - accept a packed membrane system, or say why not.

WHAT v1 MISSED, and why this version exists.

v1 checked lipid-PROTEIN contacts only. The first production pack failed on
lipid-LIPID and water-WATER overlaps it never looked at: 551 inter-molecular
pairs under 1.2 A, worst 0.090 A. GROMACS then reported an infinite force and
minimisation died at 5.2e17 kJ/mol.

v1 also never read packmol's own verdict. packmol had ALREADY said it failed -
"packing problem with the desired distance tolerance ... contains the best
solution found", STOP 173 - and that message went unread for an entire build,
a conversion and two minimisation attempts.

So this version:
  1. Reads packmol's convergence status FIRST and fails on it.
  2. Counts contacts between DIFFERENT MOLECULES, which is what packmol's
     tolerance actually governs, instead of only lipid vs protein.
  3. Still checks composition, enclosure and aromatic ring piercing.
  4. Refuses to pass when it found nothing to check (v1 once passed vacuously
     because it looked for POPC/CHL1 while Lipid21 writes PC/PA/OL/CHL).

Usage:
    03c_membrane_qc.py --dir results/03_membrane_v2 --expect-chol-frac 0.30
"""
from __future__ import annotations

import argparse
import collections
import pathlib
import sys

import numpy as np
from scipy.spatial import Delaunay, cKDTree

LIPID_HEAD = "PC"
CHOL = "CHL"
LIPIDS = {"PC", "PA", "OL", "CHL", "PE", "PS", "PGR", "OA", "ST", "LAL"}
WATERS = {"WAT", "HOH", "TIP3", "SOL"}
ION_Q = {"K+": +1, "NA+": +1, "CL-": -1}
AROMATIC = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HID": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HIE": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HIP": ["CG", "ND1", "CD2", "CE1", "NE2"],
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, type=pathlib.Path)
    ap.add_argument("--expect-chol-frac", type=float, default=0.30)
    ap.add_argument("--hard", type=float, default=1.8,
                    help="inter-molecular contacts below this are fatal")
    ap.add_argument("--warn", type=float, default=2.2)
    a = ap.parse_args()
    fail, warn = [], []

    # ---------------------------------------------------- 1. packmol verdict
    plog = a.dir / "packmol.log"
    print("=== packmol convergence ===")
    if not plog.exists():
        fail.append("packmol.log missing - cannot confirm the pack converged")
        print("  packmol.log NOT FOUND")
    else:
        txt = plog.read_text()
        ok = "Success!" in txt or "Solution written" in txt
        bad = ("best solution found" in txt
               or "desired distance tolerance" in txt
               or "STOP 173" in txt)
        if bad and not ok:
            fail.append("packmol did NOT reach its distance tolerance; it wrote "
                        "its best attempt, not a converged pack")
            print("  FAILED to reach tolerance (wrote best solution found)")
        elif ok:
            print("  reported success")
        else:
            warn.append("packmol status could not be determined from its log")
            print("  status unclear")

    # ---------------------------------------------------- 2. parse system
    sysf = a.dir / "membrane_system.pdb"
    if not sysf.exists():
        sys.exit(f"FATAL: {sysf} not found")
    names, xyz, resid = [], [], []
    for l in sysf.read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")):
            names.append(l[17:20].strip())
            resid.append((l[21], l[22:27]))
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    xyz = np.asarray(xyz)
    print(f"\n=== composition ===\natoms {len(xyz)}")

    units = collections.Counter()
    seen = set()
    for nm, rid in zip(names, resid):
        if (rid, nm) not in seen:
            seen.add((rid, nm))
            units[nm] += 1
    popc, chl = units.get(LIPID_HEAD, 0), units.get(CHOL, 0)
    wat = sum(units.get(w, 0) for w in WATERS)
    print(f"POPC {popc}   CHL {chl}   waters {wat}")
    if popc + chl == 0:
        fail.append("NO LIPIDS FOUND - residue-naming mismatch; every geometric "
                    "check below would be vacuous")
    else:
        frac = chl / (popc + chl)
        print(f"cholesterol fraction {frac:.3f} (target {a.expect_chol_frac:.2f})")
        if abs(frac - a.expect_chol_frac) > 0.03:
            fail.append(f"cholesterol fraction {frac:.3f} off target")

    # ---------------------------------------------------- 2b. PERIODIC IMAGES
    #
    # THE CHECK THAT WAS MISSING. Every earlier contact search used raw
    # Cartesian distances and was structurally blind to periodic images.
    # packmol packs into a region but knows nothing about periodic boundaries,
    # so edge lipids extend past it and overlap the opposite face. A system can
    # therefore look clean in Cartesian space - worst contact 1.767 A, no bond
    # over 3 A - while both GROMACS and sander report infinite forces.
    #
    # The box is read from packmol.inp, whose packing regions define it.
    rid_hash = np.array([hash(r) for r in resid])
    print("\n=== periodic images ===")
    box = None
    pinp = a.dir / "packmol.inp"
    if pinp.exists():
        lo = np.array([np.inf] * 3)
        hi = np.array([-np.inf] * 3)
        for line in pinp.read_text().splitlines():
            # --pbc writes a single `pbc` directive instead of packing regions.
            if line.strip().startswith("pbc "):
                v = [float(x) for x in line.split()[1:7]]
                lo = np.minimum(lo, v[:3]); hi = np.maximum(hi, v[3:])
                continue
            if "inside box" in line:
                v = [float(x) for x in line.split()[2:8]]
                lo = np.minimum(lo, v[:3])
                hi = np.maximum(hi, v[3:])
        if np.all(np.isfinite(lo)):
            box = hi - lo
            print(f"  box from packing regions: {np.round(box, 2)}")
    if box is None:
        warn.append("box could not be determined; periodic check skipped")
    else:
        ext = xyz.max(0) - xyz.min(0)
        print(f"  coordinate extent       : {np.round(ext, 2)}")
        over = ext - box
        print(f"  overhang (extent - box) : {np.round(over, 2)}")
        periodic_pack = "pbc " in pinp.read_text()
        if np.any(over > 1.0) and not periodic_pack:
            fail.append(f"coordinates overhang the box by {np.round(over,1)} A; "
                        f"atoms will overlap their own periodic images")
        elif np.any(over > 1.0):
            print("  overhang is expected for a --pbc pack: molecules straddle")
            print("  the boundary. The minimum-image census below is the test.")
        w = xyz - xyz.min(0)
        wrapped = w - np.floor(w / box) * box
        tpbc = cKDTree(wrapped, boxsize=box)
        for cut in (0.5, 1.2, a.hard):
            pr = tpbc.query_pairs(cut, output_type="ndarray")
            if not len(pr):
                print(f"  minimum-image < {cut:.1f} A : 0")
                continue
            inter = pr[rid_hash[pr[:, 0]] != rid_hash[pr[:, 1]]]
            print(f"  minimum-image < {cut:.1f} A : {len(inter)} inter-residue")
            if cut <= 1.2 and len(inter):
                fail.append(f"{len(inter)} inter-residue contacts under {cut} A "
                            f"under MINIMUM IMAGE - the pack is not periodic")

    # ---------------------------------------------------- 3. inter-molecular
    # Molecule identity from residue id, with the three Lipid21 fragments of one
    # POPC (PC + PA + OL) treated as one molecule via their shared chain+resid
    # block is not reliable here, so use residue id and report lipid-internal
    # contacts separately rather than counting them as clashes.
    print(f"\n=== inter-molecular contacts ===")
    rid_arr = rid_hash
    t = cKDTree(xyz)
    for cut in (1.2, 1.5, a.hard, a.warn):
        pr = t.query_pairs(cut, output_type="ndarray")
        if not len(pr):
            print(f"  < {cut:.1f} A : 0")
            continue
        inter = pr[rid_arr[pr[:, 0]] != rid_arr[pr[:, 1]]]
        print(f"  < {cut:.1f} A : {len(inter)} between different residues")
    pr = t.query_pairs(a.hard, output_type="ndarray")
    inter = pr[rid_arr[pr[:, 0]] != rid_arr[pr[:, 1]]] if len(pr) else np.empty((0, 2), int)
    if len(inter):
        d = np.linalg.norm(xyz[inter[:, 0]] - xyz[inter[:, 1]], axis=1)
        print(f"  worst: {d.min():.3f} A")
        kinds = collections.Counter(
            tuple(sorted((names[i], names[j]))) for i, j in inter)
        for k, v in kinds.most_common(8):
            print(f"    {k[0]:<5} <-> {k[1]:<5} {v}")
        fail.append(f"{len(inter)} inter-residue contacts under {a.hard} A "
                    f"(worst {d.min():.3f} A) - minimisation will overflow")

    # ---------------------------------------------------- 4. enclosure / rings
    prot_mask = [n not in LIPIDS | WATERS and n.upper() not in ION_Q for n in names]
    prot = xyz[np.array(prot_mask)]
    lip_idx = np.array([i for i, n in enumerate(names) if n in LIPIDS])
    print(f"\n=== geometry ===\nprotein atoms {len(prot)}   lipid atoms {len(lip_idx)}")
    if len(prot) and len(lip_idx):
        lx = xyz[lip_idx]
        tp = cKDTree(prot)
        enclosed = 0
        for i, nb in enumerate(tp.query_ball_point(lx, 8.0)):
            if len(nb) < 12:
                continue
            try:
                if Delaunay(prot[nb]).find_simplex(lx[i]) >= 0:
                    enclosed += 1
            except Exception:
                pass
        print(f"lipid atoms enclosed by protein: {enclosed}")
        if enclosed:
            warn.append(f"{enclosed} lipid atoms enclosed by protein")

        rings = collections.defaultdict(dict)
        for i, n in enumerate(names):
            if n in AROMATIC:
                pass
        # ring piercing, using atom names from the file
        byres = collections.defaultdict(dict)
        for l in sysf.read_text().splitlines():
            if l.startswith(("ATOM", "HETATM")):
                rn = l[17:20].strip()
                if rn in AROMATIC and l[12:16].strip() in AROMATIC[rn]:
                    byres[(l[21], l[22:27], rn)][l[12:16].strip()] = np.array(
                        [float(l[30:38]), float(l[38:46]), float(l[46:54])])
        tl = cKDTree(lx)
        pierced = 0
        for key, named in byres.items():
            want = AROMATIC[key[2]]
            if len(named) < len(want):
                continue
            pts = np.array([named[n] for n in want])
            c = pts.mean(0)
            _, _, vt = np.linalg.svd(pts - c)
            nrm = vt[2]
            rad = float(np.max(np.linalg.norm(pts - c, axis=1)))
            for idx in tl.query_ball_point(c, rad + 0.4):
                v = lx[idx] - c
                if abs(float(np.dot(v, nrm))) < 1.6 and \
                   float(np.linalg.norm(v - np.dot(v, nrm) * nrm)) < rad * 0.75:
                    pierced += 1
                    break
        print(f"aromatic rings pierced by lipid: {pierced}")
        if pierced:
            fail.append(f"{pierced} aromatic rings pierced - minimisation "
                        f"cannot undo a threaded tail")

    # ---------------------------------------------------- verdict
    print("\n" + "=" * 60)
    for w in warn:
        print("WARN: " + w)
    if fail:
        print("QC FAILED:")
        for f in fail:
            print("  - " + f)
        sys.exit(1)
    print("QC PASSED")


if __name__ == "__main__":
    main()
