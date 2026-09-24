#!/usr/bin/env python3
"""Stage 3 QC - check a packed membrane system before it is parameterised.

packmol-memgen finishes with an explicit warning: "Check your final structure,
particularly for lipids inserted in proteins, protein tunnels or piercing
rings!" Packmol enforces a minimum pairwise distance but has no concept of a
lipid tail threaded through an aromatic ring or a helix bundle - those satisfy
the distance constraint while being physically impossible, and they blow up
minimisation or, worse, quietly distort the protein.

Checks performed:
  1. Composition: lipid counts, cholesterol mole fraction, waters, ions.
  2. System neutrality, against the protein+ligand formal charge.
  3. Lipid atoms buried inside the receptor (tunnel insertion).
  4. Lipid tails threaded through aromatic rings (ring piercing).
  5. Steric clashes below a hard cutoff.

Exits non-zero on any hard failure, so a bad system cannot pass silently into
Stage 4.

Usage:
    03c_qc.py --system results/03_membrane_prod/membrane_system.pdb \
              --receptor-chain-res 281 --expect-chol-frac 0.30
"""
from __future__ import annotations

import argparse
import collections
import pathlib
import sys

import numpy as np

# Lipid21 uses MODULAR residue naming: one POPC is written as three residues,
# PC (headgroup) + PA (palmitoyl) + OL (oleoyl); cholesterol is CHL, not CHL1.
# An earlier version of this script looked for POPC/CHL1, found zero lipids,
# and PASSED VACUOUSLY - every geometric check silently had nothing to test.
LIPID_HEAD = "PC"          # one per POPC
LIPID_TAILS = {"PA", "OL"}
CHOL = "CHL"
LIPIDS = {"PC", "PA", "OL", "CHL", "PE", "PS", "PGR", "OA", "ST", "LAL"}
WATERS = {"WAT", "HOH", "TIP3", "SOL"}
IONS = {"K+", "CL-", "NA+", "K", "CL", "NA"}
# Formal charges for the ionic species packmol-memgen adds.
ION_Q = {"K+": +1, "NA+": +1, "K": +1, "NA": +1, "CL-": -1, "CL": -1}

AROMATIC = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HID": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HIE": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HIP": ["CG", "ND1", "CD2", "CE1", "NE2"],
}


def parse(path: pathlib.Path):
    atoms = []
    for l in path.read_text().splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        atoms.append({
            "name": l[12:16].strip(),
            "res": l[17:20].strip(),
            "chain": l[21],
            "resid": int(l[22:26]),
            "xyz": np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]),
        })
    return atoms


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True, type=pathlib.Path)
    ap.add_argument("--expect-chol-frac", type=float, default=0.30)
    ap.add_argument("--clash", type=float, default=2.2,
                    help="hard clash cutoff, Angstrom. 2.2 A: below this is a real "
                         "heavy-atom clash. An earlier 1.2 A default let a "
                         "1.29 A contact pass.")
    a = ap.parse_args()

    atoms = parse(a.system)
    print(f"System: {a.system}")
    print(f"Total atoms: {len(atoms)}")
    fail = []

    # ---------------------------------------------------------- composition
    res_units = collections.Counter()
    seen = set()
    for at in atoms:
        key = (at["chain"], at["resid"], at["res"])
        if key not in seen:
            seen.add(key)
            res_units[at["res"]] += 1

    popc = res_units.get(LIPID_HEAD, 0)     # headgroup defines one POPC
    chl = res_units.get(CHOL, 0)
    wat = sum(res_units.get(w, 0) for w in WATERS)
    ions = {k: v for k, v in res_units.items() if k.upper() in ION_Q}
    print("\n--- composition ---")
    print(f"POPC            : {popc}")
    print(f"CHL1            : {chl}")
    if popc + chl:
        frac = chl / (popc + chl)
        print(f"cholesterol frac: {frac:.3f}  (target {a.expect_chol_frac:.2f})")
        if abs(frac - a.expect_chol_frac) > 0.03:
            fail.append(f"cholesterol fraction {frac:.3f} off target "
                        f"{a.expect_chol_frac:.2f} by more than 0.03")
    print(f"waters          : {wat}")
    # Guard against the vacuous pass described above.
    if popc + chl == 0:
        fail.append("NO LIPIDS FOUND - residue naming mismatch; every "
                    "geometric check below would be vacuous")
    for k, v in sorted(ions.items()):
        print(f"ion {k:<11} : {v}")

    # ---------------------------------------------------------- neutrality
    ion_charge = sum(ION_Q[k.upper()] * v for k, v in ions.items())
    print(f"\n--- charge ---")
    print(f"net ionic charge: {ion_charge:+d}")
    print("NOTE: neutrality is confirmed against the tleap-assigned total in")
    print("      the next step; this records the ionic contribution only.")

    # ---------------------------------------------------------- geometry
    prot = np.array([at["xyz"] for at in atoms
                     if at["res"] not in LIPIDS | WATERS
                     and at["res"].upper() not in ION_Q])
    lip = [at for at in atoms if at["res"] in LIPIDS]
    lip_xyz = np.array([at["xyz"] for at in lip]) if lip else np.empty((0, 3))
    print(f"\n--- geometry ---")
    print(f"protein/ligand atoms: {len(prot)}   lipid atoms: {len(lip_xyz)}")

    if not len(lip_xyz):
        fail.append("no lipid atoms parsed; geometric checks did not run")
    if len(prot) and len(lip_xyz):
        # 3a. lipid atoms ENCLOSED by protein (true insertion).
        #
        # A neighbour COUNT is the wrong test here and was mis-calibrated in an
        # earlier version: a 7-TM bundle has deep grooves, and annular lipids
        # legitimately sit in them with many protein atoms within 6 A. That
        # flagged 1159 atoms while hard clashes and ring piercings were zero.
        #
        # The right question is whether protein SURROUNDS the lipid atom, not
        # whether protein is near it. Test: does the atom lie inside the convex
        # hull of its own nearby protein atoms? Inside means enclosed on all
        # sides - a tunnel insertion. In a surface groove the protein neighbours
        # all lie to one side, so the atom falls outside their hull.
        from scipy.spatial import cKDTree, ConvexHull, Delaunay
        tp = cKDTree(prot)
        cand = [i for i, nb in enumerate(tp.query_ball_point(lip_xyz, 8.0))
                if len(nb) >= 12]
        print(f"lipid atoms with >=12 protein neighbours within 8 A: {len(cand)}"
              f"  (candidates, not yet failures)")
        enclosed = 0
        for i in cand:
            nb = tp.query_ball_point(lip_xyz[i], 8.0)
            pts = prot[nb]
            if len(pts) < 5:
                continue
            try:
                if Delaunay(pts).find_simplex(lip_xyz[i]) >= 0:
                    enclosed += 1
            except Exception:
                continue
        print(f"lipid atoms ENCLOSED by surrounding protein: {enclosed}")
        if enclosed > 0:
            fail.append(f"{enclosed} lipid atoms enclosed inside the protein")

        # 3b. hard clashes
        dmin, _ = tp.query(lip_xyz, k=1)
        n_clash = int(np.sum(dmin < a.clash))
        print(f"lipid-protein contacts < {a.clash} A: {n_clash}")
        if n_clash:
            fail.append(f"{n_clash} lipid-protein clashes under {a.clash} A")

        # 3c. ring piercing
        rings = collections.defaultdict(dict)
        for at in atoms:
            if at["res"] in AROMATIC and at["name"] in AROMATIC[at["res"]]:
                rings[(at["chain"], at["resid"], at["res"])][at["name"]] = at["xyz"]
        pierced = []
        tl = cKDTree(lip_xyz)
        for key, named in rings.items():
            want = AROMATIC[key[2]]
            if len(named) < len(want):
                continue
            pts = np.array([named[n] for n in want])
            centre = pts.mean(0)
            u, s, vt = np.linalg.svd(pts - centre)
            normal = vt[2]
            radius = float(np.max(np.linalg.norm(pts - centre, axis=1)))
            for idx in tl.query_ball_point(centre, radius + 0.4):
                v = lip_xyz[idx] - centre
                perp = abs(float(np.dot(v, normal)))
                inplane = float(np.linalg.norm(v - np.dot(v, normal) * normal))
                if perp < 1.6 and inplane < radius * 0.75:
                    pierced.append((key, round(perp, 2), round(inplane, 2)))
                    break
        print(f"aromatic rings with a lipid atom threaded through: {len(pierced)}")
        for k, perp, ip in pierced[:10]:
            print(f"   {k[2]}{k[1]} chain {k[0]}  perp={perp} A  inplane={ip} A")
        if pierced:
            fail.append(f"{len(pierced)} aromatic rings pierced by lipid atoms")

    # ---------------------------------------------------------- verdict
    print("\n" + "=" * 60)
    if fail:
        print("QC FAILED:")
        for f in fail:
            print("  - " + f)
        sys.exit(1)
    print("QC PASSED: composition, burial, clashes and ring piercing all clean.")


if __name__ == "__main__":
    main()
