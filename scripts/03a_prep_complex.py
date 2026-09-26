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
    ap.add_argument(
        "--pose", type=pathlib.Path,
        help="Co-folded receptor+peptide, used ONLY to bring the ligand into "
             "the receptor's coordinate frame. Required when the ligand comes "
             "from co-folding; omitted when it was extracted from the "
             "receptor structure and is already in frame.")
    ap.add_argument("--pose-chain", default="A",
                    help="Receptor chain within --pose (the peptide is a "
                         "separate chain and must be excluded).")
    ap.add_argument("--pose-offset", type=int, default=0,
                    help="Added to --pose residue numbers to reach the "
                         "receptor's numbering. Co-folding numbers the "
                         "construct from 1; the crystal starts at 65.")
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

    # Every atom handed in here is relabelled as a single residue LIG, so this
    # file MUST be the peptide alone. When it was accidentally wired to a
    # co-folding output containing receptor+peptide, the receptor was silently
    # duplicated into the "ligand", the bounding box tripled in every
    # dimension, and packmol packed a 1.9M-atom system. Every step still
    # exited 0.
    #
    # The test is a backbone count measured against the receptor, so there is
    # no absolute cutoff to go stale. A digestion peptide is a handful of
    # residues; anything approaching the receptor's own chain length is the
    # receptor. A positional test would NOT work here: the duplicate arrived at
    # the co-folding model's coordinates, hundreds of angstroms away, so it
    # overlapped nothing.
    rec_ca = ca_map(rec)
    lig_ca = ca_map(lig)
    if len(lig_ca) > 0.25 * len(rec_ca):
        sys.exit(
            f"FATAL: --ligand {a.ligand} has {len(lig_ca)} CA atoms against "
            f"the receptor's {len(rec_ca)} - this is a receptor+peptide "
            f"complex, not a ligand.\n"
            f"       Pass the peptide alone (best_ligand.pdb), not the "
            f"co-folded pose (best_pose.pdb).")

    # ------------------------------------- 2b. bring the ligand into frame
    #
    # A co-folding model predicts receptor and peptide together in ITS OWN
    # arbitrary frame, near the origin. The experimental receptor is at its
    # crystal coordinates - here about 230 A away. Concatenating the two
    # directly leaves the peptide hundreds of angstroms outside the protein,
    # which is geometrically absurd but breaks nothing: the box simply grows to
    # contain both, and packmol fills it.
    #
    # So the peptide is moved by superposing the PREDICTED receptor onto the
    # experimental one and applying that transform to the peptide. The
    # prediction's own receptor is the only thing that relates the two frames.
    #
    # DAMGO never needed this - it was extracted from the receptor structure and
    # is already in frame - which is exactly why the gap went unnoticed until a
    # co-folded peptide first reached this stage.
    if a.pose:
        # Two corrections are needed before the numbering lines up, and getting
        # either wrong superposes mismatched residues rather than failing:
        #   - the peptide is a separate chain in the pose and must be excluded,
        #     or its residues 1-5 collide with the receptor's own 1-5;
        #   - the prediction numbers the construct from 1, the crystal from 65.
        pose_ca = {k + a.pose_offset: v for k, v in
                   ca_map(atoms(a.pose), chain=a.pose_chain).items()}
        shared = sorted(set(pose_ca) & set(rec_ca))
        print(f"\nPose frame: chain {a.pose_chain}, {len(pose_ca)} CA, "
              f"offset {a.pose_offset:+d}, {len(shared)} shared with receptor")
        if len(shared) < 100:
            sys.exit(f"FATAL: only {len(shared)} CA shared with --pose "
                     f"{a.pose} - cannot determine the ligand's frame")
        Rp, tp = kabsch(np.array([pose_ca[i] for i in shared]),
                        np.array([rec_ca[i] for i in shared]))
        fit_p = float(np.sqrt(np.mean(np.sum(
            ((Rp @ np.array([pose_ca[i] for i in shared]).T).T + tp
             - np.array([rec_ca[i] for i in shared])) ** 2, axis=1))))
        # A predicted receptor is not identical to the crystal one, so this is
        # a model-vs-experiment RMSD, not the ~0 of a same-structure fit.
        print(f"predicted -> experimental receptor RMSD: {fit_p:.3f} A")
        if fit_p > 5.0:
            sys.exit(f"FATAL: {fit_p:.3f} A - the predicted receptor does not "
                     f"match the experimental one well enough to place the "
                     f"ligand")
        lig = [l[:30] + "".join(f"{v:8.3f}" for v in (
                   Rp @ np.array([float(l[30:38]), float(l[38:46]),
                                  float(l[46:54])]) + tp)) + l[54:]
               for l in lig]

    # ------------------------------------- 2c. verify the MOVE, not the POSE
    #
    # A weak or non-binding peptide is a RESULT, not an error. Nothing here may
    # reject a peptide for sitting outside the pocket, and nothing anywhere in
    # the workflow pulls it towards one: the transform above is rigid, so the
    # peptide keeps exactly the pose the co-folding model predicted relative to
    # the receptor.
    #
    # What must still be caught is the transform going wrong, which looks like
    # a peptide hundreds of angstroms away. The two are told apart by checking
    # the OPERATION rather than the outcome: a rigid transform preserves
    # ligand-receptor distance, so the closest approach measured against the
    # PREDICTED receptor and against the EXPERIMENTAL one must agree to within
    # the difference between those two structures. That test is indifferent to
    # whether the peptide binds well, badly, or not at all.
    def closest(lig_lines, rec_lines):
        a1 = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                       for l in lig_lines])
        a2 = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                       for l in rec_lines])
        return float(np.min(np.linalg.norm(a1[:, None, :] - a2[None, :, :],
                                           axis=2)))

    gap = closest(lig, rec)
    if a.pose:
        pose_lines = atoms(a.pose)
        pred_rec = [l for l in pose_lines if l[21] == a.pose_chain]
        pred_lig = [l for l in pose_lines if l[21] != a.pose_chain]
        gap_pred = closest(pred_lig, pred_rec)
        print(f"closest ligand-receptor approach: {gap:.2f} A "
              f"(as predicted: {gap_pred:.2f} A)")
        if abs(gap - gap_pred) > 10.0:
            sys.exit(
                f"FATAL: closest approach changed from {gap_pred:.1f} A to "
                f"{gap:.1f} A across a RIGID transform, which is impossible - "
                f"the ligand was not placed in the receptor's frame.\n"
                f"       This is a coordinate error, not a weak binder.")
    else:
        print(f"closest ligand-receptor approach: {gap:.2f} A")

    # Recorded as data for Stage 7 to interpret, never used to gate the run.
    # A peptide parked on the surface is exactly the negative result the
    # benchmark needs in order to mean anything.
    (a.outdir / "pose_geometry.tsv").write_text(
        "metric\tvalue_A\n"
        f"closest_ligand_receptor_approach\t{gap:.3f}\n"
        + (f"closest_as_predicted\t{gap_pred:.3f}\n" if a.pose else ""))
    if gap > 5.0:
        print(f"NOTE: nearest ligand atom is {gap:.1f} A from the receptor - "
              f"this peptide is not in close contact. Recorded and carried "
              f"forward; it is a finding, not a failure.")

    lig_fixed = [l[:17] + f"{a.ligand_resname:>3}" + " L" + f"{1:>4}" + l[26:]
                 for l in lig]
    print(f"Ligand: {len(lig_fixed)} atoms as {a.ligand_resname} "
          f"({len(lig_ca)} CA)")

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
