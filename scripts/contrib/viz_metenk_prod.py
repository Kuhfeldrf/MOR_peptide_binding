#!/usr/bin/env python3
"""Build small receptor-aligned viewing sets from the surviving prod.trr frames.

Per arm: receptor + peptide only, PBC-corrected (peptide image nearest the
receptor), superposed on receptor CA so peptide motion is receptor-relative.
Writes <arm>_frame0.pdb + <arm>_prod.xtc (6 frames) and a peptide-only
multi-model PDB for overlay figures.
"""
import pathlib
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align

ARMS = ["metenk_exp", "metenk_hdk", "metenk_chai", "metenk_ptx", "metenk_scr"]
ROOT = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp/results/04_md")
OUT  = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp/tmp/metenk_viz")
OUT.mkdir(parents=True, exist_ok=True)

def nearest_image(pos, ref_pt, box):
    L = box[:3]
    return pos + np.round((ref_pt - pos.mean(0)) / L) * L

# Common reference for every arm: the experimental arm's frame-0 receptor CA.
# Aligning each arm to its own frame 0 would make the arms mutually
# incomparable, which is the whole point of the comparison.
_u0 = mda.Universe(str(ROOT / "metenk_exp" / "min.tpr"),
                   str(ROOT / "metenk_exp" / "prod.trr"))
_u0.trajectory[0]
REF_CA = _u0.select_atoms("segid seg_0_system1 and name CA").positions.copy()
REF_C = REF_CA.mean(0)

for arm in ARMS:
    d = ROOT / arm
    u = mda.Universe(str(d / "min.tpr"), str(d / "prod.trr"))
    rec = u.select_atoms("segid seg_0_system1")
    recca = u.select_atoms("segid seg_0_system1 and name CA")
    pep = u.select_atoms("segid seg_7_system3")
    keep = rec + pep
    nrec = len(rec)

    ref_ca, rc0 = REF_CA, REF_C

    frames = []
    for ts in u.trajectory:
        box = u.dimensions
        mob = recca.positions.copy()
        mc = mob.mean(0)
        R, _ = align.rotation_matrix(mob - mc, ref_ca - rc0)
        rp = (rec.positions.copy() - mc) @ R.T + rc0
        pp = (nearest_image(pep.positions.copy(), mc, box) - mc) @ R.T + rc0
        frames.append(np.vstack([rp, pp]))

    # a universe holding only the kept atoms, then write the aligned frames
    w = mda.Merge(keep)
    w.load_new(np.array(frames), order="fac")
    w.atoms.positions = frames[0]
    w.atoms.write(str(OUT / f"{arm}_frame0.pdb"))
    with mda.Writer(str(OUT / f"{arm}_prod.xtc"), n_atoms=len(w.atoms)) as W:
        for ts in w.trajectory:
            ts.time = ts.frame * 5000.0
            W.write(w.atoms)
    # peptide-only multi-model PDB for a static overlay
    pw = mda.Merge(pep)
    pw.load_new(np.array([f[nrec:] for f in frames]), order="fac")
    with mda.Writer(str(OUT / f"{arm}_pep_frames.pdb"), n_atoms=len(pw.atoms),
                    multiframe=True) as W:
        for ts in pw.trajectory:
            W.write(pw.atoms)
    print(f"{arm}: {len(w.atoms)} atoms, {len(frames)} frames")
