#!/usr/bin/env python3
"""How far is each predicted pose from the experimental one, and does MD close it?

The per-arm RMSD in pose_stability_summary.tsv is measured against each arm's
OWN first frame, so it answers "did this pose move" and not "is this pose
right". Those are different questions and the second is the one that decides
whether a pose generator is any good.

The viewing sets in md_viz/metenk are all superposed on one reference (the
experimental arm's frame-0 receptor CA), so peptide coordinates are directly
comparable between arms with no further fitting.

The scramble is GMYGF against YGGFM: same composition, different order, so
atom-wise RMSD to the experimental pose is meaningless and is reported as NaN.
Its centre of mass is still comparable.
"""
import csv
import pathlib

import numpy as np
import MDAnalysis as mda

VIZ = pathlib.Path(__file__).resolve().parents[2] / "md_viz/metenk"
OUT = pathlib.Path(__file__).resolve().parents[2] / "results/04_md_metenk"
ARMS = ["metenk_exp", "metenk_ptx", "metenk_chai", "metenk_hdk", "metenk_scr"]
REF = "metenk_exp"
SEL = "resid 898:902 and not name H*"


def frames(arm):
    u = mda.Universe(str(VIZ / f"{arm}_frame0.pdb"), str(VIZ / f"{arm}_prod.xtc"))
    sel = u.select_atoms(SEL)
    out = []
    for ts in u.trajectory:
        out.append((ts.time / 1000.0, sel.positions.copy()))
    return out, [a.name for a in sel]


ref_frames, ref_names = frames(REF)
rows = []
for arm in ARMS:
    fr, names = frames(arm)
    comparable = names == ref_names
    for (t, p), (_, e) in zip(fr, ref_frames):
        rows.append(dict(
            arm=arm, time_ns=t,
            rmsd_vs_exp_A=float(np.sqrt(((p - e) ** 2).sum(1).mean()))
            if comparable else float("nan"),
            com_dist_vs_exp_A=float(np.linalg.norm(p.mean(0) - e.mean(0))),
            atomwise_comparable=comparable))

with open(OUT / "pose_vs_experimental.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v)
                    for k, v in r.items()})

print(f"{'arm':14s} {'RMSD t=0':>9s} {'RMSD t=25':>10s} {'COM t=0':>8s} {'COM t=25':>9s}")
for arm in ARMS:
    a = [r for r in rows if r["arm"] == arm]
    print(f"{arm:14s} {a[0]['rmsd_vs_exp_A']:9.2f} {a[-1]['rmsd_vs_exp_A']:10.2f} "
          f"{a[0]['com_dist_vs_exp_A']:8.2f} {a[-1]['com_dist_vs_exp_A']:9.2f}")
print("\nwrote", OUT / "pose_vs_experimental.tsv")
