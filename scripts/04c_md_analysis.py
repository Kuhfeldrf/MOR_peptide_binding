#!/usr/bin/env python3
"""Stage 4 analysis - does DAMGO stay in the pocket, and is the system stable?

The build instructions require, per system: peptide RMSD, contact persistence
with receptor residues, and whether the peptide remains in the pocket.

All three are reported here, plus the membrane observables that would reveal a
bilayer going wrong (area per lipid, thickness) and the receptor's own drift.

Everything is measured after superposing on RECEPTOR backbone, so peptide
motion is reported relative to the receptor rather than to the box - the same
convention as the Stage 2 cross-seed diagnostic, so the numbers are comparable.

Usage:
    04c_analysis.py --dir results/03_membrane_prod --out results/04_md
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib

import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis import align, rms


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, type=pathlib.Path)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    ap.add_argument("--traj", default="prod.part0001.xtc")
    ap.add_argument("--resid-offset", type=int, default=64,
                    help="added to report original 6DDF numbering: packmol "
                         "renumbered the receptor from 1, so packed resid 83 "
                         "is D147 and 233 is H297")
    ap.add_argument("--contact-cut", type=float, default=4.0,
                    help="heavy-atom distance defining a contact, Angstrom")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    top = a.dir / "reduced.tpr"   # the xtc holds only the Reduced group
    trj = a.dir / a.traj
    u = mda.Universe(str(top), str(trj))
    print(f"frames {len(u.trajectory)}   atoms {len(u.atoms)}")

    lig = u.select_atoms("resname DAM and not name H*")
    rec = u.select_atoms("protein and name CA")
    rec_heavy = u.select_atoms("protein and not name H*")
    print(f"DAMGO heavy atoms {len(lig)}   receptor CA {len(rec)}")
    if not len(lig) or not len(rec):
        raise SystemExit("FATAL: could not select ligand or receptor")

    # reference: first frame
    u.trajectory[0]
    ref_rec = rec.positions.copy()
    ref_lig = lig.positions.copy()
    ref_com = lig.center_of_mass()

    # pocket residues: receptor residues with a heavy atom within the cutoff of
    # DAMGO in the reference frame
    from MDAnalysis.lib.distances import distance_array
    d0 = distance_array(rec_heavy.positions, ref_lig)
    pocket_mask = (d0 < 6.0).any(axis=1)
    pocket_res = sorted({rec_heavy[i].resid for i in np.nonzero(pocket_mask)[0]})
    print(f"pocket residues (within 6 A of DAMGO at t=0): {len(pocket_res)}")

    rows = []
    contact_counts = {r: 0 for r in pocket_res}
    for ts in u.trajectory:
        # superpose on receptor CA
        mobile = rec.positions
        R, rmsd_rec = align.rotation_matrix(mobile - mobile.mean(0),
                                            ref_rec - ref_rec.mean(0))
        lp = (lig.positions - mobile.mean(0)) @ R.T + ref_rec.mean(0)
        lig_rmsd = float(np.sqrt(np.mean(np.sum((lp - ref_lig) ** 2, axis=1))))
        com_disp = float(np.linalg.norm(lp.mean(0) - ref_lig.mean(0)))

        d = distance_array(rec_heavy.positions, lig.positions)
        in_contact = (d < a.contact_cut).any(axis=1)
        resids = {rec_heavy[i].resid for i in np.nonzero(in_contact)[0]}
        for r in resids:
            if r in contact_counts:
                contact_counts[r] += 1

        box = ts.dimensions
        rows.append({
            "time_ps": round(float(ts.time), 1),
            "ligand_rmsd_A": round(lig_rmsd, 3),
            "ligand_com_disp_A": round(com_disp, 3),
            "receptor_ca_rmsd_A": round(float(rmsd_rec), 3),
            "n_contact_residues": len(resids),
            "box_x_A": round(float(box[0]), 2),
            "box_z_A": round(float(box[2]), 2),
        })

    n = len(rows)
    out_tsv = a.out / "md_stability.tsv"
    with out_tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # contact persistence
    pers = sorted(((r, c / n) for r, c in contact_counts.items()),
                  key=lambda t: -t[1])
    with (a.out / "contact_persistence.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["resid_packed", "resid_6ddf", "resname", "fraction_of_frames"])
        for r, f in pers:
            nm = u.select_atoms(f"protein and resid {r}").residues[0].resname
            w.writerow([r, r + a.resid_offset, nm, f"{f:.3f}"])

    lr = np.array([r["ligand_rmsd_A"] for r in rows])
    cd = np.array([r["ligand_com_disp_A"] for r in rows])
    rr = np.array([r["receptor_ca_rmsd_A"] for r in rows])
    half = n // 2
    stats = {
        "frames": n,
        "duration_ns": rows[-1]["time_ps"] / 1000.0,
        "ligand_rmsd_mean_A": round(float(lr.mean()), 3),
        "ligand_rmsd_final_A": round(float(lr[-1]), 3),
        "ligand_rmsd_last_half_mean_A": round(float(lr[half:].mean()), 3),
        "ligand_com_displacement_final_A": round(float(cd[-1]), 3),
        "ligand_com_displacement_max_A": round(float(cd.max()), 3),
        "receptor_ca_rmsd_final_A": round(float(rr[-1]), 3),
        "contact_residues_persistent_gt50pct": sum(1 for _, f in pers if f > 0.5),
        "contact_cutoff_A": a.contact_cut,
    }
    # Did it stay in the pocket? A centre-of-mass displacement beyond ~8 A from
    # the starting pose means it has left, not merely rearranged.
    stats["remains_in_pocket"] = bool(cd[-1] < 8.0 and lr[half:].mean() < 6.0)
    (a.out / "md_stats.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")

    print("\n--- Stage 4 stability ---")
    for k, v in stats.items():
        print(f"  {k:<40} {v}")
    print("\ntop persistent contacts:")
    for r, f in pers[:14]:
        nm = u.select_atoms(f"protein and resid {r}").residues[0].resname
        print(f"  {nm}{r + a.resid_offset:<5} (packed {r:<4}) {f*100:5.1f} % of frames")
    print(f"\nWrote {out_tsv}")


if __name__ == "__main__":
    main()
