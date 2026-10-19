#!/usr/bin/env python3
"""Parse the five met-enkephalin pose-stability arms (partial, timed-out runs).

Only prod.trr survives (6 frames @ 5 ns) plus prod.edr (50 ps). Receptor-frame
peptide RMSD, COM displacement and pocket-contact persistence, PBC-corrected by
taking the peptide image nearest the t=0 pocket centre.
"""
import json, pathlib, sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.lib.distances import distance_array

ARMS = ["metenk_exp", "metenk_hdk", "metenk_chai", "metenk_ptx", "metenk_scr"]
ROOT = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp/results/04_md")
OUT  = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp/tmp/metenk_parse")
CUT  = 4.0

def nearest_image(pos, ref_pt, box):
    """Shift pos by whole box vectors so its COM is nearest ref_pt."""
    L = box[:3]
    com = pos.mean(0)
    shift = np.round((ref_pt - com) / L) * L
    return pos + shift

summary, per_frame, contacts_out = [], [], []
for arm in ARMS:
    d = ROOT / arm
    u = mda.Universe(str(d / "min.tpr"), str(d / "prod.trr"))
    pep   = u.select_atoms("segid seg_7_system3 and not name H*")
    pepall= u.select_atoms("segid seg_7_system3")
    rec   = u.select_atoms("segid seg_0_system1 and name CA")
    rech  = u.select_atoms("segid seg_0_system1 and not name H*")
    seq = "".join(mda.lib.util.convert_aa_code(r.resname) for r in pepall.residues)

    u.trajectory[0]
    box0 = u.dimensions
    ref_rec = rec.positions.copy()
    rc0 = ref_rec.mean(0)
    ref_pep = nearest_image(pep.positions.copy(), rc0, box0)
    ref_rech = rech.positions.copy()
    d0 = distance_array(ref_rech, ref_pep)
    pocket_mask = (d0 < 6.0).any(axis=1)
    pocket_res = sorted({int(rech[i].resid) for i in np.nonzero(pocket_mask)[0]})
    ref_contacts = {int(r) for r in pocket_res
                    if (distance_array(rech[rech.resids == r].positions, ref_pep) < CUT).any()}

    rows, contact_keep = [], {int(r): 0 for r in ref_contacts}
    for ts in u.trajectory:
        box = u.dimensions
        mob = rec.positions.copy()
        mob_c = mob.mean(0)
        R, rec_rmsd = align.rotation_matrix(mob - mob_c, ref_rec - rc0)
        pp = nearest_image(pep.positions.copy(), mob_c, box)
        pp = (pp - mob_c) @ R.T + rc0
        rmsd = float(np.sqrt(np.mean(np.sum((pp - ref_pep) ** 2, axis=1))))
        com  = float(np.linalg.norm(pp.mean(0) - ref_pep.mean(0)))
        rh = (rech.positions.copy() - mob_c) @ R.T + rc0
        dm = distance_array(rh, pp)
        ncon = int((dm < CUT).sum())
        now = {int(rech[i].resid) for i in np.nonzero((dm < CUT).any(axis=1))[0]}
        for r in contact_keep:
            if r in now: contact_keep[r] += 1
        jac = len(now & ref_contacts) / max(1, len(now | ref_contacts))
        rows.append(dict(arm=arm, time_ns=ts.time/1000.0, pep_rmsd_A=rmsd,
                         com_disp_A=com, rec_ca_rmsd_A=float(rec_rmsd),
                         n_atom_contacts=ncon, n_res_contacts=len(now),
                         jaccard_vs_t0=jac))
    per_frame += rows
    nf = len(rows)
    kept = [v / nf for v in contact_keep.values()]
    summary.append(dict(
        arm=arm, sequence=seq, frames=nf, ns_covered=rows[-1]["time_ns"],
        pep_rmsd_final_A=rows[-1]["pep_rmsd_A"],
        pep_rmsd_max_A=max(r["pep_rmsd_A"] for r in rows),
        pep_rmsd_mean_last3_A=float(np.mean([r["pep_rmsd_A"] for r in rows[-3:]])),
        com_disp_final_A=rows[-1]["com_disp_A"],
        rec_ca_rmsd_final_A=rows[-1]["rec_ca_rmsd_A"],
        n_pocket_res_t0=len(pocket_res), n_contact_res_t0=len(ref_contacts),
        contact_persistence_mean=float(np.mean(kept)) if kept else None,
        n_res_contacts_final=rows[-1]["n_res_contacts"],
        jaccard_final=rows[-1]["jaccard_vs_t0"]))
    for r, v in sorted(contact_keep.items()):
        contacts_out.append(dict(arm=arm, resid=r, auth_resid=r + 64,
                                 frac_frames=v / nf))
    print(f"{arm}: {seq} rmsd_final={rows[-1]['pep_rmsd_A']:.2f} com={rows[-1]['com_disp_A']:.2f}",
          file=sys.stderr)

def write_tsv(path, rows):
    keys = list(rows[0].keys())
    with open(path, "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")

write_tsv(OUT / "pose_stability_summary.tsv", summary)
write_tsv(OUT / "pose_stability_per_frame.tsv", per_frame)
write_tsv(OUT / "pose_contact_persistence.tsv", contacts_out)
print(json.dumps(summary, indent=1))
