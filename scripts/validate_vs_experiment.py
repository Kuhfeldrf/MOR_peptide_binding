#!/usr/bin/env python3
"""Validate Chai-1's Met-enkephalin pose against experiment (8F7Q).

Cross-seed convergence shows the model is CONSISTENT, not that it is CORRECT.
This asks the separate question: is the converged pose the right pose?

Reference: 8F7Q, Gi-bound mu-opioid receptor with beta-endorphin. The peptide
chain is YGGFMTSEKSQTPLVTLFKNA; its first five residues are the "message"
sequence, identical to Met-enkephalin, occupying the orthosteric pocket.

RESIDUE NUMBERING IS DETERMINED BY SEQUENCE, NOT ASSUMED. 8F7Q is numbered +2
relative to the 6DDF-derived construct. An earlier version of this script
assumed the construct's own numbering and produced 8.9% sequence identity at
matched positions - every RMSD from that mapping was meaningless. The offset
is now fitted by maximising identity and asserted to exceed 90%.

Two frames are reported, because they answer different questions:
  * global  - superpose on all matched receptor CA. Includes any
              activation-state or loop differences.
  * pocket  - superpose on orthosteric pocket residues only. Asks whether the
              peptide sits correctly in the pocket, independent of the rest.
"""
import csv
import pathlib

import gemmi
import numpy as np

ROOT = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp")
CUTOFF = 10.0
N_MSG = 5


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def rmsd(A, B):
    return float(np.sqrt(np.mean(np.sum((A - B) ** 2, axis=1))))


def one(res):
    return gemmi.find_tabulated_residue(res.name).one_letter_code.upper()


ref = gemmi.read_structure(str(ROOT / "data/raw/8f7q.cif")); ref.setup_entities()
ref_ca, ref_aa = {}, {}
for r in ref[0]["R"]:
    at = r.find_atom("CA", "*")
    if at:
        ref_ca[r.seqid.num] = np.array([at.pos.x, at.pos.y, at.pos.z])
        ref_aa[r.seqid.num] = one(r)

ref_pep = []
for r in ref[0]["P"]:
    at = r.find_atom("CA", "*")
    if at:
        ref_pep.append([at.pos.x, at.pos.y, at.pos.z])
    if len(ref_pep) >= N_MSG:
        break
ref_pep = np.asarray(ref_pep)
msg = "".join(one(r) for r in list(ref[0]["P"])[:N_MSG])
assert msg == "YGGFM", f"expected YGGFM message, got {msg}"

pocket = sorted(n for n, c in ref_ca.items()
                if np.min(np.linalg.norm(ref_pep - c, axis=1)) <= CUTOFF)

seeds = sorted(d for d in (ROOT / "results/02_cofold/real__met_enkephalin").glob("seed*") if d.is_dir())
first = gemmi.read_structure(str(sorted(seeds[0].glob("pred.model_idx_*.cif"))[0]))
first.setup_entities()
pred_aa = [one(r) for r in first[0][0] if r.find_atom("CA", "*")]

best_off, best_pct = None, -1.0
for off in range(40, 100):
    m = t = 0
    for i, aa in enumerate(pred_aa):
        n = off + i
        if n in ref_aa:
            t += 1
            m += (ref_aa[n] == aa)
    if t >= 100:
        pct = 100.0 * m / t
        if pct > best_pct:
            best_off, best_pct = off, pct
print(f"Fitted numbering offset: {best_off} ({best_pct:.1f}% identity at matched positions)")
assert best_pct > 90.0, f"numbering alignment failed at {best_pct:.1f}% identity"
print(f"Orthosteric pocket: {len(pocket)} residues within {CUTOFF} A of YGGFM")

rows = []
for sd in seeds:
    cif = sorted(sd.glob("pred.model_idx_*.cif"))
    if not cif:
        continue
    st = gemmi.read_structure(str(cif[0])); st.setup_entities()
    pred_ca = {}
    i = 0
    for r in st[0][0]:
        at = r.find_atom("CA", "*")
        if at:
            pred_ca[best_off + i] = np.array([at.pos.x, at.pos.y, at.pos.z])
            i += 1
    pep = []
    for r in st[0][1]:
        at = r.find_atom("CA", "*")
        if at:
            pep.append([at.pos.x, at.pos.y, at.pos.z])
    pep = np.asarray(pep)[:N_MSG]

    gc = sorted(set(pred_ca) & set(ref_ca))
    Pg = np.array([pred_ca[n] for n in gc]); Qg = np.array([ref_ca[n] for n in gc])
    Rg, tg = kabsch(Pg, Qg)
    rec_rmsd = rmsd((Rg @ Pg.T).T + tg, Qg)
    pep_g = rmsd((Rg @ pep.T).T + tg, ref_pep)

    pc = [n for n in pocket if n in pred_ca]
    Pp = np.array([pred_ca[n] for n in pc]); Qp = np.array([ref_ca[n] for n in pc])
    Rp, tp = kabsch(Pp, Qp)
    pocket_rmsd = rmsd((Rp @ Pp.T).T + tp, Qp)
    pep_al = (Rp @ pep.T).T + tp
    pep_p = rmsd(pep_al, ref_pep)
    centroid = float(np.linalg.norm(pep_al.mean(0) - ref_pep.mean(0)))
    per_res = [float(np.linalg.norm(a - b)) for a, b in zip(pep_al, ref_pep)]

    rows.append({
        "seed": sd.name,
        "numbering_offset": best_off,
        "receptor_ca_aligned": len(gc),
        "receptor_rmsd_global_A": round(rec_rmsd, 3),
        "peptide_rmsd_global_A": round(pep_g, 3),
        "pocket_residues": len(pc),
        "pocket_rmsd_A": round(pocket_rmsd, 3),
        "peptide_rmsd_pocket_A": round(pep_p, 3),
        "centroid_offset_A": round(centroid, 3),
        "per_residue_dev_A": ";".join(f"{d:.2f}" for d in per_res),
    })
    print(f"  {sd.name}: receptor {rec_rmsd:5.2f} | pocket {pocket_rmsd:5.2f} | "
          f"peptide(global) {pep_g:5.2f} | peptide(pocket) {pep_p:5.2f} | "
          f"centroid {centroid:5.2f} A")

out = ROOT / "results/02_cofold/validation_vs_8f7q.tsv"
with out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    for r in rows:
        w.writerow(r)

rec = [r["receptor_rmsd_global_A"] for r in rows]
pk = [r["pocket_rmsd_A"] for r in rows]
pg = [r["peptide_rmsd_global_A"] for r in rows]
pp = [r["peptide_rmsd_pocket_A"] for r in rows]
print(f"\nreceptor RMSD (global) : mean {np.mean(rec):.2f} A")
print(f"pocket RMSD            : mean {np.mean(pk):.2f} A")
print(f"peptide RMSD (global)  : mean {np.mean(pg):.2f} A")
print(f"peptide RMSD (pocket)  : mean {np.mean(pp):.2f} A  "
      f"(min {min(pp):.2f}, max {max(pp):.2f})")
v = np.mean(pp)
print("\nVERDICT: Chai-1 " + ("RECOVERS the experimental pose" if v < 2.5
      else "APPROXIMATELY recovers it" if v < 5.0 else "DOES NOT recover it")
      + f" (pocket-frame mean {v:.2f} A).")
print(f"Wrote {out}")
