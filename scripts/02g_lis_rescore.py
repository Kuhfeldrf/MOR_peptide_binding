#!/usr/bin/env python3
"""LIS / cLIS / iLIS rescoring of co-folded complexes, from saved PAE.

WHY. Every arm of this pilot has failed to separate opioid peptides from
non-opioid peptides using GLOBAL confidence metrics - ipTM, ipSAE, ranking_score
- or physical scores dominated by peptide size. The documented fix for exactly
this failure mode is the Local Interaction Score family:

  Kim, A.R. et al. "Enhanced Protein-Protein Interaction Discovery via
  AlphaFold-Multimer." bioRxiv 2024.02.19.580970.
  https://github.com/flyark/AFM-LIS

LIS restricts attention to the CONFIDENTLY PREDICTED INTERFACE rather than
averaging over a whole complex, and the authors report it "outperformed global
metrics in every pLDDT subgroup" and "uniquely rescued 47 well-established PPIs
missed by ipTM", with the gain concentrated in "flexible and small interfaces".
A 4-20mer peptide in a GPCR vestibule is precisely a small flexible interface,
and ipTM averaging over a 281-residue receptor is precisely the dilution LIS
was designed to avoid. This is a literature-documented metric, not a new one
invented to rescue a negative result.

DEFINITIONS, and an honest note on one ambiguity:
  LIA  = local interaction area = count of inter-chain residue pairs with
         PAE <= 12 A.
  LIS  = mean of the PAE values in the LIA after inversion to 0-1.
  cLIS = same, restricted to pairs ALSO in direct contact (Cb-Cb <= 8 A).
  iLIS = sqrt(LIS * cLIS).

  The paper states only that "PAE values within the LIA are inverted to a scale
  of 0 to 1" and that the AVERAGE is taken; neither it nor the repository README
  gives the arithmetic. This implementation uses the standard reading
  score = (12 - PAE) / 12, so PAE 0 -> 1.0 and PAE 12 -> 0.0. Recorded because
  a different normalisation would change absolute values (it would not change
  the RANKING, which is monotonic in mean PAE over the LIA either way).

  Both inter-chain submatrices (receptor->peptide and peptide->receptor) are
  computed and averaged, per the paper's "average ... at the interaction
  interfaces".

Cb is used where present and CA is substituted for glycine, which has no Cb.

Usage:
    02g_lis_rescore.py --cofold results/02_cofold --out docs/lis_scores.tsv
                       [--n-receptor 281] [--pae-cutoff 12] [--contact-cutoff 8]
"""
from __future__ import annotations

import argparse
import csv
import pathlib

import gemmi
import numpy as np

PAE_CUTOFF = 12.0
CONTACT_CUTOFF = 8.0


def cb_coords(st: gemmi.Structure) -> tuple[np.ndarray, np.ndarray]:
    """Cb coordinates (CA for Gly) for chain 0 and chain 1."""
    out = []
    for ch in list(st[0])[:2]:
        pts = []
        for res in ch:
            at = res.find_atom("CB", "*") or res.find_atom("CA", "*")
            if at is not None:
                pts.append([at.pos.x, at.pos.y, at.pos.z])
        out.append(np.asarray(pts, dtype=float))
    while len(out) < 2:
        out.append(np.zeros((0, 3)))
    return out[0], out[1]


def lis_from_block(block: np.ndarray, contact_mask: np.ndarray | None,
                   pae_cutoff: float) -> tuple[float, int, float, int]:
    """(LIS, LIA, cLIS, cLIA) for one inter-chain PAE submatrix."""
    lia = block <= pae_cutoff
    n_lia = int(lia.sum())
    lis = float(((pae_cutoff - block[lia]) / pae_cutoff).mean()) if n_lia else 0.0
    if contact_mask is None:
        return lis, n_lia, float("nan"), 0
    c = lia & contact_mask
    n_c = int(c.sum())
    clis = float(((pae_cutoff - block[c]) / pae_cutoff).mean()) if n_c else 0.0
    return lis, n_lia, clis, n_c


def score_one(pae: np.ndarray, n_rec: int, cb_r: np.ndarray, cb_p: np.ndarray,
              pae_cutoff: float, contact_cutoff: float) -> dict:
    rp = pae[:n_rec, n_rec:]          # receptor rows, peptide cols
    pr = pae[n_rec:, :n_rec]          # peptide rows, receptor cols

    cm_rp = None
    if len(cb_r) and len(cb_p):
        nr = min(len(cb_r), rp.shape[0])
        npep = min(len(cb_p), rp.shape[1])
        d = np.sqrt(((cb_r[:nr, None, :] - cb_p[None, :npep, :]) ** 2).sum(-1))
        cm_rp = np.zeros(rp.shape, dtype=bool)
        cm_rp[:nr, :npep] = d <= contact_cutoff

    lis1, lia1, clis1, clia1 = lis_from_block(rp, cm_rp, pae_cutoff)
    cm_pr = cm_rp.T if cm_rp is not None else None
    lis2, lia2, clis2, clia2 = lis_from_block(pr, cm_pr, pae_cutoff)

    lis = (lis1 + lis2) / 2.0
    clis = np.nanmean([clis1, clis2])
    # iLIS is undefined if either component is zero; report 0 rather than nan so
    # a complex with no confident contacts ranks last instead of dropping out.
    ilis = float(np.sqrt(lis * clis)) if (lis > 0 and clis > 0) else 0.0
    return {
        "LIS": round(lis, 5), "LIA": lia1 + lia2,
        "cLIS": round(float(clis), 5), "cLIA": clia1 + clia2,
        "iLIS": round(ilis, 5),
        "pae_interface_mean": round(float(np.mean([rp.mean(), pr.mean()])), 3),
        "pae_interface_min": round(float(min(rp.min(), pr.min())), 3),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cofold", required=True, type=pathlib.Path)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    ap.add_argument("--n-receptor", type=int, default=None,
                    help="receptor token count; inferred from the CIF if omitted")
    ap.add_argument("--pae-cutoff", type=float, default=PAE_CUTOFF)
    ap.add_argument("--contact-cutoff", type=float, default=CONTACT_CUTOFF)
    a = ap.parse_args()

    rows: list[dict] = []
    pepdirs = sorted(p for p in a.cofold.iterdir() if p.is_dir() and not p.name.startswith("_"))
    for pepdir in pepdirs:
        for seeddir in sorted(d for d in pepdir.glob("seed*") if d.is_dir()):
            meta = seeddir / "_cand_meta.npz"
            if not meta.is_file():
                continue
            pae_all = np.asarray(np.load(meta)["pae"])
            if pae_all.ndim == 2:
                pae_all = pae_all[None, ...]
            cifs = sorted(seeddir.glob("pred.model_idx_*.cif"))
            for idx in range(pae_all.shape[0]):
                pae = pae_all[idx]
                cb_r = cb_p = np.zeros((0, 3))
                n_rec = a.n_receptor
                if idx < len(cifs):
                    st = gemmi.read_structure(str(cifs[idx]))
                    st.setup_entities()
                    cb_r, cb_p = cb_coords(st)
                    if n_rec is None:
                        n_rec = len(cb_r)
                if n_rec is None or n_rec <= 0 or n_rec >= pae.shape[0]:
                    continue
                d = score_one(pae, n_rec, cb_r, cb_p, a.pae_cutoff, a.contact_cutoff)
                d.update(peptide_id=pepdir.name, seed=seeddir.name,
                         model_idx=idx, n_receptor=n_rec,
                         n_peptide=pae.shape[0] - n_rec)
                rows.append(d)
        print(f"  {pepdir.name}: {sum(1 for r in rows if r['peptide_id']==pepdir.name)} structures")

    if not rows:
        raise SystemExit("FATAL: no PAE matrices found; nothing rescored")

    cols = ["peptide_id", "seed", "model_idx", "n_receptor", "n_peptide",
            "LIS", "LIA", "cLIS", "cLIA", "iLIS",
            "pae_interface_mean", "pae_interface_min"]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {a.out}  ({len(rows)} structures, {len(pepdirs)} peptides)")


if __name__ == "__main__":
    main()
