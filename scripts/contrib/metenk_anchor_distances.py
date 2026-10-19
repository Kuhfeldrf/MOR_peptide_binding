#!/usr/bin/env python3
"""The two anchoring interactions that make the pocket view interpretable.

1. The "message" salt bridge: the peptide's N-terminal ammonium against the
   D3.32 carboxylate (D147, packed 83). In met-enkephalin that amine belongs to
   Tyr1. In the scramble GMYGF it belongs to a glycine, and the tyrosine is
   buried at position 3 - so any salt bridge the scramble forms is not the
   native one even when the distance looks identical.

2. The tyrosine phenol OH against the W293 / H297 / Y326 cluster, which is what
   distinguishes the Protenix pose from the experimental one: both anchor the
   same ammonium, but they point the phenol differently.

Reads the common-reference viewing sets in md_viz/metenk, so distances are
directly comparable between arms.
"""
import csv
import pathlib

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

VIZ = pathlib.Path(__file__).resolve().parents[2] / "md_viz/metenk"
OUT = pathlib.Path(__file__).resolve().parents[2] / "results/04_md_metenk"
ARMS = ["metenk_exp", "metenk_ptx", "metenk_chai", "metenk_hdk", "metenk_scr"]
# packed numbering; auth = packed + 64 (6DDF)
ASP = 83          # D147, D3.32
ANCHORS = "229 233 262"   # W293, H297, Y326


def main():
    rows = []
    for arm in ARMS:
        u = mda.Universe(str(VIZ / f"{arm}_frame0.pdb"),
                         str(VIZ / f"{arm}_prod.xtc"))
        pep = u.select_atoms("resid 898:902")
        # the tyrosine wherever it sits in this arm's sequence
        tyr = [r.resid for r in pep.residues if r.resname == "TYR"]
        nterm = u.select_atoms("resid 898 and name N")
        asp = u.select_atoms(f"resid {ASP} and name OD1 OD2")
        oh = u.select_atoms(f"resid {tyr[0]} and name OH") if tyr else None
        anc = u.select_atoms(f"resid {ANCHORS} and not name H*")
        seq = "".join(mda.lib.util.convert_aa_code(r.resname)
                      for r in pep.residues)
        for ts in u.trajectory:
            rows.append(dict(
                arm=arm, sequence=seq, time_ns=ts.time / 1000.0,
                nterm_to_D147_A=float(distance_array(nterm.positions,
                                                     asp.positions).min()),
                tyr_OH_to_anchors_A=float(distance_array(oh.positions,
                                                         anc.positions).min()),
                tyr_position=tyr[0] - 897))

    with open(OUT / "anchor_distances.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.3f}" if isinstance(v, float) else v)
                        for k, v in r.items()})

    times = sorted({r["time_ns"] for r in rows})
    for field, title in [("nterm_to_D147_A", "N-terminal ammonium -> D147 carboxylate"),
                         ("tyr_OH_to_anchors_A", "Tyr phenol OH -> W293/H297/Y326")]:
        print(f"\n{title}  (A)")
        print(f"{'arm':12s} {'seq':7s} " + " ".join(f"{t:>6.0f}" for t in times))
        for arm in ARMS:
            a = [r for r in rows if r["arm"] == arm]
            print(f"{arm:12s} {a[0]['sequence']:7s} "
                  + " ".join(f"{r[field]:6.2f}" for r in a))
    print("\nwrote", OUT / "anchor_distances.tsv")


if __name__ == "__main__":
    main()
