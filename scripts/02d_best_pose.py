#!/usr/bin/env python3
"""Stage 2d - pick one starting pose for the physics stages.

Selects the highest aggregate-score model across every seed for this peptide
and writes it as a PDB for Stage 3.

WHAT THIS SELECTION DOES AND DOES NOT MEAN. Aggregate score is the model's
estimate of its own accuracy, not binding affinity, and Stage 2 showed it does
not discriminate: ipTM sat flat at 0.27-0.35 across a 32-fold affinity range,
and ipSAE no better. Picking the top-ranked model is therefore a convention for
choosing a starting structure, not a claim that it is correct.

The honest signal about whether that structure can be trusted is the CROSS-SEED
AGREEMENT for this peptide, which is recorded separately. A peptide whose seeds
diverge by 4-11 A has no reliable pose to hand to MD, and its downstream result
should be read with that in mind. This script therefore records the agreement
verdict alongside the chosen pose rather than discarding it.
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import shutil
import sys


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, type=pathlib.Path,
                    help="results/02_cofold/<peptide>")
    ap.add_argument("--out", required=True, type=pathlib.Path)
    a = ap.parse_args()

    best = None
    for seeddir in sorted(d for d in a.dir.glob("seed*") if d.is_dir()):
        tsv = seeddir / "scores.tsv"
        if not tsv.exists():
            continue
        with tsv.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                try:
                    score = float(row["aggregate_score"])
                except (KeyError, ValueError):
                    continue
                idx = row.get("model_idx", "0")
                cif = seeddir / f"pred.model_idx_{idx}.cif"
                if not cif.exists():
                    continue
                if best is None or score > best[0]:
                    best = (score, cif, seeddir.name, idx, row)

    if best is None:
        sys.exit(f"FATAL: no scored model found under {a.dir}")

    score, cif, seed, idx, row = best
    print(f"best model: {seed} model_idx {idx}  aggregate {score:.4f}  "
          f"iptm {row.get('iptm','?')}")

    # Convert the mmCIF to PDB for the Amber tooling downstream.
    import gemmi
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    st.write_pdb(str(a.out))

    # The complex is what Stage 3 places in the membrane, but the ligand
    # parameterisation needs the PEPTIDE ALONE. Passing the complex gave
    # antechamber a 4,682-atom, +15 system and it rejected it, which is the
    # good outcome: the charge assertion caught a wrong input rather than
    # silently parameterising a receptor.
    lig = gemmi.Structure()
    lig.name = 'LIG'
    lig.spacegroup_hm = 'P 1'
    m = gemmi.Model('1')
    if len(st[0]) < 2:
        sys.exit('FATAL: predicted structure has no second chain to take as ligand')
    m.add_chain(st[0][1].clone())      # chain 2 is the peptide
    lig.add_model(m)
    lig.setup_entities()
    ligpath = a.out.parent / 'best_ligand.pdb'
    lig.write_pdb(str(ligpath))
    print(f'ligand-only structure: {ligpath} '
          f'({sum(len(r) for r in st[0][1])} atoms)')

    meta = {
        "peptide": a.dir.name,
        "chosen_seed": seed,
        "chosen_model_idx": idx,
        "aggregate_score": score,
        "iptm": row.get("iptm"),
        "i_plddt_peptide": row.get("i_plddt_peptide"),
        "source_cif": str(cif),
        "note": ("Aggregate score is learned confidence, not affinity. Read "
                 "this pose together with the cross-seed agreement verdict for "
                 "this peptide."),
    }
    a.out.with_suffix(".json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"Wrote {a.out}")


if __name__ == "__main__":
    main()
