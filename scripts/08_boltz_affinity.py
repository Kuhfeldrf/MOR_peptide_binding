#!/usr/bin/env python3
"""Generate Boltz-2 affinity-prediction inputs for the MOR peptide benchmark.

WHY BOLTZ-2, AND WHAT THIS CORRECTS. This project established that co-folding
CONFIDENCE metrics carry no affinity signal: on 36 peptides with ChEMBL-cited Ki
spanning 4.43 log units, with length decorrelated from label (rho=+0.021) and 88%
power, iLIS scored rho=+0.196 and ipTM +0.216, both non-significant, and seven
metrics failed together. That is a result about *confidence* metrics.

Boltz-2 (Passaro et al., bioRxiv 2025.06.14.659707) is different in kind: it has
a head TRAINED on binding affinity, reports Pearson 0.62 on the held-out FEP+
benchmark - comparable to OpenFE FEP at ~1000x lower cost - and won the CASP16
affinity challenge across 140 complexes. It is the contemporary method for this
task and it had not been tested here.

TWO HONEST CAVEATS, both of which shape how the result must be read:

1. OFF-DISTRIBUTION. Boltz-2's affinity head accepts only small-molecule ligands
   (SMILES/CCD, <=128 heavy atoms) - it explicitly rejects protein chains. Short
   peptides are therefore supplied AS SMILES, which works mechanically (34/36 of
   our set fits) but is not what the head was trained on. The published 0.62 was
   on small molecules. A weak result here may say more about domain shift than
   about the model.

2. IC50 vs Ki. Boltz-2 outputs log10(IC50 in uM). Our labels are Ki. These differ
   by an assay-dependent Cheng-Prusoff factor, Ki = IC50 / (1 + [L]/Kd_L), so
   absolute agreement is not meaningful and Pearson on raw values would be
   misleading. **SPEARMAN (rank) is the valid metric**, because a monotonic
   offset leaves ranks unchanged.

The 128-atom limit excludes 2 of 36 peptides (a 15-mer at 135 heavy atoms and a
17-mer at 153). The remaining 34 retain the full 4.43 log range and length stays
decorrelated (rho=+0.032, p=0.857), so the benchmark survives the filter.

Usage:
    08_boltz_affinity.py --inputs docs/boltz_affinity_inputs.tsv \
        --receptor results/00_receptor/mOR_clean.fasta \
        --msa <receptor a3m> --outdir <dir> [--limit N]
"""
from __future__ import annotations

import argparse
import pathlib


def read_fasta(p: pathlib.Path) -> str:
    return "".join(l.strip() for l in p.read_text().splitlines()
                   if not l.startswith(">"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True, type=pathlib.Path)
    ap.add_argument("--receptor", required=True, type=pathlib.Path)
    ap.add_argument("--msa", required=True, type=pathlib.Path,
                    help="receptor a3m; Boltz needs an MSA per protein chain")
    ap.add_argument("--outdir", required=True, type=pathlib.Path)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()

    import csv
    receptor = read_fasta(a.receptor)
    if not a.msa.is_file():
        raise SystemExit(f"FATAL: MSA {a.msa} missing. Running MSA-free is the "
                         "condition that voided arm B v1 - refusing.")
    a.outdir.mkdir(parents=True, exist_ok=True)

    n = 0
    with a.inputs.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if a.limit and n >= a.limit:
                break
            name = row["name"]
            smiles = row["smiles"]
            # Chain A = receptor (protein), chain B = peptide as a SMILES ligand.
            # The affinity head requires a LIGAND chain as binder; a proteinChain
            # is rejected, which is why the peptide is not given as a sequence.
            yml = (
                "version: 1\n"
                "sequences:\n"
                "  - protein:\n"
                "      id: A\n"
                f"      sequence: {receptor}\n"
                f"      msa: {a.msa.resolve()}\n"
                "  - ligand:\n"
                "      id: B\n"
                f"      smiles: '{smiles}'\n"
                "properties:\n"
                "  - affinity:\n"
                "      binder: B\n"
            )
            (a.outdir / f"{name}.yaml").write_text(yml)
            n += 1
    print(f"wrote {n} Boltz-2 YAML inputs to {a.outdir}")
    print(f"  receptor: {len(receptor)} aa, chain A, MSA {a.msa.name}")
    print(f"  peptides: chain B as SMILES ligand, affinity binder = B")


if __name__ == "__main__":
    main()
