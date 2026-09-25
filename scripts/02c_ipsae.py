#!/usr/bin/env python3
"""Stage 2c - compute ipSAE for the Chai-1 models (decision D21).

WHY: the instructions require ipSAE, and Chai-1 does not produce it. More
importantly, Stage 2 showed ipTM flat at 0.27-0.35 across a 32-fold affinity
range - it could not separate a memorised, correct pose (Met-enkephalin,
1.6 A cross-seed, 2.3 A from experiment) from two peptides the model plainly
could not place (4-11 A cross-seed).

ipSAE exists to fix exactly the defect that most likely explains that. ipTM
scores whole chains, so residue pairs far from the interface dilute it; with a
281-residue receptor against a 5-residue peptide, that dilution is close to
worst case. ipSAE counts only residue pairs with good interchain PAE and
rescales the TM d0 parameter to the number of residues actually at the
interface.

IMPLEMENTATION: the reference implementation is vendored (vendor/ipsae.py,
Dunbrack lab, MIT). Only a format adapter is written here. A home-rolled
scoring function that is subtly wrong would look entirely plausible in a
results table, which is the worst failure mode for this repository.

The adapter targets ipsae.py's Boltz path, which takes a PAE .npz plus a .cif -
the pair Chai-1 already produces.

Usage:
    02c_ipsae.py --cofold results/02_cofold --out results/02_cofold/ipsae.tsv
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import re
import shutil
import subprocess
import sys

import numpy as np

VENDOR = pathlib.Path("/scratch/kuhfeldr-Kuhfeld_temp/vendor/ipsae.py")
PAE_CUTOFF = 10
DIST_CUTOFF = 15


def squeeze_sample(a: np.ndarray, want_dim: int) -> np.ndarray:
    """Chai returns arrays with a leading sample axis; take sample 0."""
    a = np.asarray(a)
    while a.ndim > want_dim:
        a = a[0]
    return a


def run_one(seeddir: pathlib.Path, work: pathlib.Path) -> dict | None:
    meta = seeddir / "_cand_meta.npz"
    cifs = sorted(seeddir.glob("pred.model_idx_*.cif"))
    if not meta.exists() or not cifs:
        return None
    d = np.load(meta)
    pae = squeeze_sample(d["pae"], 2)
    plddt = squeeze_sample(d["plddt"], 1)

    work.mkdir(parents=True, exist_ok=True)
    tag = "model0"
    np.savez(work / f"pae_{tag}.npz", pae=pae.astype(np.float32))
    np.savez(work / f"plddt_{tag}.npz", plddt=plddt.astype(np.float32))
    shutil.copy(cifs[0], work / f"{tag}.cif")

    r = subprocess.run(
        [sys.executable, str(VENDOR), f"pae_{tag}.npz", f"{tag}.cif",
         str(PAE_CUTOFF), str(DIST_CUTOFF)],
        cwd=work, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"    ipsae failed ({r.returncode}): "
              f"{(r.stderr or r.stdout).strip().splitlines()[-1][:160]}")
        return None

    # ipsae writes <stem>_<pae>_<dist>.txt next to the inputs
    out = sorted(work.glob("*.txt"))
    if not out:
        print("    ipsae produced no output file")
        return None
    rows = []
    for line in out[0].read_text().splitlines():
        parts = line.split()
        if len(parts) > 6 and parts[0] not in ("Chn1", "#"):
            rows.append(parts)
    if not rows:
        return None
    header = [l.split() for l in out[0].read_text().splitlines()
              if l.startswith("Chn1")]
    cols = header[0] if header else []
    best = None
    for p in rows:
        rec = dict(zip(cols, p)) if cols else {}
        # keep the "max" row type if present, else the first
        if rec.get("Type", "") in ("max", "") and best is None:
            best = rec
    return best or dict(zip(cols, rows[0]))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cofold", required=True, type=pathlib.Path)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    a = ap.parse_args()
    if not VENDOR.exists():
        sys.exit(f"FATAL: vendored ipsae.py not found at {VENDOR}")

    results = []
    for pepdir in sorted(p for p in a.cofold.iterdir()
                         if p.is_dir() and not p.name.startswith("_")):
        for sd in sorted(d for d in pepdir.glob("seed*") if d.is_dir()):
            work = sd / "_ipsae"
            rec = run_one(sd, work)
            if rec is None:
                continue
            row = {"peptide_id": pepdir.name, "seed": sd.name}
            row.update({k: v for k, v in rec.items()})
            results.append(row)
            print(f"  {pepdir.name} {sd.name}: "
                  f"ipSAE={rec.get('ipSAE','?')}  ipTM={rec.get('ipTM_af','?')}  "
                  f"n0res={rec.get('n0res','?')}")

    if not results:
        sys.exit("FATAL: no ipSAE results produced")

    keys = sorted({k for r in results for k in r},
                  key=lambda k: (k not in ("peptide_id", "seed"), k))
    with a.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        for r in results:
            w.writerow(r)
    print(f"\nWrote {a.out}  ({len(results)} rows)")


if __name__ == "__main__":
    main()
