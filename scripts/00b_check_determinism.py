#!/usr/bin/env python3
"""Stage 0 determinism check.

The build instructions require Stage 0 to be deterministic - running twice must
produce byte-identical output - and require that to be verified and stated.
This re-runs the preparation into a scratch directory and compares SHA-256
digests against the recorded ones, so the claim is a checkable artefact rather
than a sentence in a log.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import subprocess
import sys
import tempfile


def sha256(p: pathlib.Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", required=True, type=pathlib.Path)
    ap.add_argument("--reference", required=True, type=pathlib.Path)
    ap.add_argument("--ph", type=float, default=7.4)
    ap.add_argument("--keep-chain", default="R")
    ap.add_argument("--extract-ligands", default="")
    a = ap.parse_args()

    ref = json.loads(a.reference.read_text())
    here = pathlib.Path(__file__).resolve().parent

    with tempfile.TemporaryDirectory() as tmp:
        cmd = [sys.executable, str(here / "00_receptor.py"),
               "--cif", str(a.cif.resolve()), "--outdir", tmp,
               "--ph", str(a.ph), "--keep-chain", a.keep_chain]
        if a.extract_ligands:
            cmd += ["--extract-ligands", a.extract_ligands]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-2000:])
            print(r.stderr[-2000:])
            sys.exit("FATAL: re-run of Stage 0 failed")

        rerun = {p.name: sha256(p) for p in sorted(pathlib.Path(tmp).glob("*"))
                 if p.suffix in (".pdb", ".pqr", ".fasta")}

    print("Stage 0 determinism check")
    print(f"  input          : {a.cif}")
    print(f"  reference      : {a.reference}")
    print(f"  files compared : {len(ref)}\n")

    ok = True
    for name, digest in sorted(ref.items()):
        got = rerun.get(name)
        if got is None:
            print(f"  MISSING   {name}")
            ok = False
        elif got != digest:
            print(f"  DIFFERS   {name}")
            print(f"            recorded {digest}")
            print(f"            re-run   {got}")
            ok = False
        else:
            print(f"  identical {name}  {digest[:16]}...")

    extra = set(rerun) - set(ref)
    for name in sorted(extra):
        print(f"  EXTRA     {name} (not in the recorded set)")

    print()
    if not ok:
        sys.exit("DETERMINISM CHECK FAILED: Stage 0 output is not reproducible")
    print("PASS: Stage 0 is byte-identical across runs.")


if __name__ == "__main__":
    main()
