#!/usr/bin/env python3
"""Stage 0 - receptor preparation.

Takes the 6DDF mmCIF and emits a clean single-chain mu-opioid receptor plus the
experimentally resolved DAMGO ligand, with a log of every modification made.

Determinism is a hard requirement: running twice must produce byte-identical
output. No timestamps are written into any output file; the run log carries
them instead, and is excluded from the determinism hash.

Usage:
    00_receptor.py --cif data/raw/6ddf.cif --outdir results/00_receptor
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import subprocess
import sys
from datetime import datetime, timezone

import gemmi

RECEPTOR_CHAIN = "R"
LIGAND_CHAIN = "D"
DAMGO_COMPONENTS = ["TYR", "DAL", "GLY", "MEA", "ETA"]
INSPECT = [147, 297]


def log(lines: list[str], msg: str) -> None:
    print(msg, flush=True)
    lines.append(msg)


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def extract(cif: pathlib.Path, outdir: pathlib.Path, lines: list[str]):
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    model = st[0]
    present = [c.name for c in model]
    log(lines, f"Input: {cif} ({st.name})")
    log(lines, f"Chains present: {present}")

    # --- what the instructions expect to strip, versus what is actually here
    if RECEPTOR_CHAIN not in present:
        sys.exit(f"FATAL: receptor chain {RECEPTOR_CHAIN!r} absent from {cif}")
    stripped = [c for c in present if c not in (RECEPTOR_CHAIN, LIGAND_CHAIN)]
    log(lines, f"Stripped chains (Gi heterotrimer): {stripped}")
    log(lines,
        "NOTE: 6DDF contains no nanobody and no scFv16. The build instructions "
        "call for stripping them; they are absent from this structure, so that "
        "step is a no-op here and is recorded rather than silently skipped.")

    rec = model[RECEPTOR_CHAIN]
    nums = [r.seqid.num for r in rec]
    log(lines, f"Receptor chain {RECEPTOR_CHAIN}: {len(rec)} residues, "
               f"span {nums[0]}-{nums[-1]}")

    gaps = [(a, b) for a, b in zip(nums, nums[1:]) if b != a + 1]
    if gaps:
        for a, b in gaps:
            log(lines, f"  GAP: {a} -> {b} ({b - a - 1} residues missing)")
    else:
        log(lines, "  No chain breaks: residue numbering is continuous.")
        log(lines,
            "  CONSEQUENCE: ECL2 is fully resolved in 6DDF. The instructions "
            "call for modelling missing residues, particularly ECL2. There are "
            "none to model, so NO loop modelling is performed. This is a "
            "property of the structure, not a skipped step.")

    for t in INSPECT:
        hit = [r for r in rec if r.seqid.num == t]
        log(lines, f"  residue {t}: {hit[0].name if hit else 'ABSENT'}")

    # --- write receptor-only structure
    st_rec = gemmi.Structure()
    st_rec.name = "mOR"
    st_rec.spacegroup_hm = "P 1"
    m = gemmi.Model("1")
    m.add_chain(rec.clone())
    st_rec.add_model(m)
    st_rec.setup_entities()
    st_rec.remove_ligands_and_waters()
    raw_rec = outdir / "mOR_chainR_raw.pdb"
    st_rec.write_pdb(str(raw_rec))
    log(lines, f"Wrote receptor (pre-protonation): {raw_rec.name}")

    # --- write DAMGO, verbatim from experiment
    if LIGAND_CHAIN in present:
        lig = model[LIGAND_CHAIN]
        comps = [r.name for r in lig]
        log(lines, f"DAMGO chain {LIGAND_CHAIN}: {comps}")
        if comps != DAMGO_COMPONENTS:
            sys.exit(f"FATAL: chain {LIGAND_CHAIN} is {comps}, "
                     f"expected DAMGO {DAMGO_COMPONENTS}")
        n_atoms = sum(len(r) for r in lig)
        log(lines, f"  {n_atoms} heavy atoms (no hydrogens in the deposited model)")
        st_lig = gemmi.Structure()
        st_lig.name = "DAMGO"
        st_lig.spacegroup_hm = "P 1"
        ml = gemmi.Model("1")
        ml.add_chain(lig.clone())
        st_lig.add_model(ml)
        st_lig.setup_entities()
        out_lig = outdir / "damgo_ref.pdb"
        st_lig.write_pdb(str(out_lig))
        log(lines, f"Wrote DAMGO (experimental pose, not predicted): {out_lig.name}")
    else:
        sys.exit(f"FATAL: DAMGO chain {LIGAND_CHAIN!r} absent")

    return raw_rec


def protonate(raw: pathlib.Path, outdir: pathlib.Path, ph: float,
              lines: list[str]) -> pathlib.Path:
    """Assign protonation at the given pH with PROPKA via pdb2pqr."""
    pqr = outdir / "mOR_clean.pqr"
    pdb = outdir / "mOR_clean.pdb"
    cmd = ["pdb2pqr30", "--ff=AMBER", "--keep-chain",
           "--with-ph", str(ph), "--titration-state-method", "propka",
           "--pdb-output", str(pdb), str(raw), str(pqr)]
    log(lines, f"Protonation: {' '.join(cmd)}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        log(lines, r.stdout[-3000:])
        log(lines, r.stderr[-3000:])
        sys.exit(f"FATAL: pdb2pqr30 failed with code {r.returncode}")

    # Record the states of the residues the instructions name explicitly.
    st = gemmi.read_structure(str(pdb))
    chain = st[0][0]
    log(lines, f"Protonation assigned at pH {ph} (PROPKA):")
    for t in INSPECT:
        hit = [r for r in chain if r.seqid.num == t]
        if not hit:
            log(lines, f"  residue {t}: ABSENT after protonation")
            continue
        res = hit[0]
        names = {a.name for a in res}
        # Determine state from the hydrogens actually placed, not from the
        # residue name: pdb2pqr's --pdb-output keeps the generic HIS/ASP names,
        # so the name alone cannot distinguish HID/HIE/HIP or ASP/ASH.
        if res.name in ("HIS", "HID", "HIE", "HIP"):
            hd1 = "HD1" in names        # on ND1
            he2 = "HE2" in names        # on NE2
            if hd1 and he2:
                state, amber = "PROTONATED, +1 (both ring N)", "HIP"
            elif hd1:
                state, amber = "neutral, delta-protonated", "HID"
            elif he2:
                state, amber = "neutral, epsilon-protonated", "HIE"
            else:
                sys.exit(f"FATAL: H{t} has neither HD1 nor HE2; "
                         f"protonation could not be determined")
            log(lines, f"  H{t}: {amber} - {state}  (HD1={hd1}, HE2={he2})")
        elif res.name in ("ASP", "ASH"):
            hd2 = "HD2" in names        # on OD2
            if hd2:
                log(lines, f"  D{t}: ASH - PROTONATED, neutral  (HD2 present)")
            else:
                log(lines, f"  D{t}: ASP - deprotonated, -1  (no HD2)")
        elif res.name in ("GLU", "GLH"):
            he2 = "HE2" in names
            log(lines, f"  E{t}: {'GLH - PROTONATED, neutral' if he2 else 'GLU - deprotonated, -1'}")
        else:
            log(lines, f"  residue {t}: {res.name} ({len(names)} atoms)")

    return pdb


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", required=True, type=pathlib.Path)
    ap.add_argument("--outdir", required=True, type=pathlib.Path)
    ap.add_argument("--ph", type=float, default=7.4)
    a = ap.parse_args()

    a.outdir.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    log(lines, "=" * 70)
    log(lines, "Stage 0 - receptor preparation")
    log(lines, "=" * 70)

    raw = extract(a.cif, a.outdir, lines)
    clean = protonate(raw, a.outdir, a.ph, lines)

    # --- determinism record
    outputs = sorted(p for p in a.outdir.glob("*")
                     if p.suffix in (".pdb", ".pqr"))
    digests = {p.name: sha256(p) for p in outputs}
    log(lines, "")
    log(lines, "Output digests (SHA-256):")
    for n, d in digests.items():
        log(lines, f"  {d}  {n}")

    (a.outdir / "checksums.json").write_text(
        json.dumps(digests, indent=2, sort_keys=True) + "\n")

    # Log carries the timestamp; the structural outputs deliberately do not,
    # so that the determinism check compares chemistry, not clocks.
    hdr = (f"# Stage 0 run log\n"
           f"# generated {datetime.now(timezone.utc).isoformat()}\n"
           f"# input {a.cif}\n\n")
    (a.outdir / "receptor_prep.log").write_text(hdr + "\n".join(lines) + "\n")
    print(f"\nWrote {clean} and {len(outputs)} structural outputs.")


if __name__ == "__main__":
    main()
