#!/usr/bin/env python3
"""Stage 3b - GAFF2 / AM1-BCC parameters for a peptide ligand.

The ligand is treated as a SINGLE MOLECULE rather than as residues (D10). That
avoids a mixed-force-field junction mid-peptide, and means no non-standard
residue library is needed: GAFF2 covers arbitrary organic molecules, including
the N-methylated amide and the C-terminal alcohol in DAMGO.

TWO TRAPS HANDLED.

1. Atom names must be UNIQUE within the residue. A peptide assembled from
   several source residues contributes several atoms called N, CA, C and O;
   tleap then keeps "the first occurrence and ignores the rest" and the
   connectivity silently collapses.

2. The ligand must be ONE residue. antechamber's -rn renames residues but does
   not merge them, and tleap instantiates the unit once per residue - which
   produced eight copies of DAMGO and a system charge of +7.976 instead of +1.

Both are asserted, not assumed.
"""
from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ligand", required=True, type=pathlib.Path)
    ap.add_argument("--outdir", required=True, type=pathlib.Path)
    ap.add_argument("--net-charge", type=int, default=None,
                    help="expected formal charge at the given pH. When given "
                         "it is ASSERTED, which is how a wrong input gets "
                         "caught; when omitted the value obabel computes is "
                         "used and reported.")
    ap.add_argument("--resname", default="LIG")
    ap.add_argument("--ph", type=float, default=7.4)
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------ protonate at pH
    # obabel IGNORES -p when -h is also given, which silently produces a
    # NEUTRAL N-terminus. For an opioid peptide that removes the ammonium that
    # forms the salt bridge with D147 - the defining interaction of binding.
    prot = a.outdir / "ligand_ph.pdb"
    r = subprocess.run(["obabel", str(a.ligand), "-O", str(prot),
                        "-p", str(a.ph)], capture_output=True, text=True)
    if r.returncode != 0 or not prot.exists():
        sys.exit(f"FATAL: obabel failed: {r.stderr[-400:]}")

    from openbabel import pybel
    m = next(pybel.readfile("pdb", str(prot)))
    print(f"formula {m.formula}   charge {m.charge}   atoms {len(m.atoms)}")
    if a.net_charge is None:
        a.net_charge = m.charge
        print(f"net charge not specified; using computed {m.charge:+d} at pH {a.ph}")
    elif m.charge != a.net_charge:
        sys.exit(f"FATAL: protonated charge {m.charge}, expected "
                 f"{a.net_charge}. Check the pH and the expected state.")

    # ------------------------------------------------ unique names, one residue
    lines = [l for l in prot.read_text().splitlines()
             if l.startswith(("ATOM", "HETATM"))]
    counts: dict[str, int] = {}
    names = []
    for l in lines:
        el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        counts[el] = counts.get(el, 0) + 1
        names.append(f"{el}{counts[el]}")
    if len(set(names)) != len(names):
        sys.exit("FATAL: generated atom names are not unique")

    uniq = a.outdir / "ligand_unique.pdb"
    uniq.write_text("\n".join(
        l[:12] + f"{n:<4}" + l[16:17] + f"{a.resname:>3}" + " L" + "   1" + l[26:]
        for l, n in zip(lines, names)) + "\nEND\n")
    print(f"unique names: {len(set(names))}; collapsed to one residue")

    # ------------------------------------------------ antechamber + parmchk2
    r = subprocess.run(
        ["antechamber", "-i", "ligand_unique.pdb", "-fi", "pdb",
         "-o", "ligand.mol2", "-fo", "mol2", "-c", "bcc",
         "-nc", str(a.net_charge), "-at", "gaff2",
         "-rn", a.resname, "-s", "2", "-pf", "y"],
        cwd=a.outdir, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:])
        sys.exit("FATAL: antechamber failed")

    subprocess.run(["parmchk2", "-i", "ligand.mol2", "-f", "mol2",
                    "-o", "ligand.frcmod", "-s", "gaff2"],
                   cwd=a.outdir, check=True)

    # ------------------------------------------------ verify
    mol2 = (a.outdir / "ligand.mol2").read_text()
    block = mol2.split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0].strip().splitlines()
    m2names = [l.split()[1] for l in block]
    rids = {l.split()[6] for l in block}
    q = sum(float(l.split()[-1]) for l in block)
    attn = [l for l in (a.outdir / "ligand.frcmod").read_text().splitlines()
            if "ATTN" in l]

    print(f"\nmol2: {len(block)} atoms, {len(set(m2names))} unique names, "
          f"{len(rids)} residue(s), charge {q:+.4f}")
    print(f"frcmod ATTN lines: {len(attn)}")
    if len(set(m2names)) != len(m2names):
        sys.exit("FATAL: mol2 contains duplicate atom names")
    if len(rids) != 1:
        sys.exit(f"FATAL: mol2 has {len(rids)} residues; tleap would build "
                 f"that many copies of the ligand")
    if abs(q - a.net_charge) > 0.01:
        sys.exit(f"FATAL: charges sum to {q:+.4f}, expected {a.net_charge:+d}")
    if attn:
        print("NOTE: ATTN lines mark parameters parmchk2 guessed by analogy; "
              "they are the first place to look if Stage 6 misbehaves.")
    print("PASS")


if __name__ == "__main__":
    main()
