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
    ap.add_argument("--force-gaff2", action="store_true",
                    help="Parameterise as a GAFF2 small molecule even if every "
                         "residue is canonical. Only for comparison.")
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------ complete the chemistry
    # Co-folding output carries no C-terminal OXT: the backbone carbonyl is
    # then read as an ALDEHYDE rather than a carboxylate, so the peptide is
    # short one oxygen and its net charge is wrong by +1. PDBFixer adds the
    # missing heavy atoms, including OXT, before protonation.
    #
    # This does not apply to ligands taken from experiment - DAMGO genuinely
    # ends in an alcohol - so the check below reports what was added rather
    # than assuming anything.
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    fixed_path = a.outdir / "ligand_complete.pdb"
    fixer = PDBFixer(filename=str(a.ligand))
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # do not build loops into a ligand
    fixer.findMissingAtoms()
    added = {str(k): [x.name for x in v] for k, v in fixer.missingAtoms.items()}
    fixer.addMissingAtoms()
    with open(fixed_path, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    # Count heavy atoms before and after rather than trusting missingAtoms.
    # Terminal atoms are tracked in a SEPARATE attribute (missingTerminals), so
    # missingAtoms is empty when the only thing added is the C-terminal OXT -
    # and this printed "no missing heavy atoms" while silently adding it. That
    # is precisely the atom D33 is about, and D33 claims this step reports what
    # it added, so the report has to be derived from the structure itself.
    def _heavy(path):
        return sum(1 for l in path.read_text().splitlines()
                   if l.startswith(("ATOM", "HETATM"))
                   and (l[76:78].strip() or l[12:16].strip()[0]).upper() != "H")
    n_before, n_after = _heavy(a.ligand), _heavy(fixed_path)
    term = getattr(fixer, "missingTerminals", {}) or {}
    term = {str(k): list(v) for k, v in term.items()}
    if n_after != n_before or added or term:
        print(f"PDBFixer: heavy atoms {n_before} -> {n_after}"
              + (f"; missingAtoms {added}" if added else "")
              + (f"; terminals {term}" if term else ""))
    else:
        print(f"PDBFixer: no heavy atoms added ({n_before} unchanged)")
    source = fixed_path

    # ------------------------------------- canonical peptide? use ff19SB
    #
    # D10 chose whole-molecule GAFF2 because "at 513 Da DAMGO is
    # small-molecule sized". That reasoning does not extend: a 20-mer is around
    # 2200 Da and is almost entirely backbone, and GAFF2's torsions are not
    # trained on peptide backbones the way ff19SB's are. D10's own revisit
    # trigger names exactly this.
    #
    # It does not scale either. AM1-BCC failed on a 6-mer and a 7-mer in this
    # library while succeeding on a 10-mer, so sqm convergence is already
    # unreliable and sequence-dependent at the short end.
    #
    # A peptide of standard L-amino acids needs none of it: ff19SB covers every
    # residue, tleap builds it directly, and the result is both faster and more
    # accurate. GAFF2 remains for ligands that genuinely are not standard
    # peptides - DAMGO, with its D-Ala, N-methyl-Phe and C-terminal Gly-ol.
    CANON = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
             "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
             "TYR", "VAL"}
    src_lines = [l for l in fixed_path.read_text().splitlines()
                 if l.startswith(("ATOM", "HETATM"))]
    res_seq, seen = [], set()
    for l in src_lines:
        key = (l[21], l[22:27])
        if key not in seen:
            seen.add(key)
            res_seq.append(l[17:20].strip())
    noncanon = sorted(set(res_seq) - CANON)

    if noncanon:
        print(f"non-canonical residues present: {noncanon} -> GAFF2")
    elif a.force_gaff2:
        print(f"{len(res_seq)} canonical residues, but --force-gaff2 given")
    else:
        # ff19SB path. Hydrogens are stripped because tleap rebuilds them from
        # its own templates - keeping externally placed ones is what made the
        # receptor fail with "Atom .R<NMET 65>.A<H 20> does not have a type".
        # Residues are renumbered from 901 so a residue-range mask can separate
        # ligand from receptor later; LIG-as-one-residue does not exist here.
        out, ridx, last = [], 900, None
        for l in src_lines:
            el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
            if el == "H":
                continue
            key = (l[21], l[22:27])
            if key != last:
                last, ridx = key, ridx + 1
            out.append(l[:21] + "L" + f"{ridx:>4}" + " " + l[27:])
        pep = a.outdir / "ligand_peptide.pdb"
        pep.write_text("\n".join(out) + "\nTER\nEND\n")

        # Charge from composition at pH 7.4, matching what tleap will build
        # with default residue names: Asp/Glu -1, Lys/Arg +1, His neutral,
        # plus a free N-terminal ammonium and C-terminal carboxylate.
        q = (sum(1 for r in res_seq if r in ("LYS", "ARG"))
             - sum(1 for r in res_seq if r in ("ASP", "GLU")) + 1 - 1)
        (a.outdir / "ligand_charge.txt").write_text(f"{q}\n")
        (a.outdir / "params_mode.txt").write_text("ff19SB\n")
        print(f"{len(res_seq)} canonical residues -> ff19SB, no antechamber")
        print(f"  sequence: {'-'.join(res_seq)}")
        print(f"  residues renumbered 901-{ridx}, {len(out)} heavy atoms")
        print(f"  net charge {q:+d} at pH {a.ph}")
        print(f"Wrote {pep}")
        return

    # ------------------------------------------------ protonate at pH
    # obabel IGNORES -p when -h is also given, which silently produces a
    # NEUTRAL N-terminus. For an opioid peptide that removes the ammonium that
    # forms the salt bridge with D147 - the defining interaction of binding.
    prot = a.outdir / "ligand_ph.pdb"
    r = subprocess.run(["obabel", str(source), "-O", str(prot),
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

    (a.outdir / "params_mode.txt").write_text("gaff2\n")
    (a.outdir / "ligand_charge.txt").write_text(f"{a.net_charge}\n")

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
