#!/usr/bin/env python3
"""How good is a co-folded receptor, measured the same way for every arm?

WHY THIS EXISTS. Arm B v1 was void because its receptor never folded (16.2 A
from experiment), and that was only discovered after 105 tasks had run. A gate
on receptor accuracy is now a precondition for reading any discrimination
statistic. But the first version of that gate was itself unsound: it matched CA
atoms BY INDEX, pairing prediction residue i with experimental residue i. The
experimental chain has 287 residues and the prediction 281, and
docs/decisions.md D-numbering records that 8F7Q is numbered +2 against the 6DDF
construct and that arm A's own 2.72 A needed `resid_offset: 64` - without it the
same comparison read 24.6 A. Index matching can therefore inflate RMSD by
misalignment alone.

THREE THINGS THIS FIXES:

1. CORRESPONDENCE IS DERIVED, NOT ASSUMED. The query and reference sequences are
   aligned with a real pairwise alignment, and only aligned, identical-residue
   pairs are superposed. No offsets are hand-set, so no numbering convention has
   to be trusted.

2. EVERY ARM IS MEASURED IDENTICALLY. Arm A's published 2.72 A came from a
   different script with a manual offset. Comparing arm B under a stricter
   measurement than arm A ever faced would be rigging the comparison, so this
   script re-measures both.

3. THE POCKET IS DEFINED STRUCTURALLY. Global CA RMSD over a 7TM bundle
   conflates fold quality with activation state - 8F7Q is an ACTIVE, Gi-bound
   structure (decisions.md D30 restricts the benchmark to active state), and an
   unconditioned prediction may land in a different state and score 4-5 A while
   being perfectly well folded. What peptide docking actually needs is a correct
   ORTHOSTERIC POCKET. So the pocket is defined as the reference receptor
   residues with any heavy atom within a cutoff of the reference peptide, read
   off 8F7Q itself - no residue list in any construct's numbering.

Reported per structure: all-CA RMSD, pocket-CA RMSD, and the number of pairs
each is computed over. Pocket RMSD is the decision-relevant number.

Usage:
    02f_receptor_accuracy.py --reference data/raw/8f7q.cif \\
        --ref-receptor-chain R --ref-peptide-chain P \\
        --pred <cif> [<cif> ...] [--pocket-cutoff 5.0] [--label arm]
"""
from __future__ import annotations

import argparse
import pathlib

import gemmi
import numpy as np

try:
    from Bio import Align
except ImportError:  # pragma: no cover
    Align = None


def chain_by_name(st: gemmi.Structure, name: str) -> gemmi.Chain:
    for ch in st[0]:
        if ch.name == name:
            return ch
    raise SystemExit(f"FATAL: chain {name!r} not in {st.name}")


def ca_records(ch: gemmi.Chain) -> list[tuple[str, np.ndarray]]:
    """(one-letter code, CA coord) for residues that have a CA."""
    out = []
    for res in ch:
        at = res.find_atom("CA", "*")
        if at is None:
            continue
        code = gemmi.find_tabulated_residue(res.name)
        letter = code.one_letter_code.upper() if code else "X"
        if not letter.isalpha():
            letter = "X"
        out.append((letter, np.array([at.pos.x, at.pos.y, at.pos.z])))
    return out


def align_indices(seq_a: str, seq_b: str) -> list[tuple[int, int]]:
    """Index pairs of aligned, IDENTICAL residues between two sequences."""
    if Align is None:
        raise SystemExit("FATAL: biopython required for sequence alignment")
    al = Align.PairwiseAligner()
    al.open_gap_score = -11.0
    al.extend_gap_score = -1.0
    al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    best = al.align(seq_a, seq_b)[0]
    pairs = []
    for (a0, a1), (b0, b1) in zip(*best.aligned):
        for k in range(a1 - a0):
            ia, ib = a0 + k, b0 + k
            if seq_a[ia] == seq_b[ib]:
                pairs.append((ia, ib))
    return pairs


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray, subset: np.ndarray | None = None):
    """Superpose on ALL pairs, then report RMSD over `subset` of those pairs.

    The superposition always uses the full aligned set: fitting on the pocket
    alone would hide pocket error by absorbing it into the fit.
    """
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    rot = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    Pa = (rot @ (P - pc).T).T + qc
    sel = slice(None) if subset is None else subset
    diff = Pa[sel] - Q[sel]
    n = len(diff)
    return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))) if n else float("nan"), n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True, type=pathlib.Path)
    ap.add_argument("--ref-receptor-chain", default="R")
    ap.add_argument("--ref-peptide-chain", default="P")
    ap.add_argument("--pred", required=True, nargs="+", type=pathlib.Path)
    ap.add_argument("--pocket-cutoff", type=float, default=5.0,
                    help="A from any reference-peptide heavy atom")
    ap.add_argument("--label", default="")
    a = ap.parse_args()

    ref = gemmi.read_structure(str(a.reference))
    ref.setup_entities()
    ref.remove_alternative_conformations()
    ref.remove_hydrogens()
    ref.remove_waters()
    rec = chain_by_name(ref, a.ref_receptor_chain)
    pep = chain_by_name(ref, a.ref_peptide_chain)

    # Pocket: reference receptor residues near the reference peptide. Defined
    # from the structure, so no numbering convention is relied on.
    pep_atoms = np.array([[at.pos.x, at.pos.y, at.pos.z]
                          for res in pep for at in res])
    ref_ca = ca_records(rec)
    ref_seq = "".join(c for c, _ in ref_ca)
    pocket_flag = []
    ci = 0
    for res in rec:
        if res.find_atom("CA", "*") is None:
            continue
        ra = np.array([[at.pos.x, at.pos.y, at.pos.z] for at in res])
        dmin = np.sqrt(((ra[:, None, :] - pep_atoms[None, :, :]) ** 2).sum(-1)).min()
        pocket_flag.append(dmin <= a.pocket_cutoff)
        ci += 1
    pocket_flag = np.array(pocket_flag)
    print(f"reference {a.reference.name}: receptor chain {a.ref_receptor_chain} "
          f"{len(ref_ca)} CA, peptide chain {a.ref_peptide_chain} {len(pep)} res")
    print(f"pocket = {int(pocket_flag.sum())} receptor residues within "
          f"{a.pocket_cutoff:.1f} A of the peptide\n")

    print(f"{'structure':52s} {'allCA':>8s} {'n':>5s} {'pocket':>8s} {'n':>4s}")
    rows = []
    for p in a.pred:
        st = gemmi.read_structure(str(p))
        st.setup_entities()
        st.remove_alternative_conformations()
        st.remove_hydrogens()
        # Receptor is the first chain in both arms' output.
        pred_ca = ca_records(st[0][0])
        pred_seq = "".join(c for c, _ in pred_ca)
        pairs = align_indices(pred_seq, ref_seq)
        if len(pairs) < 50:
            print(f"{p.name:52s}  ALIGNMENT FAILED ({len(pairs)} pairs)")
            continue
        P = np.array([pred_ca[i][1] for i, _ in pairs])
        Q = np.array([ref_ca[j][1] for _, j in pairs])
        pocket_sub = np.array([k for k, (_, j) in enumerate(pairs) if pocket_flag[j]])
        all_r, n_all = kabsch_rmsd(P, Q)
        pk_r, n_pk = kabsch_rmsd(P, Q, pocket_sub)
        print(f"{p.name:52s} {all_r:8.2f} {n_all:5d} {pk_r:8.2f} {n_pk:4d}")
        rows.append((all_r, pk_r))

    if rows:
        aa = np.array(rows)
        tag = f" [{a.label}]" if a.label else ""
        print(f"\nsummary{tag}: all-CA best {aa[:,0].min():.2f} "
              f"median {np.median(aa[:,0]):.2f} | pocket best {aa[:,1].min():.2f} "
              f"median {np.median(aa[:,1]):.2f}  (n={len(rows)} structures)")


if __name__ == "__main__":
    main()
