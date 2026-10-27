#!/usr/bin/env python
"""Receptor CA RMSD of a co-folded prediction against a deposited structure.

WHY THIS FILE EXISTS. The gate (jobs/anchor_ptx_gate.sbatch) carried this
calculation inline, and it paired CA atoms BY LIST INDEX:

    n = min(len(P), len(Q)); P, Q = P[:n], Q[:n]

That is only valid when prediction index i and reference index i are the same
residue. For Protenix output they are not: Protenix renumbers its chains from 1
while 8F7Q chain R runs 66-352, so predicted residue 1 (MET) corresponds to
8F7Q residue 67 (MET) - a one-residue shift. Index pairing therefore superposed
each residue onto its neighbour one step down the helix. On job 242496 it
reported 4.47 A and FAILED the 4.00 A gate; the same five structures are
2.59-2.81 A once the correspondence is right. Residue identity at the offset
the old code used was 11/281 (3.9%); at the correct offset it is 276/281
(98.2%).

THE DESIGN RULE HERE: a misalignment must never again be reportable as a
folding result. This script establishes the correspondence FIRST and exits
non-zero without printing an RMSD if no offset reaches --min-identity. A gate
that cannot identify which residue is which has not measured anything, and
silence is the correct output. Reporting a number that is 1.8 A too high is
worse than reporting none, because it looks like evidence.

The THRESHOLD IS NOT CHANGED by this file and is passed in by the caller. The
fix is to the measurement, not to the criterion.

SCOPE OF THE OFFSET SCAN. Both chains here are contiguous single domains, so
the correspondence is one integer shift and a scan over shifts is exact,
inspectable, and needs no alignment parameters to argue about. It is NOT a
general sequence aligner: a prediction with internal insertions or deletions
relative to the reference would not be handled, and would correctly fail the
identity check rather than be silently mispaired.

RESIDUE IDENTITY IS NOT EXPECTED TO BE 100%. The receptor is built from 6DDF
and compared against 8F7Q; 5 of 281 positions differ (species/construct), which
is why --min-identity defaults to 0.90 and not 1.0.
"""
import argparse
import pathlib
import sys

import gemmi
import numpy as np


def ca_list(st, chain_name):
    """(seqid, one-letter, xyz) for every residue in the chain that has a CA."""
    out = []
    for res in st[0][chain_name]:
        atom = res.find_atom("CA", "*")
        if atom is None:
            continue
        info = gemmi.find_tabulated_residue(res.name)
        code = info.one_letter_code.upper() if info else "X"
        out.append((res.seqid.num, code, np.array([atom.pos.x, atom.pos.y, atom.pos.z])))
    return out


def kabsch(P, Q):
    """Superpose P onto Q; return (rmsd, per-residue deviation, P_aligned)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    Pa = (R @ (P - pc).T).T + qc
    dev = np.sqrt(np.sum((Pa - Q) ** 2, axis=1))
    return float(np.sqrt(np.mean(dev ** 2))), dev, Pa


def pairs_for_offset(n_pred, n_ref, off):
    """Index pairs (pred, ref) under a shift of `off`, clipped to both chains."""
    if off >= 0:
        n = min(n_pred, n_ref - off)
        return (list(range(n)), list(range(off, off + n))) if n > 0 else ([], [])
    n = min(n_pred + off, n_ref)
    return (list(range(-off, -off + n)), list(range(n))) if n > 0 else ([], [])


def establish_correspondence(pred, ref, min_identity, max_shift, min_overlap):
    """Find the offset maximising residue identity. Returns (pi, ei, frac, off).

    Exits non-zero if the best offset does not reach min_identity: see the
    design rule in the module docstring.
    """
    print("=== residue correspondence ===")
    print(f"  prediction: {len(pred)} CA, numbered {pred[0][0]}-{pred[-1][0]}")
    print(f"  reference:  {len(ref)} CA, numbered {ref[0][0]}-{ref[-1][0]}")
    print(f"  pred seq[:40]: {''.join(r[1] for r in pred[:40])}")
    print(f"  ref  seq[:40]: {''.join(r[1] for r in ref[:40])}")

    scored = []
    for off in range(-max_shift, max_shift + 1):
        pi, ei = pairs_for_offset(len(pred), len(ref), off)
        if len(pi) < min_overlap:
            continue
        ident = sum(1 for a, b in zip(pi, ei) if pred[a][1] == ref[b][1])
        scored.append((ident / len(pi), ident, off, pi, ei))
    if not scored:
        sys.exit(f"FATAL: no offset gives at least {min_overlap} paired residues.")

    scored.sort(key=lambda t: -t[0])
    print(f"\n  best offsets by identity (scan +/-{max_shift}):")
    for frac, ident, off, pi, _ in scored[:3]:
        print(f"    offset {off:+3d}: {ident:3d}/{len(pi)} identical ({100 * frac:5.1f}%)")

    frac, ident, off, pi, ei = scored[0]
    if frac < min_identity:
        print(f"\n  runner-up {100 * scored[1][0]:.1f}%" if len(scored) > 1 else "")
        sys.exit(
            f"FATAL: best residue identity {100 * frac:.1f}% at offset {off:+d} is below "
            f"the {100 * min_identity:.0f}% required to trust the pairing.\n"
            "       NO RMSD IS REPORTED. The prediction and the reference could not be\n"
            "       put into residue correspondence, so nothing has been measured. Do not\n"
            "       read this as a folding failure - it is a comparison failure. Check\n"
            "       chain selection, then whether the two sequences are the same protein."
        )

    print(f"\n  ACCEPTED offset {off:+d}: {ident}/{len(pi)} identical ({100 * frac:.1f}%)")
    print(f"    pred {pred[pi[0]][1]}{pred[pi[0]][0]} <-> ref {ref[ei[0]][1]}{ref[ei[0]][0]}")
    mism = [(pred[a][1], pred[a][0], ref[b][1], ref[b][0])
            for a, b in zip(pi, ei) if pred[a][1] != ref[b][1]]
    if mism:
        print(f"    {len(mism)} differing position(s), expected for 6DDF vs 8F7Q: "
              + ", ".join(f"{m[0]}{m[1]}->{m[2]}{m[3]}" for m in mism[:12]))
    return pi, ei, frac, off


def pick_ligand_chain(st, receptor_chain, requested, max_len):
    """The peptide chain actually bound to this receptor copy.

    8F7Q holds two receptor copies (R, M), two peptides (P, Q) and the Gi
    heterotrimer. Taking every non-receptor chain as 'the ligand' pulls in the
    G protein and the second receptor and makes any pocket definition
    meaningless - which is exactly the error an earlier ad-hoc version of this
    analysis made, reporting a 4.60 A 'pocket' RMSD over 183 of 281 residues.
    The peptide is chosen by proximity to THIS receptor copy.
    """
    if requested:
        return requested
    recv = np.array([r[2] for r in ca_list(st, receptor_chain)])
    best, best_d = None, None
    for ch in st[0]:
        if ch.name == receptor_chain or not 2 < len(ch) <= max_len:
            continue
        pts = np.array([[a.pos.x, a.pos.y, a.pos.z] for r in ch for a in r])
        if not len(pts):
            continue
        d = float(np.min(np.linalg.norm(recv[:, None, :] - pts[None, :, :], axis=2)))
        if best_d is None or d < best_d:
            best, best_d = ch.name, d
    if best:
        print(f"  ligand chain {best} ({len(st[0][best])} res), {best_d:.2f} A "
              f"from receptor chain {receptor_chain}")
    return best


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pred-dir", required=True,
                   help="directory searched recursively for *_sample_*.cif")
    p.add_argument("--ref-cif", required=True)
    p.add_argument("--ref-chain", default="R")
    p.add_argument("--threshold", type=float, default=4.0,
                   help="PASS if the best sample is at or below this CA RMSD")
    p.add_argument("--min-identity", type=float, default=0.90)
    p.add_argument("--max-shift", type=int, default=25)
    p.add_argument("--min-overlap", type=int, default=100)
    p.add_argument("--ligand-chain", default=None,
                   help="reference chain defining the pocket; auto-detected if omitted")
    p.add_argument("--max-ligand-len", type=int, default=60)
    p.add_argument("--pocket-radius", type=float, default=10.0)
    p.add_argument("--baseline", type=float, default=None,
                   help="prior arm's RMSD, printed for context only")
    a = p.parse_args()

    cifs = sorted(pathlib.Path(a.pred_dir).glob("**/*_sample_*.cif"))
    if not cifs:
        sys.exit(f"FATAL: no *_sample_*.cif under {a.pred_dir}")
    print(f"{len(cifs)} predicted structure(s) under {a.pred_dir}\n")

    ref_st = gemmi.read_structure(a.ref_cif)
    ref_st.setup_entities()
    print(f"=== {pathlib.Path(a.ref_cif).name} chains ===")
    print("  " + ", ".join(f"{c.name}({len(c)})" for c in ref_st[0]))
    lig_chain = pick_ligand_chain(ref_st, a.ref_chain, a.ligand_chain, a.max_ligand_len)
    lig_pts = None
    if lig_chain:
        lig_pts = np.array([[at.pos.x, at.pos.y, at.pos.z]
                            for r in ref_st[0][lig_chain] for at in r])
    else:
        print("  WARNING: no ligand chain found; pocket RMSD will be skipped.")

    ref_only = gemmi.read_structure(a.ref_cif)
    ref_only.setup_entities()
    ref_only.remove_ligands_and_waters()
    ref = ca_list(ref_only, a.ref_chain)

    first = gemmi.read_structure(str(cifs[0]))
    first.setup_entities()
    pred_chain = first[0][0].name
    pi, ei, frac, off = establish_correspondence(
        ca_list(first, pred_chain), ref, a.min_identity, a.max_shift, a.min_overlap)
    Q = np.array([ref[b][2] for b in ei])

    print(f"\n=== receptor CA RMSD vs {pathlib.Path(a.ref_cif).stem} "
          f"chain {a.ref_chain} ({len(pi)} CA) ===")
    results = []
    for cif in cifs:
        st = gemmi.read_structure(str(cif))
        st.setup_entities()
        pred = ca_list(st, st[0][0].name)
        if len(pred) != len(ca_list(first, pred_chain)):
            sys.exit(f"FATAL: {cif.name} has {len(pred)} CA, expected "
                     f"{len(ca_list(first, pred_chain))}; correspondence not transferable.")
        P = np.array([pred[k][2] for k in pi])
        r, dev, _ = kabsch(P, Q)
        n_naive = min(len(pred), len(ref))
        r_naive, _, _ = kabsch(np.array([x[2] for x in pred[:n_naive]]),
                               np.array([x[2] for x in ref[:n_naive]]))
        print(f"  {cif.name}: {r:.2f} A   (index-paired, the old bug: {r_naive:.2f} A)")
        results.append((r, dev, cif.name))

    best_r, best_dev, best_name = min(results, key=lambda t: t[0])
    worst_r = max(r for r, _, _ in results)
    print(f"\n  best {best_r:.2f} A ({best_name}) | worst {worst_r:.2f} A"
          f" | median residue deviation {np.median(best_dev):.2f} A")

    if lig_pts is not None and len(lig_pts):
        pocket = [k for k, b in enumerate(ei)
                  if np.min(np.linalg.norm(lig_pts - ref[b][2], axis=1)) <= a.pocket_radius]
        print(f"\n=== orthosteric pocket: {len(pocket)} CA within "
              f"{a.pocket_radius:.0f} A of chain {lig_chain} ===")
        if pocket:
            st = gemmi.read_structure(str(pathlib.Path(a.pred_dir) / best_name)
                                      if (pathlib.Path(a.pred_dir) / best_name).exists()
                                      else str([c for c in cifs if c.name == best_name][0]))
            st.setup_entities()
            pred = ca_list(st, st[0][0].name)
            P = np.array([pred[k][2] for k in pi])
            print(f"  pocket RMSD, global superposition: "
                  f"{np.sqrt(np.mean(best_dev[pocket] ** 2)):.2f} A")
            pr, _, _ = kabsch(P[pocket], Q[pocket])
            print(f"  pocket RMSD, superposed on pocket: {pr:.2f} A")
            print(f"  pocket median deviation:           {np.median(best_dev[pocket]):.2f} A")

    print(f"\n  gate threshold {a.threshold:.2f} A"
          + (f" | prior arm baseline {a.baseline:.2f} A" if a.baseline else ""))
    if best_r <= a.threshold:
        print(f"\n  GATE PASSED - {best_r:.2f} A <= {a.threshold:.2f} A; the 7TM bundle is "
              "folded and the array may run.")
        return 0
    print(f"\n  GATE FAILED - {best_r:.2f} A > {a.threshold:.2f} A. Do NOT launch the array "
          "and do NOT\n  compute any ranking from this condition.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
