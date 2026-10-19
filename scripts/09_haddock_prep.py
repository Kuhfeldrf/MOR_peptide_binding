#!/usr/bin/env python3
"""Prepare a POCKET-DIRECTED HADDOCK3 run: met-enkephalin + length-matched
opioids vs the user-supplied non-opioid peptides.

WHAT IS NEW HERE VERSUS THE PRIOR HADDOCK WORK. The July 2026 run used
`ranair = true` - RANDOM ambiguous restraints, i.e. it docked over the whole
receptor surface. That is why its finding #6 reads the way it does: near-native
poses were found (i-RMSD 0.99 A, DockQ 0.64) but cluster rank 1 had 0/10 pocket
contacts in 20/20 runs. With surface-wide sampling that is expected behaviour,
not a scoring mystery. **Directing HADDOCK at the orthosteric pocket with real
AIRs has never been tested in this project.**

RESIDUE NUMBERING - the trap this project has hit twice. `receptor_clean.pdb`
is numbered 66-352 and is in **8F7Q numbering, +2 relative to the 6DDF/MD
construct** the pocket list came from. Verified by residue identity, not
assumed: 14/14 match at +2, 1/14 at +0. Feeding raw MD numbering would point the
restraints at the wrong residues and the run would look completely healthy.

  MD/6DDF : 144 147 148 151 217 218 219 221 236 296 297 300 318 322
  THIS pdb: 146 149 150 153 219 220 221 223 238 298 299 302 320 324

FAIRNESS OF THE RESTRAINTS. Every peptide gets IDENTICAL restraint treatment:
receptor active = the 14 pocket residues, peptide = all residues passive. The
peptide's Tyr1/Phe4 are deliberately NOT marked active - doing so would encode
the opioid pharmacophore into the restraints and rig the very comparison being
made. The only thing that differs between peptides is their chemistry.

LENGTH MATCHING. Everything in this project has turned out to read peptide
length (HADDOCK included - the prior classification result found "no signal
beyond peptide length and run-to-run noise"). The opioid set is therefore chosen
to match the non-opioid lengths pair-for-pair at 7, 8, 10 and 11 residues.
Met-enkephalin (5 aa) is included because it was explicitly asked for, but it
has no length partner and must not carry the conclusion on its own.

CONFORMER ENSEMBLES. Short peptides are flexible and HADDOCK's rigidbody stage
cannot change backbone conformation, so each peptide is supplied as a 3-model
ensemble from RDKit ETKDG + MMFF minimisation. The prior run used ESMFold
structures, which are not meaningful for a 5-11mer.

Usage:
    09_haddock_prep.py --outdir <dir> [--n-conf 3]
"""
from __future__ import annotations

import argparse
import pathlib

# (sequence, name, class, length). Opioids carry cited pKi at MOR where known.
PEPTIDES = [
    # --- opioid, length-matched to the non-opioids at 7/8/10/11 ---
    ("YGGFM",       "op_met_enkephalin",  "opioid",     8.70),
    ("YGGFLRR",     "op_YGGFLRR",         "opioid",     7.42),
    ("YGGFLRRI",    "op_YGGFLRRI",        "opioid",     7.45),
    ("YGGFLRRIR",   "op_YGGFLRRIR",       "opioid",     7.62),
    ("YGGFLRKYPK",  "op_alpha_neoendorphin", "opioid",  8.14),
    ("YGGFLRRIRPK", "op_YGGFLRRIRPK",     "opioid",     8.25),
    # --- user-supplied non-opioid comparators ---
    ("DPYKLRP",     "non_DPYKLRP",        "non_opioid", None),
    ("REKVLASS",    "non_REKVLASS",       "non_opioid", None),
    ("EGPKLVAS",    "non_EGPKLVAS",       "non_opioid", None),
    ("QNALIVRYTR",  "non_QNALIVRYTR",     "non_opioid", None),
    ("FKCRRWQWRMK", "non_lactoferricinB", "non_opioid", None),
]

# Pocket residues in THIS receptor's numbering (8F7Q frame, verified by identity)
POCKET = [146, 149, 150, 153, 219, 220, 221, 223, 238, 298, 299, 302, 320, 324]


def build_ensemble(seq: str, n_conf: int, out: pathlib.Path) -> int:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSequence(seq)
    if m is None:
        raise SystemExit(f"FATAL: RDKit could not build {seq}")
    m = Chem.AddHs(m)
    # Plain ETKDGv3 fails to embed peptides beyond ~7 residues (distance-geometry
    # struggles with the extended backbone). useRandomCoords plus more attempts
    # fixes it; the ladder below escalates rather than silently returning zero
    # conformers, which would have dropped peptides from the comparison.
    cids = []
    for attempt, (rand, tries) in enumerate(
            [(False, 100), (True, 500), (True, 5000), (True, 20000)], start=1):
        p = AllChem.ETKDGv3()
        p.randomSeed = 20260929 + attempt
        p.useRandomCoords = rand
        p.maxIterations = tries
        p.enforceChirality = True
        # Disable RMSD pruning. With pruning on, RDKit silently returns FEWER
        # conformers for peptides whose sampled conformations are similar -
        # REKVLASS got 1 while everything else got 3. An unequal ensemble size
        # handicaps that peptide in docking (fewer starting conformations =
        # fewer chances at a good pose), which would bias the very comparison
        # being made. Every peptide must get exactly n_conf.
        p.pruneRmsThresh = -1.0
        cids = list(AllChem.EmbedMultipleConfs(m, numConfs=n_conf, params=p))
        # Escalate until we have the FULL ensemble, not merely a non-empty one.
        # Breaking on `if cids` stopped at 1 conformer for REKVLASS and produced
        # an unequal ensemble that would have biased the comparison.
        if len(cids) >= n_conf:
            if attempt > 1:
                print(f"    (embedded {seq} on attempt {attempt}, "
                      f"useRandomCoords={rand})")
            break
    if not cids:
        raise SystemExit(f"FATAL: only {len(cids)} conformers for {seq} after 4 escalating attempts")
    AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
    m = Chem.RemoveHs(m)
    # HADDOCK reads a multi-MODEL PDB as a conformational ensemble.
    with out.open("w") as fh:
        for i, cid in enumerate(cids, 1):
            fh.write(f"MODEL     {i:4d}\n")
            block = Chem.MolToPDBBlock(m, confId=cid)
            for line in block.splitlines():
                if line.startswith(("ATOM", "HETATM")):
                    # force chain B, segid B - HADDOCK needs the second molecule distinct
                    line = line[:21] + "B" + line[22:]
                    fh.write(line.ljust(72)[:72] + "B   \n")
            fh.write("ENDMDL\n")
        fh.write("END\n")
    return len(cids)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, type=pathlib.Path)
    ap.add_argument("--n-conf", type=int, default=3)
    a = ap.parse_args()
    inp = a.outdir / "inputs"
    inp.mkdir(parents=True, exist_ok=True)

    counts = []
    print(f"{'name':24s} {'seq':14s} {'len':>3s} {'class':11s} conformers")
    for seq, name, cls, pki in PEPTIDES:
        n = build_ensemble(seq, a.n_conf, inp / f"{name}.pdb")
        counts.append((name, n))
        print(f"{name:24s} {seq:14s} {len(seq):3d} {cls:11s} {n}")

    bad = [(n, c) for n, c in counts if c != a.n_conf]
    if bad:
        raise SystemExit(f"FATAL: unequal ensembles {bad} - would bias the "
                         "comparison; every peptide must get exactly "
                         f"{a.n_conf} conformers")
    print(f"\nall {len(counts)} peptides have exactly {a.n_conf} conformers - "
          "ensembles are equal, so no peptide is handicapped")

    # Ambiguous Interaction Restraints: receptor pocket active, peptide passive.
    # Written directly rather than via haddock3-restraints so the fairness
    # constraint above is explicit and auditable in one file.
    lines = []
    for r in POCKET:
        lines.append(
            "assign ( resid %d and segid A )\n"
            "       (\n"
            "        ( resid 1:40 and segid B )\n"
            "       )  2.0 2.0 0.0\n" % r
        )
    (a.outdir / "ambig.tbl").write_text("\n".join(lines))
    print(f"\nwrote {a.outdir/'ambig.tbl'}: {len(POCKET)} receptor AIRs -> peptide (all residues passive)")
    print(f"pocket residues (this pdb's numbering): {','.join(map(str,POCKET))}")


if __name__ == "__main__":
    main()
