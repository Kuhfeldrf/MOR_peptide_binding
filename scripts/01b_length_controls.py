#!/usr/bin/env python3
"""Build length-varied controls that decouple PEPTIDE LENGTH from the OPIOID MESSAGE.

WHY THE SCRAMBLED DECOYS ARE NOT ENOUGH
---------------------------------------
Scrambling preserves amino-acid composition, which is exactly what makes it a
good control for size and composition bias - and exactly what makes it a weak
control for BINDING.

The opioid "message" is an N-terminal tyrosine with a free, protonated
alpha-amine: the amine forms the salt bridge to Asp147 that the 50 ns DAMGO run
found in 100% of frames, and the phenol reaches His297. A scramble of YPFP into
PYPF still carries a tyrosine, a free N-terminal amine, both aromatics and the
identical hydrophobic surface. It moves Tyr off position 1 but removes none of
the chemistry.

So scrambles are plausibly WEAK BINDERS rather than non-binders, and an AUROC
near 0.5 against them is ambiguous: it may mean the score cannot discriminate,
or that the decoys genuinely bind about as well.

WHY LENGTH ALONE IS NOT ENOUGH EITHER
-------------------------------------
Beta-endorphin is 31 residues and binds muOR with high affinity - that is 8F7Q,
already used in this repository for Stage 2 validation. Its message is the
N-terminal YGGFM; the remainder extends out of the pocket. Long peptides are
therefore not automatically non-binders, and a control set built on length
alone would confound "too long to bind" with "lacks the pharmacophore".

THE DESIGN
----------
Two matched series at 10, 15 and 20 residues:

  NEGATIVES  real fragments of milk proteins that do NOT begin with tyrosine.
             Real digest products, so realistic screening candidates, and they
             lack the message.

  POSITIVES  N-terminal fragments of beta-endorphin, which DO begin YGGFM.
             Real, known-active, and length-matched to the negatives.

Length is varied in both arms and the message only in one, so the two are
separable. If the long positives score well and the long negatives do not, the
pipeline discriminates at length. If everything fails at 20 residues, that is a
capability limit rather than a discrimination result - a different finding, and
one this design can tell apart.

An internal tyrosine is NOT excluded from the negatives: the message requires
Tyr at position 1 specifically, and beta-endorphin itself carries an internal
Tyr at 27. Which negatives contain one is recorded so the question can be
asked later.
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import random
import sys
import urllib.request

# Milk proteins the benchmark peptides actually derive from, plus the whey
# proteins, so the negatives share provenance with the real peptides.
PROTEINS = {
    "P02662": ("alphas1-casein", "bovine"),
    "P02663": ("alphas2-casein", "bovine"),
    "P02666": ("beta-casein", "bovine"),
    "P02668": ("kappa-casein", "bovine"),
    "P02754": ("beta-lactoglobulin", "bovine"),
    "P00711": ("alpha-lactalbumin", "bovine"),
    "P05814": ("beta-casein", "human"),
}

# Beta-endorphin's precursor. The message sequence YGGFM marks its N-terminus.
POMC = ("P01189", "proopiomelanocortin", "human")

LENGTHS = (10, 15, 20)


def fetch(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url, timeout=30) as r:
        body = r.read().decode()
    seq = "".join(l.strip() for l in body.splitlines() if not l.startswith(">"))
    if not seq:
        sys.exit(f"FATAL: no sequence returned for {acc}")
    return seq


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=pathlib.Path)
    ap.add_argument("--per-length", type=int, default=6,
                    help="negatives per length, spread across proteins")
    ap.add_argument("--seed", type=int, default=0,
                    help="fixed so the control set is reproducible")
    a = ap.parse_args()
    rng = random.Random(a.seed)

    rows = []

    # ------------------------------------------------ positives: has message
    pomc = fetch(POMC[0])
    i = pomc.find("YGGFM")
    if i < 0:
        sys.exit("FATAL: YGGFM not found in POMC - cannot locate beta-endorphin")
    print(f"beta-endorphin starts at POMC residue {i + 1}")
    for L in LENGTHS:
        seq = pomc[i:i + L]
        if len(seq) < L:
            sys.exit(f"FATAL: POMC too short for a {L}-mer from {i}")
        rows.append({
            "peptide_id": f"pos__beta_endorphin_1_{L}",
            "sequence": seq,
            "peptide_name": f"beta-endorphin(1-{L})",
            "parent_protein": POMC[1],
            "species": POMC[2],
            "activity": "length_positive",
            "ic50_um": "",
            "reference": f"UniProt:{POMC[0]} residues {i+1}-{i+L}",
            "scramble_of": "",
            "n_term": seq[0],
            "has_message": "YES",
            "internal_tyr": "YES" if "Y" in seq[1:] else "NO",
        })

    # ------------------------------------------------ negatives: no message
    #
    # Sampled from real milk proteins, excluding any window that begins with
    # tyrosine. Signal peptides are skipped (the first 15 residues) because they
    # are cleaved and never appear in a digest.
    pool: dict[int, list[tuple]] = {L: [] for L in LENGTHS}
    for acc, (name, species) in PROTEINS.items():
        seq = fetch(acc)
        for L in LENGTHS:
            for s in range(15, len(seq) - L + 1):
                w = seq[s:s + L]
                if w[0] == "Y":
                    continue                     # would carry the message
                if "X" in w or "U" in w:
                    continue
                pool[L].append((acc, name, species, s, w))

    for L in LENGTHS:
        cands = pool[L]
        if len(cands) < a.per_length:
            sys.exit(f"FATAL: only {len(cands)} candidate {L}-mers")
        # One per protein first, so the set is not dominated by the longest
        # sequence, then fill at random from the remainder.
        by_prot: dict[str, list] = {}
        for c in cands:
            by_prot.setdefault(c[0], []).append(c)
        picked = [rng.choice(v) for v in by_prot.values()]
        rng.shuffle(picked)
        picked = picked[:a.per_length]
        while len(picked) < a.per_length:
            c = rng.choice(cands)
            if c not in picked:
                picked.append(c)
        for acc, name, species, s, w in picked:
            rows.append({
                "peptide_id": f"neg__{name.replace('-', '_')}_{species}_{s+1}_{L}",
                "sequence": w,
                "peptide_name": f"{name} ({species}) {s+1}-{s+L}",
                "parent_protein": name,
                "species": species,
                "activity": "length_negative",
                "ic50_um": "",
                "reference": f"UniProt:{acc} residues {s+1}-{s+L}",
                "scramble_of": "",
                "n_term": w[0],
                "has_message": "NO",
                "internal_tyr": "YES" if "Y" in w else "NO",
            })

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wtr.writeheader()
        wtr.writerows(rows)

    npos = sum(1 for r in rows if r["activity"] == "length_positive")
    nneg = len(rows) - npos
    print(f"\n{npos} positives (message present), {nneg} negatives (no message)")
    for L in LENGTHS:
        p = [r for r in rows if r["activity"] == "length_positive"
             and len(r["sequence"]) == L]
        n = [r for r in rows if r["activity"] == "length_negative"
             and len(r["sequence"]) == L]
        print(f"  {L:>2} aa: {len(p)} positive, {len(n)} negative")
    assert all(r["sequence"][0] != "Y" for r in rows
               if r["activity"] == "length_negative"), \
        "a negative begins with tyrosine"
    print(f"\nWrote {a.out}")


if __name__ == "__main__":
    main()
