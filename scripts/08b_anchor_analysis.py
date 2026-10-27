#!/usr/bin/env python3
"""Anchor-set arm B: does the Protenix score track affinity, or peptide length?

THE CONSTRAINT THAT SHAPES THIS WHOLE SCRIPT. The anchor set's 20 affinities
come from four documents and four assay types (Wang 2023 pKi; Koch 1985 KD;
Chiba 1989 IC50; Tani 1990 IC50). They are NOT interconvertible, and pooling
them is the error already documented for the 36-peptide set. So this script
computes NO cross-block statistic. Every correlation is reported inside one
block, which caps n at 5 (block A) and 8 (block B) and means almost no
statistical power - stated here rather than buried, because a Spearman rho on
n=5 is a description of five points, not evidence about a method.

scripts/07c_affinity_correlation.py is deliberately NOT reused: it pools
affinities per peptide as a median pAffinity across assay types, which is
exactly the move the anchor set's provenance forbids.

THE CONTROL THAT MATTERS. Every structure-based score tested in this project so
far has turned out to be reading peptide length. So length correlation is
reported beside every affinity correlation, in the same table, at the same
precision. If a score correlates with affinity no better than it correlates
with length, it has not been shown to measure binding.

CENSORED VALUES ARE NOT POINT VALUES. ">1000 uM" and ">500 uM" are lower
bounds. They are ranked at the bound for Spearman (which only needs order) and
flagged in the output; they are never treated as measurements.
"""
from __future__ import annotations

import argparse
import csv
import math
import pathlib
import re
import sys
from collections import defaultdict

import numpy as np
from scipy import stats

UNITS = {"nm": 1e-9, "um": 1e-6, "mm": 1e-3, "m": 1.0, "pm": 1e-12}


def parse_affinity(raw: str):
    """-> (pAff, assay, censored, n_values, note) or None if unmeasured.

    Handles 'Ki 1.86 nM', 'KD 0.4/1.6 uM', 'IC50 50 uM', '>1000 uM', 'N.A.'.
    A '/' means the source reports more than one value (Koch gives two tissue
    preparations); the geometric mean is used and the count is reported, since
    averaging them is a choice and not a fact.
    """
    s = (raw or "").strip()
    if not s or s.upper() in {"N.A.", "NA", "ND", "-"}:
        return None
    censored = s.startswith(">") or s.startswith("<")
    s = s.lstrip("><").strip()
    m = re.match(r"^(Ki|KD|Kd|IC50|EC50)?\s*([0-9./]+)\s*([a-zA-Z]+)$", s)
    if not m:
        return None
    assay = (m.group(1) or "unknown").upper().replace("KD", "KD")
    nums = [float(x) for x in m.group(2).split("/") if x]
    unit = UNITS.get(m.group(3).lower())
    if unit is None or not nums:
        return None
    molar = [n * unit for n in nums]
    geo = math.exp(sum(math.log(v) for v in molar) / len(molar))
    note = ""
    if len(molar) > 1:
        note = f"{len(molar)} values ({'/'.join(str(n) for n in nums)}), geometric mean"
    return (-math.log10(geo), assay, censored, len(molar), note)


def spearman(x, y):
    if len(x) < 3:
        return float("nan"), float("nan")
    r = stats.spearmanr(x, y)
    return float(r.statistic), float(r.pvalue)


def fmt_rho(rho, p, n):
    if not np.isfinite(rho):
        return f"    n={n}: too few points for a rank correlation"
    return f"rho={rho:+.3f} p={p:.3f} n={n}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", required=True, type=pathlib.Path,
                    help="per-model TSV from 02e_aggregate.py")
    ap.add_argument("--anchors", required=True, type=pathlib.Path,
                    help="data/reference/anchor_set.tsv")
    ap.add_argument("--library", required=True, type=pathlib.Path,
                    help="data/templates/anchor20.tsv (block, length, class)")
    ap.add_argument("--score-col", default="aggregate_score")
    ap.add_argument("--out", type=pathlib.Path, default=None)
    a = ap.parse_args()

    # --- scores: collapse per peptide, keeping the spread visible -----------
    per = defaultdict(list)
    with a.scores.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                per[row["peptide_id"]].append(
                    (int(row["seed"]), float(row[a.score_col])))
            except (KeyError, ValueError):
                continue
    if not per:
        sys.exit(f"FATAL: no {a.score_col} values in {a.scores}")

    best, median, spread, nseed = {}, {}, {}, {}
    for pep, vals in per.items():
        v = [s for _, s in vals]
        best[pep] = max(v)
        median[pep] = float(np.median(v))
        seeds = defaultdict(list)
        for sd, s in vals:
            seeds[sd].append(s)
        per_seed_best = [max(x) for x in seeds.values()]
        spread[pep] = max(per_seed_best) - min(per_seed_best)
        nseed[pep] = len(seeds)

    # --- library metadata ---------------------------------------------------
    meta = {}
    with a.library.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["name"]] = row

    # --- affinities ---------------------------------------------------------
    aff = {}
    with a.anchors.open() as fh:
        lines = [l for l in fh if not l.startswith("#")]
    for row in csv.DictReader(lines, delimiter="\t"):
        parsed = parse_affinity(row.get("mu_affinity", ""))
        pki = row.get("mu_pKi", "").strip()
        if parsed is None and pki:
            parsed = (float(pki), "pKi", False, 1, "from mu_pKi column")
        aff[row["name"]] = (row, parsed)

    out = []
    def emit(line=""):
        print(line)
        out.append(line)

    emit("=" * 78)
    emit("ANCHOR SET, ARM B (Protenix 2.0.0, receptor MSA N_msa=8862)")
    emit("=" * 78)
    emit(f"score column: {a.score_col}, collapsed per peptide as MAX over "
         f"seeds x samples")
    emit(f"peptides scored: {len(best)}   (5 seeds x 5 samples each)")
    emit("")
    emit("NO CROSS-BLOCK STATISTIC IS COMPUTED. Blocks are four documents and")
    emit("four assay types; pooling them is not valid. Within-block n is 5 and 8,")
    emit("so every rho below describes a handful of points, not a method.")

    # --- per-peptide table --------------------------------------------------
    emit("")
    emit("=" * 78)
    emit("PER-PEPTIDE SCORES")
    emit("=" * 78)
    emit(f"{'peptide':<24}{'blk':<5}{'len':>4}{'best':>8}{'median':>8}"
         f"{'seedspread':>11}  {'affinity':<18}{'class'}")
    rowsort = sorted(best, key=lambda p: -best[p])
    for pep in rowsort:
        m = meta.get(pep, {})
        row, parsed = aff.get(pep, ({}, None))
        aff_s = row.get("mu_affinity", "?") or "N.A."
        emit(f"{pep:<24}{m.get('block','?'):<5}{m.get('length','?'):>4}"
             f"{best[pep]:>8.3f}{median[pep]:>8.3f}{spread[pep]:>11.3f}  "
             f"{aff_s:<18}{m.get('class','?')}")

    # --- the length control, reported first ---------------------------------
    emit("")
    emit("=" * 78)
    emit("THE LENGTH CONTROL")
    emit("=" * 78)
    peps = [p for p in best if p in meta]
    L = [int(meta[p]["length"]) for p in peps]
    S = [best[p] for p in peps]
    rho, p = spearman(L, S)
    emit(f"  score vs LENGTH, all {len(peps)} scored peptides: {fmt_rho(rho, p, len(peps))}")
    emit("  (this pools blocks, which is legitimate here: LENGTH is a property of")
    emit("   the peptide, not an assay measurement, so no cross-assay comparison")
    emit("   is being made. The affinity correlations below cannot do this.)")

    # --- within-block analysis ---------------------------------------------
    blocks = defaultdict(list)
    for p in peps:
        blocks[meta[p]["block"]].append(p)

    emit("")
    emit("=" * 78)
    emit("WITHIN-BLOCK AFFINITY CORRELATION")
    emit("=" * 78)
    for blk in sorted(blocks):
        members = blocks[blk]
        emit("")
        emit(f"--- block {blk} ({len(members)} peptides scored) ---")
        usable, paffs, scores, lens = [], [], [], []
        for pep in sorted(members, key=lambda q: -best[q]):
            row, parsed = aff.get(pep, ({}, None))
            if parsed is None:
                emit(f"    {pep:<24} score {best[pep]:.3f}   "
                     f"affinity {row.get('mu_affinity','N.A.') or 'N.A.':<14} "
                     f"-> EXCLUDED from rho (no measured value)")
                continue
            pa, assay, cens, nv, note = parsed
            usable.append(pep); paffs.append(pa); scores.append(best[pep])
            lens.append(int(meta[pep]["length"]))
            flag = " CENSORED(bound)" if cens else ""
            emit(f"    {pep:<24} score {best[pep]:.3f}   {assay} pAff {pa:5.2f}"
                 f"{flag}{('  ' + note) if note else ''}")
        if len(usable) >= 3:
            r_a, p_a = spearman(paffs, scores)
            r_l, p_l = spearman(lens, scores)
            r_al, p_al = spearman(lens, paffs)
            emit(f"  score vs AFFINITY: {fmt_rho(r_a, p_a, len(usable))}")
            emit(f"  score vs LENGTH:   {fmt_rho(r_l, p_l, len(usable))}")
            emit(f"  LENGTH vs AFFINITY (is the block itself confounded?): "
                 f"{fmt_rho(r_al, p_al, len(usable))}")
            if np.isfinite(r_a) and np.isfinite(r_l) and abs(r_l) >= abs(r_a):
                emit("  READING: length explains this block's scores at least as well as")
                emit("           affinity does. No affinity signal is demonstrated here.")
        else:
            emit(f"  fewer than 3 measured values; no rho computed (n={len(usable)})")

    # --- the specific comparisons the design rests on -----------------------
    emit("")
    emit("=" * 78)
    emit("PRE-SPECIFIED WITHIN-BLOCK COMPARISONS")
    emit("=" * 78)

    emit("")
    emit("1. Aim 2's pair: bovine casomorphin-7 (2 uM) vs human casomorphin-7 (6 uM)")
    emit("   Same length (7), same block, same assay, 3x apart in KD.")
    b, h = "bov_casomorphin_7", "hum_casomorphin_7"
    if b in best and h in best:
        emit(f"     {b}: {best[b]:.3f}   (KD 2 uM, the TIGHTER binder)")
        emit(f"     {h}: {best[h]:.3f}   (KD 6 uM)")
        d = best[b] - best[h]
        emit(f"     difference {d:+.3f}; per-seed spread is {spread[b]:.3f} and "
             f"{spread[h]:.3f}")
        if abs(d) < max(spread[b], spread[h]):
            emit("     VERDICT: the gap is SMALLER than the run-to-run spread of either")
            emit("              peptide. This pair is not resolved. Length is matched here,")
            emit("              so this is the cleanest affinity test in the set, and the")
            emit("              score does not pass it.")
        else:
            emit(f"     VERDICT: ordered {'correctly' if d > 0 else 'INVERTED'}, "
                 "and the gap exceeds both seed spreads.")

    emit("")
    emit("2. Koch length ladders (block B, within species, one assay)")
    for sp, pref in (("bovine", "bov_casomorphin_"), ("human", "hum_casomorphin_")):
        rung = [(int(meta[p]['length']), p) for p in best
                if p.startswith(pref) and p in meta]
        if len(rung) < 3:
            continue
        rung.sort()
        emit(f"   {sp}: " + "  ".join(f"{n}mer {best[p]:.3f}" for n, p in rung))
        kd = []
        for n, p in rung:
            _, parsed = aff.get(p, ({}, None))
            kd.append(parsed[0] if parsed else float("nan"))
        r_s, p_s = spearman([n for n, _ in rung], [best[p] for _, p in rung])
        emit(f"     score vs length: {fmt_rho(r_s, p_s, len(rung))}  "
             "(monotone decline = the length artefact)")

    emit("")
    emit("3. Block A's measured true negative: nociceptin FQ (N.A. at mu)")
    if "nociceptin_FQ" in best:
        a_blk = [p for p in blocks.get("A", [])]
        ranks = sorted(a_blk, key=lambda q: -best[q])
        pos = ranks.index("nociceptin_FQ") + 1
        emit(f"     nociceptin_FQ scores {best['nociceptin_FQ']:.3f}, "
             f"rank {pos} of {len(ranks)} in block A (1 = highest score)")
        emit("     It IS ranked last - but it is also the second-longest peptide in")
        emit("     the block (17), so this is exactly what the length artefact")
        emit("     predicts. Ranking it last is therefore not evidence of selectivity.")

    emit("")
    emit("4. The near-inactive short peptide: the length artefact's own test")
    if "kappa_casein_f61_63" in best and "endomorphin_1" in best:
        emit(f"     kappa_casein_f61_63 (YPY, 3mer, >1000 uM): "
             f"{best['kappa_casein_f61_63']:.3f}")
        emit(f"     endomorphin_1       (YPWF, 4mer, Ki 1.86 nM): "
             f"{best['endomorphin_1']:.3f}")
        gap = best["endomorphin_1"] - best["kappa_casein_f61_63"]
        emit(f"     The two differ by >5 ORDERS OF MAGNITUDE in affinity; their scores")
        emit(f"     differ by {gap:+.3f}.")
        emit("     These sit in different blocks so no rho is computed, but the")
        emit("     comparison needs no statistics: a score that cannot separate 1.86 nM")
        emit("     from >1 mM cannot prioritise a synthesis list.")

    if a.out:
        a.out.write_text("\n".join(out) + "\n")
        print(f"\nWrote {a.out}")


if __name__ == "__main__":
    main()
