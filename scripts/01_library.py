#!/usr/bin/env python3
"""Stage 1 - peptide library ingestion.

Reads the curated peptide library, joins the literature-sourced affinity data
onto it by sequence, applies the filters from config.yaml, and emits a metadata
table plus a FASTA.

Nothing is silently dropped: every input row appears in the output table with an
explicit `included` flag and, when excluded, a reason.

Usage:
    01_library.py --library data/reference/known_opioid_peptides.csv \
                  --reference data/reference/reference_peptides.tsv \
                  --config config/config.yaml \
                  --outdir results/01_library
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys
from datetime import datetime, timezone

import yaml

CANONICAL = set("ACDEFGHIKLMNPQRSTVWY")


def load_reference(path: pathlib.Path) -> dict[str, dict]:
    """Affinity data keyed by sequence. Only unmodified C-termini can be
    matched to a bare library sequence; modified ones are keyed separately so
    they cannot be joined onto the wrong peptide."""
    ref: dict[str, dict] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["c_term"] != "free":
                continue          # amide / methyl ester - not the same molecule
            seq = row["sequence"].strip().upper()
            if seq == "NON_CANONICAL":
                continue
            ref.setdefault(seq, row)
    return ref


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", required=True, type=pathlib.Path)
    ap.add_argument("--reference", required=True, type=pathlib.Path)
    ap.add_argument("--config", required=True, type=pathlib.Path)
    ap.add_argument("--outdir", required=True, type=pathlib.Path)
    a = ap.parse_args()

    cfg = yaml.safe_load(a.config.read_text())
    lib_cfg = cfg["library"]
    min_len = int(lib_cfg["min_length"])
    max_len = int(lib_cfg["max_length"])
    min_score = float(lib_cfg["min_id_score"])
    dedupe = bool(lib_cfg["deduplicate"])

    a.outdir.mkdir(parents=True, exist_ok=True)
    ref = load_reference(a.reference)
    print(f"Reference affinities loaded: {len(ref)} unmodified sequences")

    rows: list[dict] = []
    seen: dict[str, str] = {}
    with a.library.open() as fh:
        for r in csv.DictReader(fh):
            seq = (r.get("sequence") or "").strip().upper()
            pid = r.get("peptide_id") or ""
            reasons: list[str] = []

            if not seq:
                reasons.append("empty_sequence")
            else:
                bad = sorted(set(seq) - CANONICAL)
                if bad:
                    reasons.append(f"non_canonical_residues:{''.join(bad)}")
                if len(seq) < min_len:
                    reasons.append(f"length_below_{min_len}")
                if len(seq) > max_len:
                    reasons.append(f"length_above_{max_len}")

            score = r.get("id_score") or r.get("ic50_um") or ""
            try:
                if score != "" and float(score) < min_score:
                    reasons.append(f"id_score_below_{min_score}")
            except ValueError:
                pass

            if dedupe and seq:
                if seq in seen:
                    reasons.append(f"duplicate_of:{seen[seq]}")
                else:
                    seen[seq] = pid

            m = ref.get(seq, {})
            # CANONICAL NAME. The screening library keys peptides as
            # real__casoxin_C while the curated reference table uses
            # casoxin_C. Carrying both conventions meant the workflow, which
            # keys on the reference name, could not find its own peptides.
            # Stage 1 is where the two tables meet, so it emits one canonical
            # name: the reference name where a row matched, otherwise the
            # library id with its real__ / scr__ prefix stripped.
            # Strip only the real__ prefix. Stripping scr__ as well would make
            # scr__casoxin_C collapse onto casoxin_C and collide with the real
            # peptide, so a lookup could silently return a scrambled control in
            # place of the peptide it is the control FOR.
            if pid.startswith("real__"):
                name = m.get("name") or pid[len("real__"):]
            else:
                name = pid
            rows.append({
                "name": name,
                "peptide_id": pid,
                "sequence": seq,
                "length": len(seq),
                "peptide_name": r.get("peptide_name", ""),
                "parent_protein": r.get("parent_protein", ""),
                "species": r.get("species", ""),
                "activity": r.get("activity", ""),
                "is_control": "yes" if pid.startswith("scr__") else "no",
                "scramble_of": r.get("scramble_of", ""),
                # joined from the literature-sourced reference set
                "ic50_um": m.get("ic50_um", ""),
                "assay": m.get("assay", ""),
                "pmid": m.get("pmid", ""),
                "reference": m.get("reference", ""),
                "affinity_status": m.get("affinity_status", "NO_REFERENCE_MATCH"),
                "benchmark_include": m.get("benchmark_include", "NO"),
                "included": "no" if reasons else "yes",
                "exclusion_reason": ";".join(reasons),
            })

    fields = list(rows[0].keys())
    tsv = a.outdir / "peptides.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["peptide_id"]):
            w.writerow(r)

    kept = [r for r in rows if r["included"] == "yes"]
    fasta = a.outdir / "peptides.fasta"
    with fasta.open("w") as fh:
        for r in sorted(kept, key=lambda x: x["peptide_id"]):
            fh.write(f">{r['peptide_id']}\n{r['sequence']}\n")

    bench = [r for r in kept if r["benchmark_include"] == "YES"]
    stats = {
        "input_rows": len(rows),
        "included": len(kept),
        "excluded": len(rows) - len(kept),
        "controls_scrambled": sum(1 for r in kept if r["is_control"] == "yes"),
        "with_affinity": sum(1 for r in kept if r["ic50_um"]),
        "benchmark_set": len(bench),
        "filters": {"min_length": min_len, "max_length": max_len,
                    "min_id_score": min_score, "deduplicate": dedupe},
    }
    (a.outdir / "library_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n")

    print(f"\nInput rows        : {stats['input_rows']}")
    print(f"Included          : {stats['included']}")
    print(f"Excluded          : {stats['excluded']}")
    print(f"Scrambled controls: {stats['controls_scrambled']}")
    print(f"With affinity     : {stats['with_affinity']}")
    print(f"Benchmark set     : {stats['benchmark_set']}")
    if stats["excluded"]:
        print("\nExclusions:")
        for r in rows:
            if r["included"] == "no":
                print(f"  {r['peptide_id']:<32} {r['exclusion_reason']}")
    print(f"\nWrote {tsv} and {fasta}")

    if not kept:
        sys.exit("FATAL: every peptide was filtered out; nothing to run.")


if __name__ == "__main__":
    main()
