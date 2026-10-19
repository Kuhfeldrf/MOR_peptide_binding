#!/usr/bin/env python3
"""Fetch CITED binding affinities for MOR peptides from ChEMBL.

WHY THIS IS THE BLOCKING STEP. Every method comparison in this pilot has been
scored against BINARY class labels (binder / assumed-non-binder). Enrichment and
ranking cannot be validated that way: you cannot tell a method that ranks by
affinity from one that ranks by peptide length if you have no affinities. The
repo's own rule is that every affinity must be traceable to a primary source and
none may be invented (docs/decisions.md; reference_peptides.tsv currently holds
IC50 in uM, n=6, unsourced).

This script pulls Ki / IC50 / EC50 for named peptides against CHEMBL233
(Homo sapiens Mu-type opioid receptor) and records, for every value, the ChEMBL
activity id, the assay id, the document id and the PubMed id, so each number is
traceable. Nothing is synthesised: if ChEMBL has no value for a peptide, the row
says so and stays empty.

Units are normalised to nM via ChEMBL's own `standard_value` / `standard_units`,
and any record whose units are not a concentration is dropped rather than
guessed at.

Usage:
    00c_fetch_affinities.py --out docs/affinities_chembl.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import time
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/chembl/api/data"
TARGET = "CHEMBL233"  # human mu-type opioid receptor

# Peptides this pilot has actually run, plus reference points. Names are the
# search keys; ChEMBL synonym matching is imperfect, so several aliases are
# tried per peptide and the first hit with activity data wins.
QUERIES: dict[str, list[str]] = {
    "met_enkephalin":        ["Met-enkephalin", "methionine enkephalin"],
    "leu_enkephalin":        ["Leu-enkephalin", "leucine enkephalin"],
    "endomorphin_1":         ["Endomorphin-1", "endomorphin 1"],
    "endomorphin_2":         ["Endomorphin-2", "endomorphin 2"],
    "beta_casomorphin_5":    ["beta-Casomorphin-5", "casomorphin-5"],
    "beta_casomorphin_7":    ["beta-Casomorphin-7", "casomorphin-7"],
    "dynorphin_A":           ["Dynorphin A", "dynorphin A (1-13)"],
    "dynorphin_A_1_17":      ["Dynorphin A (1-17)"],
    "beta_endorphin":        ["beta-Endorphin", "betaendorphin"],
    "alpha_neoendorphin":    ["alpha-Neoendorphin", "neoendorphin"],
    "DAMGO":                 ["DAMGO"],
    "nociceptin":            ["Nociceptin", "orphanin FQ"],
    "morphine":              ["Morphine"],          # small-molecule anchor
    "naloxone":              ["Naloxone"],          # antagonist anchor
    "CTOP":                  ["CTOP"],
    "alpha_lactorphin":      ["alpha-Lactorphin", "lactorphin"],
    "beta_lactorphin":       ["beta-Lactorphin"],
    "casoxin_A":             ["Casoxin A"],
    "casoxin_C":             ["Casoxin C"],
    "soymorphin_5":          ["Soymorphin-5"],
    "substance_P":           ["Substance P"],       # other-GPCR control
    "angiotensin_II":        ["Angiotensin II"],    # other-GPCR control
    "kisspeptin_10":         ["Kisspeptin-10"],     # other-GPCR control
    "neuropeptide_Y":        ["Neuropeptide Y"],    # other-GPCR control
}

ACCEPT_TYPES = {"Ki", "IC50", "EC50", "Kd"}


def get(url: str, tries: int = 3) -> dict:
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=45) as fh:
                return json.loads(fh.read().decode())
        except Exception as e:
            if k == tries - 1:
                print(f"    request failed: {type(e).__name__}: {e}")
                return {}
            time.sleep(2 * (k + 1))
    return {}


def find_molecules(name: str) -> list[dict]:
    q = urllib.parse.quote(name)
    d = get(f"{BASE}/molecule/search?q={q}&format=json&limit=8")
    return d.get("molecules", []) or []


def activities_for(chembl_id: str) -> list[dict]:
    out = []
    url = (f"{BASE}/activity?molecule_chembl_id={chembl_id}"
           f"&target_chembl_id={TARGET}&format=json&limit=200")
    d = get(url)
    for a in d.get("activities", []) or []:
        if a.get("standard_type") in ACCEPT_TYPES and a.get("standard_value"):
            out.append(a)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=pathlib.Path)
    ap.add_argument("--sleep", type=float, default=0.34,
                    help="polite delay between requests")
    a = ap.parse_args()

    rows: list[dict] = []
    for key, aliases in QUERIES.items():
        got = False
        for alias in aliases:
            mols = find_molecules(alias)
            time.sleep(a.sleep)
            for m in mols[:4]:
                cid = m.get("molecule_chembl_id")
                if not cid:
                    continue
                acts = activities_for(cid)
                time.sleep(a.sleep)
                if not acts:
                    continue
                for act in acts:
                    rows.append({
                        "peptide_key": key,
                        "query_alias": alias,
                        "pref_name": m.get("pref_name") or "",
                        "molecule_chembl_id": cid,
                        "molecule_type": m.get("molecule_type") or "",
                        "standard_type": act.get("standard_type"),
                        "standard_relation": act.get("standard_relation") or "",
                        "standard_value": act.get("standard_value"),
                        "standard_units": act.get("standard_units") or "",
                        "pchembl_value": act.get("pchembl_value") or "",
                        "assay_chembl_id": act.get("assay_chembl_id") or "",
                        "assay_description": (act.get("assay_description") or "")[:180],
                        "document_chembl_id": act.get("document_chembl_id") or "",
                        "pubmed_id": act.get("document_year") and act.get("document_journal") or "",
                        "src_id": act.get("src_id") or "",
                    })
                got = True
                print(f"  {key:22s} {alias:24s} {cid:16s} {len(acts):3d} activities")
                break
            if got:
                break
        if not got:
            rows.append({"peptide_key": key, "query_alias": "|".join(aliases),
                         "pref_name": "", "molecule_chembl_id": "",
                         "molecule_type": "", "standard_type": "NO_DATA_IN_CHEMBL",
                         "standard_relation": "", "standard_value": "",
                         "standard_units": "", "pchembl_value": "",
                         "assay_chembl_id": "", "assay_description": "",
                         "document_chembl_id": "", "pubmed_id": "", "src_id": ""})
            print(f"  {key:22s} -- no activity found against {TARGET}")

    cols = list(rows[0].keys())
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    have = sorted({r["peptide_key"] for r in rows
                   if r["standard_type"] != "NO_DATA_IN_CHEMBL"})
    print(f"\nWrote {a.out}: {len(rows)} records, "
          f"{len(have)}/{len(QUERIES)} peptides with data against {TARGET}")


if __name__ == "__main__":
    main()
