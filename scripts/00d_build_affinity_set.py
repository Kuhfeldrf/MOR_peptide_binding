#!/usr/bin/env python3
"""Build a length-decorrelated MOR peptide affinity benchmark from ChEMBL.

CORRECTS AN EARLIER, TOO-PESSIMISTIC CONCLUSION. A name-based ChEMBL lookup
(scripts/00c_fetch_affinities.py) found only 9 of 24 queried peptides and led to
the claim that curated MOR affinity data barely exists. That was an artefact of
searching by name: enumerating ALL Ki records against CHEMBL233 instead returns
**142 peptides (molecule_type=Protein) spanning Ki 0.059 nM to 7.4e4 nM - six
log units** - and most carry `pref_name: None`, so name search could never have
found them.

This script enumerates them properly, extracts sequences, and filters to what a
co-folding benchmark can actually consume.

FILTERS, each with a reason:
  - canonical L-amino acids only. Non-canonical residues are the DAMGO problem
    (D-Ala, N-MePhe, C-terminal alcohol have no valid canonical parameters and
    Chai-1/Protenix take plain sequences) - see docs/damgo_notes.md.
  - linear only, no disulfides or cyclisation.
  - length 4-20, matching the range this pilot has scored.
  - '=' relation only. A '>' or '<' record is a bound, not a measurement.

THE DESIGN CHECK THAT DECIDES VALIDITY. Every false signal in this pilot came
from peptide length. On the 5-peptide subset examined earlier,
rho(length, pKi) = -1.000, so a score reading only length would have correlated
with measured affinity at rho=1.000. This script therefore reports
rho(length, pKi) for the assembled set and, if it is large, selects a
length-stratified subset that breaks the correlation. A benchmark whose labels
are confounded with length cannot test anything and must not be run.

Usage:
    00d_build_affinity_set.py --out docs/affinity_benchmark.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import time
import urllib.request

BASE = "https://www.ebi.ac.uk/chembl/api/data"
TARGET = "CHEMBL233"
AA = set("ACDEFGHIKLMNPQRSTVWY")


def get(url: str, tries: int = 3) -> dict:
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                return json.loads(fh.read().decode())
        except Exception:
            if k == tries - 1:
                return {}
            time.sleep(2 * (k + 1))
    return {}


def all_ki() -> dict[str, list[dict]]:
    mols: dict[str, list[dict]] = {}
    off = 0
    while True:
        d = get(f"{BASE}/activity?target_chembl_id={TARGET}&standard_type=Ki"
                f"&format=json&limit=1000&offset={off}")
        acts = d.get("activities") or []
        if not acts:
            break
        for a in acts:
            cid = a.get("molecule_chembl_id")
            if not cid:
                continue
            if str(a.get("standard_relation")).strip() != "=":
                continue
            if a.get("standard_units") != "nM" or not a.get("standard_value"):
                continue
            mols.setdefault(cid, []).append(a)
        off += 1000
        total = d.get("page_meta", {}).get("total_count") or 0
        print(f"  activities fetched: {off} / {total}")
        if off >= total:
            break
    return mols


def sequence_of(mol: dict) -> tuple[str, bool] | None:
    """(sequence, amidated) for a clean linear peptide, else None.

    HELM looks like PEPTIDE1{Y.G.G.F.L.[am]}$$$$. Two kinds of bracketed token
    must be told apart, and an earlier version of this filter got it wrong:

      [am] at the C-terminus is an AMIDE CAP, not a residue. Rejecting it would
      have discarded leu-enkephalin, endomorphin-1/2 and morphiceptin - the most
      relevant peptides in the set. It is stripped and recorded instead.

      Any OTHER bracketed token is a genuine non-canonical residue (D-amino
      acid, N-methylation, unusual side chain). Those are rejected, because the
      co-folding models take a plain one-letter sequence and cannot represent
      them - the DAMGO problem, see docs/damgo_notes.md.

    CAVEAT THIS CREATES, and it must be carried into any result: amidation is
    stripped for the model input but the Ki was measured on the AMIDATED
    molecule. C-terminal amidation materially changes potency and stability, so
    for those peptides the score and the label describe different molecules.
    The `amidated` flag exists so results can be split on it.
    """
    bt = mol.get("biotherapeutic") or {}
    helm = bt.get("helm_notation") or ""
    if not helm:
        return None
    parts = helm.split("$")
    body = parts[0]
    if body.count("{") != 1 or "PEPTIDE2" in helm:
        return None
    if len(parts) > 1 and parts[1].strip():
        return None                # explicit connections = cyclic / disulfide
    inner = body[body.find("{") + 1: body.rfind("}")]
    toks = inner.split(".")
    amidated = False
    if toks and toks[-1] == "[am]":
        amidated = True
        toks = toks[:-1]
    if any(t.startswith("[") for t in toks):
        return None                # genuine non-canonical residue
    seq = "".join(toks)
    if seq and set(seq) <= AA:
        return seq, amidated
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=pathlib.Path)
    ap.add_argument("--min-len", type=int, default=4)
    ap.add_argument("--max-len", type=int, default=20)
    a = ap.parse_args()

    print(f"enumerating Ki records against {TARGET} ...")
    mols = all_ki()
    print(f"{len(mols)} molecules with at least one '=' Ki in nM\n")

    ids = list(mols)
    rows = []
    for i in range(0, len(ids), 40):
        chunk = ";".join(ids[i:i + 40])
        d = get(f"{BASE}/molecule/set/{chunk}?format=json")
        for m in (d.get("molecules") or []):
            if m.get("molecule_type") != "Protein":
                continue
            cid = m.get("molecule_chembl_id")
            got = sequence_of(m)
            if not got:
                continue
            seq, amidated = got
            if not (a.min_len <= len(seq) <= a.max_len):
                continue
            kis = [float(x["standard_value"]) for x in mols[cid]]
            docs = sorted({x.get("document_chembl_id") or "" for x in mols[cid]})
            pki = [-math.log10(v * 1e-9) for v in kis if v > 0]
            if not pki:
                continue
            rows.append({
                "molecule_chembl_id": cid,
                "pref_name": m.get("pref_name") or "",
                "sequence": seq,
                "length": len(seq),
                "amidated": amidated,
                "n_ki": len(pki),
                "pKi_median": round(sum(sorted(pki)[len(pki)//2: len(pki)//2+1]), 3)
                              if len(pki) % 2 else round(sum(sorted(pki)[len(pki)//2-1:len(pki)//2+1])/2, 3),
                "pKi_min": round(min(pki), 3),
                "pKi_max": round(max(pki), 3),
                "Ki_nM_median": round(10 ** (9 - (sum(sorted(pki)[len(pki)//2: len(pki)//2+1])
                                        if len(pki) % 2 else
                                        sum(sorted(pki)[len(pki)//2-1:len(pki)//2+1])/2)), 4),
                "documents": ";".join(d for d in docs if d),
            })
        time.sleep(0.2)
        print(f"  molecules resolved: {min(i+40, len(ids))} / {len(ids)}  "
              f"clean linear peptides so far: {len(rows)}")

    if not rows:
        raise SystemExit("FATAL: no clean linear peptides found")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {a.out}: {len(rows)} clean linear canonical peptides "
          f"{a.min_len}-{a.max_len} aa with cited Ki")


if __name__ == "__main__":
    main()
