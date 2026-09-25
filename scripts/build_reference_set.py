#!/usr/bin/env python3
"""Emit data/reference/reference_peptides.tsv from the curated NIH R21 Table 1.

Source: Table1_Updated_GPI_IC50.xlsx, NIH R21 Rusty Jun 2026 submission.
Every value is transcribed verbatim from that table, including its assay type
and citation. Nothing is inferred, converted, or filled in.
"""
import csv, pathlib

SRC = ("Table1_Updated_GPI_IC50.xlsx (NIH R21 Rusty Jun 2026); "
       "values transcribed verbatim 2026-09-24")

# ---------------------------------------------------------------------------
# INTER-LABORATORY CAVEAT - read before interpreting any Stage 7 correlation.
#
# Koch 1985 Table 1 also reports GPI IC50 values for the BOVINE casomorphins,
# and they disagree sharply with the Brantl 1981 values used in this set:
#
#   peptide              Brantl 1981 (used here)   Koch 1985 Table 1
#   beta-casomorphin-4   22 uM                     3.60 uM
#   beta-casomorphin-5   6.5 uM                    0.53 uM
#   beta-casomorphin-7   57 uM                     5.14 uM
#
# Same peptide, same assay type (guinea-pig ileum), roughly an order of
# magnitude apart. The benchmark set therefore mixes laboratories, and the
# spread between labs is comparable to the spread the model is being asked to
# reproduce. Stage 7 must state this: it bounds how much any rank correlation
# over n=10 can be claimed to mean.
# ---------------------------------------------------------------------------

# name, seq, cterm, role, species, ic50_um, assay, reference, pmid, notes
ROWS = [
 ("alphas1_casein_exorphin_7","RYLGYLE","free","agonist","bovine","30","MVD",
  "Loukas et al. 1983, Biochemistry 22:4567-4573","6326810",
  "MVD assay, not GPI - not directly comparable to GPI values."),
 ("alphas1_casein_exorphin_6","RYLGYL","free","agonist","bovine","70","MVD",
  "Loukas et al. 1983, Biochemistry 22:4567-4573","6326810",
  "MVD assay, not GPI. Shorter fragment."),
 ("beta_casomorphin_4_bov","YPFP","free","agonist","bovine","22","GPI",
  "Brantl et al. 1981, Life Sci 28:1903-1909","7278645",
  "Per Brantl 1981 / Garg 2016 Table 3."),
 ("beta_casomorphin_5_bov","YPFPG","free","agonist","bovine","6.5","GPI",
  "Brantl et al. 1981, Life Sci 28:1903-1909","6265721",
  "CORRECTED from 5.2 in source table. Garg 2016 Table 3."),
 ("beta_casomorphin_7_bov","YPFPGPI","free","agonist","bovine","57","GPI",
  "Brantl et al. 1981, Life Sci 28:1903-1909","6265721",
  "CORRECTED from 14 in source table; prior value was a misattribution "
  "from human BCM-5."),
 ("neocasomorphin_6","YPVEPF","free","agonist","bovine","59","GPI",
  "Jinsmaa & Yoshikawa 1999, Peptides 20:957-962","10503774",""),
 ("casoxin_A","YPSYGLNY","free","antagonist","bovine","200","GPI",
  "Chiba et al. 1989, J Dairy Res 56:363-366","2760234",""),
 ("casoxin_B","YPYY","free","antagonist","bovine","100","GPI",
  "Chiba et al. 1989, J Dairy Res 56:363-366","2760234",""),
 ("casoxin_C","YIPIQYVLSR","free","antagonist","bovine","5","GPI",
  "Chiba et al. 1989, J Dairy Res 56:363-366","2760234",""),
 ("casoxin_D","SRYPSY","free","antagonist","bovine","","",
  "Yoshikawa et al. 1994, beta-Casomorphins and Related Peptides, VCH, 43-48","",
  "Low affinity; no IC50 reported in source."),
 ("alpha_lactorphin","YGLF","amide","agonist","bovine","50","GPI",
  "Yoshikawa et al. 1986, Agric Biol Chem 50:2419-2421","",
  "C-TERMINAL AMIDE (YGLF-NH2). No PMID in source table."),
 ("beta_lactorphin","YLLF","amide","agonist","bovine","160","GPI",
  "Yoshikawa et al. 1986, Agric Biol Chem 50:2419-2421","",
  "C-TERMINAL AMIDE (YLLF-NH2). Low potency per Teschemacher 1997. "
  "No PMID in source table."),
 ("lactoferroxin_A_bov","YLGSRY","free","antagonist","bovine","15","Radioreceptor",
  "Tani et al. 1990, Agric Biol Chem 54:1803-1810","1369293",
  "Radioreceptor assay with [3H]naloxone, not GPI."),
 ("beta_casomorphin_5_hum","YPFVE","free","agonist","human","13.5","GPI",
  "Koch, Wiedemann & Teschemacher 1985, Naunyn-Schmiedeberg's Arch Pharmacol "
  "331:351-354, Table 1","3005882",
  "VERIFIED against the paper PDF 2026-09-24. Source table said 14; Koch "
  "Table 1 reports 13.50 umol/l. Mean of 12 determinations, SD <16%."),
 ("beta_casomorphin_7_hum","YPFVEPI","free","agonist","human","29.0","GPI",
  "Koch, Wiedemann & Teschemacher 1985, Naunyn-Schmiedeberg's Arch Pharmacol "
  "331:351-354, Table 1","3005882",
  "CORRECTED. VERIFIED against the paper PDF 2026-09-24. Source table said 25; "
  "Koch Table 1 reports 29.00 umol/l. Mean of 12 determinations, SD <16%."),
 ("casoxin_B_hum","YPYY","free","antagonist","human","100","GPI",
  "Chiba et al. 1989, J Dairy Res 56:363-366","2760234",
  "Sequence identical to bovine casoxin B (YPYY). Source workbook had YPFP in its duplicate sequence column; corrected 2026-09-24."),
 ("lactoferroxin_A_hum","YLGSGY","methyl_ester","antagonist","human","","Radioreceptor",
  "Tani et al. 1990, Agric Biol Chem 54:1803-1810","1369293",
  "C-TERMINAL METHYL ESTER (YLGSGY-OCH3). No IC50 value in source. "
  "Radioreceptor assay, not GPI."),
 ("met_enkephalin","YGGFM","free","agonist","human","0.2","GPI",
  "Gacel et al. 1981, Life Sci 29:2483-2488","6276540",
  "Internal standard used by Brantl 1981."),
 ("CTOP","NON_CANONICAL","amide","antagonist","synthetic","0.003","Radioligand",
  "Gulya et al. 1986, Life Sci 38:2221-2229","2872570",
  "D-Phe-Cys-Tyr-D-Trp-Orn-Thr-Pen-Thr-NH2. IC50 = 2.80 nM. "
  "Radioligand binding, not GPI. Non-canonical: no FASTA representation."),
]

def benchmark(name, c_term, ic50, assay, role):
    """Benchmark-set membership, per the rules set 2026-09-24.

    GPI only; no C-terminally modified peptides; no CTOP; AGONISTS ONLY. No
    length cutoff - the 4-mer YPFP is retained deliberately.

    Antagonists are excluded because the receptor is 6DDF, the Gi-bound ACTIVE
    state. Antagonists preferentially bind the inactive conformation, so
    scoring them here would test the wrong receptor state and would show up in
    Stage 7 as apparent failure for a reason unrelated to the pipeline. The
    honest options were to drop them, to add an inactive-state receptor, or to
    keep them and stratify; dropping them gives the smaller but cleaner claim.
    """
    if name == "CTOP":
        return "NO", "excluded_by_instruction_non_canonical"
    if c_term != "free":
        return "NO", f"c_terminal_modification_{c_term}"
    if not ic50:
        return "NO", "no_ic50_reported"
    if assay != "GPI":
        return "NO", f"assay_not_GPI_{assay}"
    if name == "casoxin_B_hum":
        return "NO", "duplicate_sequence_of_casoxin_B"
    if role != "agonist":
        return "NO", "antagonist_excluded_active_state_receptor"
    return "YES", ""


out = pathlib.Path("data/reference/reference_peptides.tsv")
with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["name","sequence","c_term","role","species","ic50_um","assay",
                "reference","pmid","affinity_status","ki_nM","ki_status",
                "benchmark_include","exclusion_reason","source_table","notes"])
    for (n,s,c,r,sp,ic,a,ref,pmid,note) in ROWS:
        status = "SOURCED" if ic else "NO_VALUE_REPORTED"
        inc, why = benchmark(n, c, ic, a, r)
        w.writerow([n,s,c,r,sp,ic,a,ref,pmid,status,"","NOT_REPORTED",
                    inc,why,SRC,note])

print(f"wrote {out} ({len(ROWS)} rows)")
inc = [r for r in ROWS if benchmark(r[0], r[2], r[5], r[6], r[3])[0] == "YES"]
print(f"  benchmark set: n={len(inc)}")
for r in inc:
    print(f"    {r[1]:<11} {r[5]:>6} uM  {r[3]:<10} {r[0]}")
print("  excluded:")
for r in ROWS:
    ok, why = benchmark(r[0], r[2], r[5], r[6], r[3])
    if ok == "NO":
        print(f"    {r[0]:<28} {why}")
