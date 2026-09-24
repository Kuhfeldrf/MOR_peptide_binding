#!/usr/bin/env python3
"""Emit data/reference/reference_peptides.tsv from the curated NIH R21 Table 1.

Source: Table1_Updated_GPI_IC50.xlsx, NIH R21 Rusty Jun 2026 submission.
Every value is transcribed verbatim from that table, including its assay type
and citation. Nothing is inferred, converted, or filled in.
"""
import csv, pathlib

SRC = ("Table1_Updated_GPI_IC50.xlsx (NIH R21 Rusty Jun 2026); "
       "values transcribed verbatim 2026-09-24")

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
 ("beta_casomorphin_5_hum","YPFVE","free","agonist","human","14","GPI",
  "Koch et al. 1985; Garg Table 3 / Yoshikawa 1986","3005882",""),
 ("beta_casomorphin_7_hum","YPFVEPI","free","agonist","human","25","GPI",
  "Koch et al. 1985; Garg Table 3","3005882",""),
 ("casoxin_B_hum","YPYY","free","antagonist","human","100","GPI",
  "Chiba et al. 1989, J Dairy Res 56:363-366","2760234",
  "Sequence identical to bovine casoxin B."),
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

out = pathlib.Path("data/reference/reference_peptides.tsv")
with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["name","sequence","c_term","role","species","ic50_um","assay",
                "reference","pmid","affinity_status","ki_nM","ki_status",
                "source_table","notes"])
    for (n,s,c,r,sp,ic,a,ref,pmid,note) in ROWS:
        status = "SOURCED" if ic else "NO_VALUE_REPORTED"
        w.writerow([n,s,c,r,sp,ic,a,ref,pmid,status,"","NOT_REPORTED",SRC,note])

print(f"wrote {out} ({len(ROWS)} rows)")
n_val = sum(1 for r in ROWS if r[5])
gpi = sum(1 for r in ROWS if r[5] and r[6] == "GPI")
print(f"  with IC50: {n_val}   GPI-only: {gpi}")
