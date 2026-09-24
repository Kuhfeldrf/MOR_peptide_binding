#!/usr/bin/env python3
"""Write the Stage 2 checkpoint report into docs/ (results/ is gitignored)."""
import csv
import pathlib
import statistics as st

root = pathlib.Path(".")
scores = list(csv.DictReader((root / "results/02_cofold/scores.tsv").open(), delimiter="\t"))
agree = list(csv.DictReader((root / "results/02_cofold/cross_seed_agreement.tsv").open(), delimiter="\t"))

by_pep: dict[str, list[dict]] = {}
for r in scores:
    by_pep.setdefault(r["peptide_id"], []).append(r)

agree_by = {a["peptide_id"]: a for a in agree}

META = {
    "real__met_enkephalin":        ("YGGFM", "0.2", "agonist", "HIGH_OVERLAP"),
    "real__beta_casomorphin_5_bov": ("YPFPG", "6.5", "agonist", "MODERATE"),
    "real__casoxin_C":             ("YIPIQYVLSR", "5", "antagonist", "MODERATE/NOVEL"),
}

lines = []
w = lines.append
w("# Stage 2 checkpoint - co-folding confidence and cross-seed agreement")
w("")
w("Three test peptides, five seeds each, Chai-1 single-sequence mode.")
w("15 inference runs producing 75 structures (5 diffusion samples per run).")
w("")
w("## Headline finding")
w("")
w("**ipTM could not distinguish the three peptides. Cross-seed agreement")
w("separated them cleanly.**")
w("")
w("Across a 32-fold affinity range and two pharmacologies, every ipTM value")
w("fell in a narrow band of 0.27-0.35. The same three peptides gave peptide")
w("RMSD across seeds of 1.6 A, 4.2 A and 7.1 A - converged, divergent and very")
w("divergent. **The confidence score was blind to a difference the geometry")
w("made obvious.**")
w("")
w("## Per-peptide results")
w("")
w("| Peptide | Seq | IC50 (uM) | Role | Overlap | ipTM (min-max) | aggregate | i-pLDDT | seed RMSD mean / max | Verdict |")
w("|---|---|---|---|---|---|---|---|---|---|")
for pid in ("real__met_enkephalin", "real__beta_casomorphin_5_bov", "real__casoxin_C"):
    rs = by_pep[pid]
    ip = [float(r["iptm"]) for r in rs]
    ag = [float(r["aggregate_score"]) for r in rs]
    il = [float(r["i_plddt_peptide"]) for r in rs]
    a = agree_by[pid]
    seq, ic50, role, ov = META[pid]
    w(f"| {pid.replace('real__','')} | {seq} | {ic50} | {role} | {ov} | "
      f"{min(ip):.3f}-{max(ip):.3f} | {st.mean(ag):.3f} | {st.mean(il):.3f} | "
      f"{a['peptide_rmsd_mean_A']} / {a['peptide_rmsd_max_A']} A | "
      f"**{a['cross_seed_agreement']}** |")
w("")
w("## Reading this honestly")
w("")
w("**1. The one peptide that converged is the one the model has seen.**")
w("Met-enkephalin (YGGFM) is an exact sequence match to 8F7Q, deposited 2023,")
w("and is flagged HIGH_OVERLAP by Stage 2b. It is the only peptide whose five")
w("seeds agree (1.64 A mean, 1.98 A max). The two peptides with no deposited")
w("mu-opioid complex diverge by 4-11 A - different sub-sites, not just")
w("different rotamers.")
w("")
w("The most economical reading is that the model reproduces a memorised pose")
w("consistently, and genuinely does not know where the other two bind. That is")
w("exactly the failure mode Stage 2b exists to expose, and it appeared on the")
w("first three peptides tested.")
w("")
w("**2. ipTM did not detect this.** Met-enkephalin's ipTM (0.286-0.304) is")
w("indistinguishable from casoxin C's (0.277-0.352), despite one being")
w("converged and the other scattered across 11 A. **A confidence score that")
w("cannot separate a memorised answer from no answer cannot be used to rank")
w("candidates.** This is direct evidence for the pipeline's stated premise:")
w("learned confidence is not affinity, and the physics stages are not optional.")
w("")
w("**3. The absolute ipTM values should not be over-interpreted yet.** All are")
w("low (0.27-0.35), but ipTM is known to be diluted when one chain is large and")
w("mostly not at the interface - a 281-residue receptor against a 5-residue")
w("peptide is close to the worst case. **ipSAE exists precisely to correct")
w("this, and is not yet computed (decision D21).** Until it is, these numbers")
w("cannot be read as evidence that the poses are bad - only that ipTM is")
w("uninformative here.")
w("")
w("**4. i-pLDDT tracks length, not affinity.** The 10-mer casoxin C scores")
w("lowest (0.459-0.497) and the two 5-mers higher (~0.52-0.54). The ordering")
w("follows peptide length, not the 32-fold affinity range. No evidence of")
w("affinity discrimination.")
w("")
w("## Runtime and budget")
w("")
w("| Quantity | Value |")
w("|---|---|")
w("| Per inference run (cached weights) | **~78 s** |")
w("| Per run including one-off weight download | 349 s |")
w("| This checkpoint: 3 peptides x 5 seeds | 19.7 min |")
w("| Structures produced | 75 (5 per run) |")
w("| **Projected full library: 40 peptides x 5 seeds** | **~4.3 h on one L40S** |")
w("")
w("Comfortably inside budget. The instructions ask for a flag if runtime")
w("exceeds the stated budget by more than ~2x; it does not.")
w("")
w("## What this does not show")
w("")
w("- Nothing here is validated. Three peptides, one receptor, one model.")
w("- Cross-seed convergence is **not** evidence of correctness. A model can")
w("  converge confidently on a wrong pose; convergence only shows consistency.")
w("- No comparison against an experimental pose has been made. Met-enkephalin's")
w("  converged pose could be checked against 8F7Q, which would say whether it is")
w("  recalling the right answer - that is worth doing and has not been done.")
w("")
pathlib.Path("docs/stage2_checkpoint.md").write_text("\n".join(lines) + "\n")
print("wrote docs/stage2_checkpoint.md")
print("\n".join(lines[:40]))
