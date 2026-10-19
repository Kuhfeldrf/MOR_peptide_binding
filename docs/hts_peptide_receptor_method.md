# High-throughput screening of dairy peptides against enteroendocrine receptors

**Reframing, 2026-09-29.** The deliverable is not a MOR result. It is a
**validated, general method for screening a peptidome against a receptor panel**,
demonstrated on the one system where ground truth exists. Everything in
`pose_generation_comparison.md` and `validated_screening_framework.md` is
recast here as *calibration of the method*, not as the finding.

## 1. The actual problem

| axis | scale |
|---|---|
| Major bovine dairy proteins | 12 chains, ~3,700 residues |
| All contiguous 4-20mers | 61,149 |
| Realistic 3-protease in-silico digest (<=2 missed cleavages) | ~2,000-8,000 unique peptides |
| Enteroendocrine receptor panel | ~12 core, ~35 including TAS2R subtypes |
| **Complexes to evaluate** | **~24,000 to ~280,000** |

Receptor panel: OPRM1/OPRD1/OPRK1 (opioid, gut motility), **CASR** (the major gut
peptide sensor), GPR93/LPAR5 (protein hydrolysate), GPRC6A, GPR142,
TAS1R1/TAS1R3 (umami), CCKAR/CCKBR, NPY1R/NPY2R, GLP1R, GHSR, MRGPRX1/2, and the
TAS2R bitter family (~25 subtypes, expressed in gut).

This is a **matrix**, not a list, and that turns out to matter methodologically
(section 4).

## 2. Why MOR is the calibration anchor, not the target

MOR is the only receptor in the panel with all four of:

1. **Cited affinity ground truth** - 36 peptides, pKi 4.71-9.14, ChEMBL-sourced
   (`data/reference/affinity_set.tsv`).
2. **A known pharmacophore** - the Schwyzer message-address model.
3. **Held-out experimental structures** - 8F7Q (beta-endorphin, released
   2022-12-14, after the 2021-09-30 model cutoff) and 6DDF (DAMGO).
4. **A validated membrane MD system** - 150k atoms, 139.7 ns/day on one L40S.

So MOR is where you find out whether a method works *before* applying it to
receptors where you could not tell. Every negative result obtained there is
calibration data: it tells you which filters to **not** deploy across 280,000
complexes.

## 3. What transfers to other receptors, and what does not

This is the central scientific content of the calibration work.

### TRANSFERS - receptor-general, because it is a property of the scoring functions

**Every score tested is dominated by peptide length.** Seven independent methods:

| method | rho(score, length) |
|---|---|
| Protenix ipTM | -0.879 |
| Boltz-2 affinity (MW correction) | -0.857 |
| Chai-1 iLIS | -0.834 |
| ADCP docking energy | -0.718 |
| HADDOCK score (pocket-directed) | -0.647 |
| HADDOCK BSA term | **+0.961** |
| HADDOCK (undirected, prior work) | "no signal beyond peptide length" |
| FlexPepDock total_score | receptor-repacking dominated |

The mechanism is generic: buried surface area scales with peptide length
(rho=+0.961 measured), and BSA drives van der Waals and desolvation in every
empirical scoring function. **A peptidome screen spans 4-20+ residues, so any
un-normalised score will rank by size across the whole matrix, for every
receptor.** This is the single most important transferable finding and it is not
specific to opioids.

Also transferring:
- **Co-folding requires an MSA.** MSA-free Protenix gave a 16.2 A receptor with
  no pocket; with a 4,801-sequence MSA it reproduced the held-out 8F7Q
  orthosteric site to **1.03 A**. Budget one MSA per receptor (cheap, 7 s each).
- **Pose generation works; ranking does not.** Two independent routes agree -
  Protenix 1.03 A pocket, pocket-directed HADDOCK 13/14 pocket contacts.
- **Cross-seed agreement is a reliability filter.** It was the only signal that
  worked when ipTM and ipSAE failed.
- **Interaction fingerprints are pose DIAGNOSTICS, not affinity filters.** 5 of 6
  real opioids in HADDOCK poses failed to form the Asp147 salt bridge that is
  present in 100% of MD frames - while scoring -121 to -141, i.e. "good". Pocket
  occupancy is not binding mode, and no score reported the discrepancy.
- **Confidence metrics are not affinity.** Powered negative: n=36, 4.43 log
  units of Ki, length decorrelated (rho=+0.021), 88% power for rho=0.5; iLIS
  +0.196, ipTM +0.216, both ns. A single Gly3->Pro substitution worth 1600x in
  measured affinity is invisible (iLIS AUROC 0.409, i.e. backwards).

### DOES NOT TRANSFER - and this is why the expensive tiers are unavoidable

**The pharmacophore is receptor-specific.** `Tyr1 + aromatic at 3/4` reaches
AUROC 0.967 on MOR because it encodes the opioid message. There is **no
equivalent published motif** for CASR, GPR93, GPRC6A, GPR142, TAS1R1/TAS1R3, or
the TAS2R family - several of which are the most interesting dairy-peptide
sensors in the gut.

So the cheapest, best-performing filter found in the calibration work **cannot be
applied to most of the panel**. For a receptor with no known pharmacophore,
structure-based prediction is the only available approach, and the calibration
tells us exactly which parts of it to trust (pose) and which to discard
(score-as-affinity).

## 4. The matrix design solves the length artefact - and only a matrix can

A single-receptor screen cannot separate "this peptide scores high because it
binds this receptor" from "this peptide scores high everywhere because it is
long and charged." A matrix can:

**Score each peptide against all N receptors, then z-score per peptide across
receptors.** Peptide-intrinsic properties - length, charge, buried surface area,
hydrophobicity - inflate a peptide's score against *every* receptor, so they
cancel in the per-peptide z-score. What survives is **selectivity**: peptide i
scoring high for receptor j *relative to its own baseline across the panel*.

This converts the screen's biggest liability into a control, and it is testable
cheaply before committing to the full matrix (pilot: ~200 peptides x 5
receptors). Predictions to register in advance:
- lactoferricin B (+5 charge) should score high against **all** receptors and
  z-score to ~0 - it is a promiscuous artefact, and HADDOCK already ranked it
  best of 11 on MOR;
- beta-casomorphin-7 should z-score high for OPRM1 specifically.

## 5. The cascade, with measured costs

All timings measured here: Protenix 436 s/complex with MSA; GROMACS 139.7 ns/day
per L40S on a 150k-atom membrane system. Rubin assumed a conservative 2x L40S.

| tier | filter | basis | cost, full matrix |
|---|---|---|---|
| **1** | Pharmacophore, where one is known | AUROC 0.967 on MOR | free (CPU) |
| **2** | Co-folding + MSA, 5 seeds, pose only | 1.03 A pocket validated | **169,556 GPU-h** = 49 d on 144 |
| **2b** | Cross-seed agreement; interaction fingerprint vs the receptor's known anchor | reliability + pose QC | free (analysis) |
| **2c** | **Per-peptide z-score across receptors** | removes size/charge artefacts | free (analysis) |
| **3** | 20 ns membrane MD stability on survivors | binary, physical | 9,600-24,000 GPU-h = 3-7 d |
| **4** | **RBFE / pmx on the survivors** | the only thing that ranks potency | 66,000-330,000 GPU-h = 19-95 d |

**Hero run = tier 2 + tier 4: ~499,000 GPU-h = 144 days on 144 Rubin, against
4.8 YEARS on ORCA's 24 L40S.** That ratio is the resource argument, and it rests
on measured throughput rather than vendor figures.

## 6. Why each tier is on the right hardware

- **Tier 1** - CPU. Vera ARM. Trivial.
- **Tier 2** - GPU inference. This is the AI regime where Rubin gives ~10x
  rather than ~2x. Protenix/Boltz on aarch64 is **unverified** and is a named
  risk; fallback is the RTX PRO 6000 Blackwell x86 partition.
- **Tier 3/4** - GPU MD. GROMACS has NVIDIA-maintained ARM builds via NGC.
  **NVLink at 3.6 TB/s all-to-all is the specific architectural requirement**:
  lambda windows are independent, but Hamiltonian replica exchange across lambda
  needs fast all-to-all, which a conventional cluster cannot serve.
- **HADDOCK stays on ORCA.** CNS ships an `aarch64-linux.bin` so it *would* run,
  but it is CPU-only and the Huang Complex has one Vera CPU per two Rubin GPUs -
  running it there would consume the scarce resource and idle the abundant one.
  Two-site workflow, with PSU/ORCA as the CPU and staging site: this also makes
  the multi-Oregon-institution collaboration structural rather than incidental.

## 7. Honest gaps

- **Tier 4 has never been run.** No FEP/pmx in this project. Validation target
  is fixed in advance: 16 clean side-chain mutation pairs from matched-study data
  spanning 0.03-2.82 kcal/mol, including near-zero pairs that test for false
  positives. Pass = correct sign on all and MUE <= 1.5 kcal/mol.
- **The z-score design is a hypothesis, not a result.** It needs the pilot.
- **Nothing has run on aarch64.**
- **Only MOR has ground truth.** For the rest of the panel, the method can be
  validated on pose (against experimental structures where they exist) but not
  on affinity. That limit should be stated, not hidden.
- **Dairy peptidome needs real provenance** - published milk peptidomics
  datasets, not only in-silico digests, and the 36-peptide MOR affinity set is
  enkephalin/endorphin analogues, not dairy-derived.
