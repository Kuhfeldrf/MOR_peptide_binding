# Decision register

Every tool and parameter choice made in building this pipeline, with the
reasoning, the alternatives considered, the cost accepted, and the condition
that should cause it to be revisited.

**Purpose:** this is an audit trail. A reviewer should be able to check any
choice without reading the commit history, and we should be able to revisit a
choice when its assumptions change. **Reversals are recorded, not overwritten**
— a decision log that hides its own errors is not evidence of anything.

**Status key:** `FIRM` settled, revisit only on the stated trigger ·
`PROVISIONAL` working choice, expected to be re-examined ·
`OPEN` not yet decided, blocking something ·
`REVERSED` superseded, kept for the record.

---

## Index

| # | Decision | Status | Date |
|---|----------|--------|------|
| [D1](#d1) | Working root = existing scratch allocation | FIRM | 2026-09-24 |
| [D2](#d2) | MD engine: GROMACS | FIRM | 2026-09-24 |
| [D3](#d3) | GROMACS built via spack, not CMake or conda | FIRM | 2026-09-24 |
| [D4](#d4) | `cuda_arch=80,89` | FIRM | 2026-09-24 |
| [D5](#d5) | Two conda environments, not one | FIRM | 2026-09-24 |
| [D6](#d6) | Force field: CHARMM36m + CGenFF | **REVERSED** | 2026-09-24 |
| [D7](#d7) | Force field: AMBER ff19SB + LIPID21 + GAFF2 | FIRM | 2026-09-24 |
| [D8](#d8) | Ligand charges: AM1-BCC | PROVISIONAL | 2026-09-24 |
| [D9](#d9) | DAMGO extracted from 6DDF, not predicted | FIRM | 2026-09-24 |
| [D10](#d10) | DAMGO as one molecule, not residues | FIRM | 2026-09-24 |
| [D11](#d11) | Benchmark: GPI only, unmodified, no CTOP | FIRM | 2026-09-24 |
| [D12](#d12) | Human casoxin B dropped as duplicate | PROVISIONAL | 2026-09-24 |
| [D13](#d13) | IC50 and Ki kept separate, never converted | FIRM | 2026-09-24 |
| [D14](#d14) | Koch 1985 values corrected against source PDF | FIRM | 2026-09-24 |
| [D15](#d15) | Protonation: PROPKA via pdb2pqr, read from placed H | FIRM | 2026-09-24 |
| [D16](#d16) | No ECL2 loop modelling | FIRM | 2026-09-24 |
| [D17](#d17) | Chai-1 single-sequence mode, 5 seeds | PROVISIONAL | pending |
| [D18](#d18) | CGenFF academic licence | OPEN (not blocking) | — |
| [D19](#d19) | LICENSE copyright holder | DEFERRED | 2026-09-24 |
| [D20](#d20) | OSU target architecture | OPEN (narrowed) | 2026-09-24 |
| [D21](#d21) | ipSAE computed; unusable for short peptides | RESOLVED | 2026-09-24 |
| [D22](#d22) | i-pLDDT taken as peptide-chain pLDDT | PROVISIONAL | 2026-09-24 |
| [D23](#d23) | Cross-seed agreement by receptor-superposed peptide RMSD | FIRM | 2026-09-24 |
| [D24](#d24) | Stage 2 test-peptide selection | FIRM | 2026-09-24 |
| [D25](#d25) | Structure comparisons fit numbering by sequence | FIRM | 2026-09-24 |
| [D26](#d26) | Bilayer composition for intestinal muOR | FIRM | 2026-09-24 |
| [D27](#d27) | Box size raised for charged-ligand ABFE | FIRM | 2026-09-24 |
| [D28](#d28) | QC must read the tool's own verdict, and be tested against known-bad input | FIRM | 2026-09-24 |
| [D29](#d29) | Membrane must be packed periodically; contacts checked under minimum image | FIRM | 2026-09-25 |
| [D30](#d30) | Benchmark restricted to agonists; receptor is active-state only | FIRM | 2026-09-25 |
| [D31](#d31) | Bulk cation Na+, not K+ | FIRM | 2026-09-25 |
| [D32](#d32) | The workflow is the deliverable | FIRM | 2026-09-25 |

---

<a name="d1"></a>
## D1 — Working root = `/scratch/kuhfeldr-Kuhfeld_temp` · FIRM

`/scratch` is not writable at top level on ORCA; per-user project directories
are admin-provisioned as `<user>-<project>`. An allocation already existed,
empty, owned by the user. Used it rather than requesting a new `kuhfeld-temp`.
Home (92 G quota) is far too small for the 50–75 G footprint.
**Revisit if:** the allocation is reclaimed, or a quota is imposed.

<a name="d2"></a>
## D2 — MD engine: GROMACS · FIRM

Instructions permit GROMACS or OpenMM. Chosen for **portability to the OSU
NVIDIA target** ([D20](#d20)), which may be Grace/aarch64. GROMACS is C++/CMake
and builds from source on any architecture; OpenMM's CUDA support arrives via
conda-forge binaries whose aarch64 coverage could not be verified. NVIDIA also
maintains optimised ARM/GH200 GROMACS builds in NGC, matching the BioNeMo
containers already in use here.

**An earlier draft argued against GROMACS on the grounds that Stage 6 would
need hand-assembled lambda schedules. That was wrong**: GROMACS FEP is
first-class (`couple-moltype`, `init-lambda-state`, Boresch restraints via
`intermolecular_interactions`), emitting `dhdl.xvg` that `alchemlyb`/`pymbar`
consume directly. Correcting it removed OpenMM's main advantage.

**Cost accepted:** a CUDA source build on ORCA (~5 min, see [D3](#d3)).
**Revisit if:** the OSU system is confirmed x86_64 — OpenMM then wins on
`openmmtools` convenience for Stage 6.

<a name="d3"></a>
## D3 — GROMACS via spack · FIRM

Alternatives: hand-rolled CMake, or conda-forge binary. Spack records the full
dependency graph and a spec hash (`45hxd2evh2mxcjqopqa56cr3atcbwlrv`), which is
exactly the provenance artefact `docs/provenance.md` must carry — a conda
binary would record nothing about how it was built. Build took ~5 min on 16
cores, clean on the first attempt.

**Known imperfection, not corrected:** spack targeted `linux-zen4` but emitted
`-march=skylake-avx512`. GROMACS selected `AVX_512`, which Zen4 supports
natively, so the build is **correct but not tuned** — `znver4` would match, and
GROMACS on Zen4 is sometimes faster with `AVX2_256`. Left as-is because the
pilot is not a benchmark.
**Revisit before:** any scaling decision, since MD throughput drives the
GPU-hour estimate.

<a name="d4"></a>
## D4 — `cuda_arch=80,89` · FIRM

Covers both GPU types on ORCA: A30 (sm_80) and L40S (sm_89). Verified on an
L40S node — compute capability 8.9, cuFFT and NBNxM GPU setup active.
**Revisit if:** the OSU target has a different compute capability (GH200 is
sm_90, GB200 sm_100) — add it to the list rather than replacing.

<a name="d5"></a>
## D5 — Two conda environments · FIRM

`chai_lab` pulls `torch 2.6.0+cu124`, which downgrades numpy to 1.26.4 and
breaks every conda-forge package built against numpy 2.x — MDAnalysis, pymbar
and packmol-memgen all failed at import, killing Stages 3, 5 and 6. **The
conflict is structural, not a bad solve:** Chai-1 needs numpy<2, the modern
conda-forge stack needs numpy>=2, and one environment cannot satisfy both.

Resolution: `mor-pilot` (numpy>=2, all stages but 2) and `mor-chai` (numpy<2,
Stage 2 only), selected by per-rule `conda:` directives with a shared
`conda-prefix`. Also added `pdb2pqr`, which `packmol-memgen` requires but does
not pull in transitively.

**Cost accepted:** two environment files to pin rather than one.
**Revisit if:** chai_lab gains numpy>=2 support — then re-merge.

<a name="d6"></a>
## D6 — CHARMM36m + CGenFF · **REVERSED**

Chosen on literature precedent: published DAMGO–μOR MD uses CHARMM36m with
CGenFF in GROMACS. **Reversed the same day.** The supporting argument was
incomplete — it ruled out AMBER because no N-methyl amino-acid *residue*
library exists for ff19SB, then proposed treating DAMGO as a single whole
molecule under CGenFF, without applying that same reframing back to AMBER.
GAFF2 parameterises arbitrary organic molecules including N-methylated amides
and alcohols, so the cited blocker never applied to a whole-molecule treatment.

Superseded by [D7](#d7). Kept because the reversal is itself evidence about how
the decision was reached.

<a name="d7"></a>
## D7 — AMBER ff19SB + LIPID21 + GAFF2 · FIRM

| Component | Choice |
|-----------|--------|
| Protein | ff19SB |
| Membrane | LIPID21 (POPC + cholesterol) |
| Water | OPC (ff19SB's matched model) |
| DAMGO | GAFF2, single molecule |
| Builder | `packmol-memgen`, default Amber output |

**Deciding principle, set by the user: "do what we can actually run and
reproduce."** CGenFF generation is licence-gated and interactive; `antechamber`
is already installed and verified. For a repository whose purpose is
demonstrating reproducibility, a licence gate mid-pipeline is a worse defect
than less direct precedent.

Hard constraint: **ligand FF must match system FF** — GAFF2 with CHARMM36m is
invalid (different Lennard-Jones combination rules). Choosing the ligand force
field chooses the whole stack.

**Costs accepted:**
1. Amber→GROMACS topology conversion returns as a silent-error risk. Stage 3
   must verify atom count and total charge across it and fail loudly.
2. Less direct precedent for this exact receptor.

**Revisit if:** a cross-check against the published CHARMM36m work is wanted —
`packmol-memgen --charmm` exists, so that path stays open.

<a name="d8"></a>
## D8 — Ligand charges: AM1-BCC · PROVISIONAL

Standard and well validated for ligand work, far cheaper than RESP, and runs
via `antechamber`/`sqm` with no QM package.
**Revisit if:** Stage 6 convergence is poor or DAMGO's ΔG is badly off the
published affinity. **Escalation is RESP at HF/6-31G(d,p)** — the level Amber's
own charges were derived at, and the level used in published Amber-based
μ-opioid ligand work. Check charges before blaming sampling.

<a name="d9"></a>
## D9 — DAMGO extracted from 6DDF, not predicted · FIRM

DAMGO has no canonical FASTA form, so it cannot be handed to Chai-1 as a
sequence. It is **experimentally resolved in 6DDF as entity 5** (chain D,
`TYR-DAL-GLY-MEA-ETA`, 37 heavy atoms), which is already the Stage 0 receptor
source. Verified directly against the mmCIF.

Better than predicting it: ABFE starts from the experimental pose rather than
compounding pose error into the free energy, and receptor and ligand share a
coordinate frame so no docking or alignment step is needed.

**Cost, which must be stated in the README:** DAMGO's ABFE validates the
free-energy machinery only — decoupling, restraints, MBAR, convergence — and
**does not test pose prediction**.

<a name="d10"></a>
## D10 — DAMGO as one GAFF2 molecule · FIRM

Alternative was ff19SB standard residues plus two custom ones (`MEA`, `ETA`).
Rejected: that creates a **mixed-force-field junction mid-peptide** where
charges must be reconciled across the boundary, a known source of quiet error,
and it needs an N-methyl amino-acid residue library that does not exist for
ff19SB. At 513 Da DAMGO is small-molecule sized, so whole-molecule treatment is
appropriate.
**Revisit if:** conformational sampling of the peptide backbone looks wrong —
GAFF2 torsions are not trained on peptide backbones the way ff19SB is.

<a name="d11"></a>
## D11 — Benchmark: GPI only, unmodified, no CTOP · FIRM

Set by the user. **No length cutoff** — the 4-mers YPFP and YPYY are retained
deliberately, which is what lifts the antagonist count from 2 to 3 (the source
table's own ≥5 AA rule had made antagonists "insufficient for correlation").

Excluded: αs1-casein exorphins (MVD assay), lactoferroxin A (radioreceptor),
α-/β-lactorphin (C-terminal amide), human lactoferroxin A (methyl ester),
casoxin D (no IC50), CTOP (non-canonical, per instruction).

Result: **n=10**, 7 agonists / 3 antagonists, 0.2–200 µM. Membership is
precomputed into a `benchmark_include` column with a reason recorded for every
exclusion, so the filter is auditable rather than implicit in Stage 7.

<a name="d12"></a>
## D12 — Human casoxin B dropped as duplicate · PROVISIONAL

Not specified by the user; my call. It carries the **same sequence (YPYY), same
100 µM, and same Chiba 1989 citation** as bovine casoxin B — including both
would count one measurement twice in the Spearman. One flag in the TSV reverses
it.
**Revisit if:** the user wants species-level replication retained.

<a name="d13"></a>
## D13 — IC50 and Ki kept in separate columns, never converted · FIRM

The instructions ask for Ki; the available data is IC50 in µM. **No
Cheng–Prusoff conversion is applied** — it requires the radioligand Kd and
concentration, which the sources do not report. `ki_nM` stays empty and flagged
`NOT_REPORTED` throughout.

<a name="d14"></a>
## D14 — Koch 1985 values corrected against the source PDF · FIRM

Verified against the paper rather than the summary table. Koch, Wiedemann &
Teschemacher 1985, *Naunyn-Schmiedeberg's Arch Pharmacol* 331:351–354, Table 1
(GPI, µmol/l, means of 12 determinations, SD <16%):

| Peptide | Summary table | Koch Table 1 |
|---|---|---|
| Human β-CM-5 | 14 | **13.50** (rounding) |
| Human β-CM-7 | 25 | **29.00** (wrong by 16%) |

**Inter-laboratory caveat, recorded in the reference set:** Koch also reports
bovine β-CM-4/-5/-7 at 3.60 / 0.53 / 5.14 µM against the Brantl 1981 values of
22 / 6.5 / 57 used here. Same peptides, same assay type, ~10× apart. **The
benchmark mixes laboratories, and the between-lab spread is comparable to the
signal Stage 7 is meant to detect.** This bounds what any n=10 correlation can
claim and belongs in the figure caption.
**Revisit:** consider standardising on a single laboratory's values.

<a name="d15"></a>
## D15 — Protonation: PROPKA via pdb2pqr, state read from placed hydrogens · FIRM

`pdb2pqr30 --ff=AMBER --with-ph 7.4 --titration-state-method propka`.

**State is determined from the hydrogens actually placed, not from residue
names.** `pdb2pqr --pdb-output` retains generic `HIS`/`ASP` names and cannot
distinguish HID/HIE/HIP or ASP/ASH — an early version read names and returned
"unknown" for H297. Test: `HD1`/`HE2` on the histidine ring nitrogens, `HD2` on
aspartate OD2. (`HD2`/`HE1` on His sit on carbons and are not diagnostic.)

Result at pH 7.4: **D147 deprotonated (−1)** — chemically expected for the
conserved TM3 aspartate that salt-bridges the ligand amine — and **H297 = HID**,
neutral, δ-protonated.
**Revisit if:** a different pH is simulated, or if D147's protonation is
suspected of affecting ligand binding thermodynamics.

<a name="d16"></a>
## D16 — No ECL2 loop modelling · FIRM

The instructions call for modelling missing residues "particularly ECL2".
**6DDF has no chain breaks** — residues 65–345 are continuous — so ECL2 is
fully resolved and there is nothing to model. Running a modelling step anyway
would add a gratuitous source of non-determinism.

Related: **6DDF contains no nanobody and no scFv16**; only the Gi heterotrimer
(A/B/C). The instruction to strip them is a no-op here. Both facts are logged
by Stage 0 rather than silently skipped.
**Revisit if:** the 8F7Q alternates are processed — they may have gaps.

<a name="d17"></a>
## D17 — Chai-1: single-sequence mode, 5 seeds · PROVISIONAL

From the instructions: single-sequence (no MSA) avoids multi-terabyte database
hosting; 5 seeds per peptide with the full ensemble retained; cross-seed
agreement reported as an explicit diagnostic. Not yet executed.
**Revisit if:** cross-seed agreement is so poor that ranking is meaningless, in
which case the no-MSA choice is the first suspect.

<a name="d18"></a>
## D18 — CGenFF academic licence · OPEN, not blocking

Superseded as a dependency by [D7](#d7), so nothing is blocked. Retained as an
option: the web app at `app.cgenff.com/signup` is free and instant for
academics; the scriptable binary needs an institutional signatory and takes
weeks. **Only becomes relevant** if a CHARMM cross-check is wanted, or if
scaling makes a per-molecule web step untenable.

<a name="d19"></a>
## D19 — LICENSE copyright holder · DEFERRED by user

`LICENSE` retains its placeholder. The user chose on 2026-09-24 to leave it as
is for now. Blocks nothing until the repository is distributed; a one-line
change when decided.

Worth re-checking before distribution: for NIH-funded work, copyright often
sits with the institution under its IP policy rather than with the individual.

<a name="d20"></a>
## D20 — OSU target architecture · OPEN (narrowed 2026-09-24)

**The question is one fact: are the target machine CPUs ARM (Grace/Vera) or
x86?**

### What was established by public sources

- The facility is the **Jen-Hsun and Lori Huang Collaborative Innovation
  Complex** at OSU, a $200M centre funded partly by a $50M gift from NVIDIA's
  founder and his spouse, both OSU graduates.
- The machine is an **NVIDIA DGX SuperPOD plus OVX SuperPOD**, reported as
  roughly **60 DGX and OVX systems**.
- At least one report describes **Vera Rubin GPUs** and **NVL144** rack-scale
  configuration.

### Why that still does not settle it

**NVIDIA official material says only "next-generation CPUs" and names no
architecture.** The two possibilities point opposite ways:

- **NVL144 / Vera Rubin** is a rack-scale coherent CPU+GPU design, and the
  **Vera CPU is ARM** (successor to Grace). If this is the configuration, the
  target is **aarch64**.
- **DGX SuperPOD** has historically also shipped with **x86** host CPUs.

Since reports mention both DGX/OVX systems *and* NVL144, the configuration is
genuinely ambiguous from outside. **This should not be resolved by inference.**

### What would resolve it

Internal OSU sources, which the user has access to and public reporting does
not: research computing / the Huang Complex technical documentation, or the
**Huang Complex Supercomputing Seed Fund** programme pages, which typically
document the system available to applicants.

### Consequence, and why nothing is blocked

If ARM, every Python wheel in the stack needs an availability check before
porting, and `arch_notes.md` becomes a real porting assessment. If x86, the
port is close to free and [D2](#d2) (GROMACS over OpenMM) could be reopened.

**The lean toward Vera/ARM, if it holds, vindicates the GROMACS choice** —
GROMACS builds from source on any architecture, whereas OpenMM's CUDA support
depends on conda-forge aarch64 coverage that could not be verified.

Until confirmed, `docs/arch_notes.md` stays explicitly labelled **x86_64
evidence only**, and the README scaling section records the porting assessment
as pending. That is a legitimate state for the document; implying an assessment
of a machine whose architecture is unknown would not be.

Sources: OSU Newsroom, "$50 Million Gift by NVIDIA Founder and Spouse Helps
Launch Oregon State University Research Center"; NVIDIA blog, "AI Supercomputer
to Power $200 Million Oregon State University Innovation Complex".


<a name="d21"></a>
## D21 - ipSAE: computed, and found unusable at this peptide length · RESOLVED

The instructions require ipSAE. Chai-1 does not produce it, so the reference
implementation was **vendored** (`vendor/ipsae.py`, Dunbrack lab, MIT licence)
and only a format adapter written (`scripts/02c_ipsae.py`), targeting the
script's Boltz path, which takes a PAE `.npz` plus a `.cif` - the pair Chai-1
already emits. A home-rolled scoring function that was subtly wrong would look
entirely plausible in a results table.

### The expectation

Stage 2 found ipTM flat at 0.27-0.35 across a 32-fold affinity range, unable to
separate a memorised and correct pose from two the model could not place.
ipSAE exists to fix precisely that defect: ipTM scores whole chains, so pairs
far from the interface dilute it, and 281 residues against 5 is close to worst
case. The expectation was that ipSAE would restore the discrimination.

### The result

| Peptide | Interface residues | ipSAE (mean) |
|---------|-------------------|--------------|
| Met-enkephalin (converged, correct) | 5 | **0.046** |
| beta-casomorphin-5 (divergent) | 5 | **0.046** |
| Casoxin C (divergent) | 7-10 | **0.031** |

**ipSAE did not separate them either.** The memorised, experimentally validated
pose scores identically to one the model placed 4-11 A apart across seeds.

### Why - a structural limit, not a bad run

The input is sound: 286 tokens, interchain PAE mean 7.33 A with a minimum of
3.76 A, which is a reasonable interface.

The cause is in `calc_d0`: for an interface of **27 residues or fewer, d0 is
clamped to 1.0 A**. Each pair contribution is then roughly 1/(1+(PAE/d0)^2), so
with PAE around 7 A every pair contributes about 0.02 and the total collapses
to the 0.03-0.05 band observed. **A 4-10 residue peptide interface cannot score
well on ipSAE regardless of how correct the pose is.**

### Conclusion carried into the README

**Both learned-confidence metrics are unfit for ranking peptides of this length
against a GPCR, for opposite structural reasons:**

- **ipTM** is *diluted* by the chain-size mismatch - too many irrelevant pairs.
- **ipSAE** is *compressed* by its d0 floor at small interface size - too few
  relevant ones.

The residual ordering ipSAE does show (casoxin C below the 5-mers) tracks
interface residue count, not pose quality.

**Cross-seed agreement remains the only Stage 2 diagnostic that separated a
correct pose from a failed one**, and it did so decisively (1.6 A against
4-11 A). That is the column Stage 7 should stratify on, and the confidence
scores should be reported as what they are rather than used for ranking.

ipSAE is still computed and reported, because the instructions ask for it and
because its failure here is informative. It is not used to rank.

<a name="d22"></a>
## D22 — i-pLDDT taken as mean pLDDT over the peptide chain · PROVISIONAL

Chai-1 returns a per-token pLDDT. "Interface pLDDT" has no single canonical
definition. Taken here as **the mean over the peptide chain tokens**, because
averaging over the whole complex would be swamped by a 281-residue well-folded
receptor and would not report on the interface at all. The smoke test makes the
scale of that concern concrete: receptor pTM 0.868 against peptide pTM 0.137.

**A stricter definition** — mean pLDDT over residues within a distance cutoff
of the partner chain — is defensible and arguably better. Not used yet because
it introduces a cutoff parameter needing its own justification.
**Revisit if:** i-pLDDT is used for ranking rather than as a diagnostic.

<a name="d23"></a>
## D23 — Cross-seed agreement: receptor-superposed peptide RMSD · FIRM

Seeds are compared by superposing **receptor CA atoms only**, then measuring
the peptide displacement in that frame. Superposing on the whole complex would
be dominated by the receptor and would report agreement even when the peptide
lands in a different sub-site. The question is whether the peptide lands in the
same place on the same receptor, and this measures that directly.

Reported as bands, not pass/fail: `CONVERGED` under 2 A, `PARTIAL` 2-5 A,
`DIVERGENT` over 5 A, plus centroid shift. Roughly 2 A is the scale of a
well-converged pose and beyond 5 A generally means a different sub-site.

**Why it matters:** disagreement across seeds flags an untrustworthy prediction
independently of the confidence the model states. High ipTM with poor
cross-seed agreement is a stronger warning than either number alone.

<a name="d24"></a>
## D24 — Stage 2 test-peptide selection · FIRM

Three peptides chosen to span the axes that could confound the checkpoint,
rather than to give a flattering result:

| Peptide | Seq | Len | IC50 | Role | Overlap |
|---------|-----|-----|------|------|---------|
| Met-enkephalin | YGGFM | 5 | 0.2 uM | agonist | **HIGH_OVERLAP** (exact match, 8F7Q) |
| beta-casomorphin-5 (bov) | YPFPG | 5 | 6.5 uM | agonist | MODERATE |
| Casoxin C | YIPIQYVLSR | 10 | 5 uM | **antagonist** | NOVEL/MODERATE |

Spans a 32x affinity range, both pharmacologies, 5-mer against 10-mer, and
memorised against not. Including the HIGH_OVERLAP case is deliberate: if
Met-enkephalin predicts far better than the others, that is a memorisation
signal, and it is better to see it at this checkpoint than in Stage 7.


<a name="d25"></a>
## D25 — Residue numbering is fitted by sequence, never assumed · FIRM

**Recorded because it produced wrong published-looking numbers before it was
caught.**

The first version of `validate_vs_experiment.py` mapped the Chai-1 prediction
onto 8F7Q by assuming the prediction used the 6DDF construct numbering, which
starts at residue 65. **8F7Q is numbered +2 relative to that construct.** The
assumed mapping gave **8.9% sequence identity** at matched positions; the
correct offset of 67 gives **98.2%**.

Every RMSD from the bad mapping was meaningless, and none of them looked
obviously wrong:

| Quantity | Wrong mapping (offset 65) | Correct (offset 67) |
|---|---|---|
| Receptor RMSD, global | 6.09 A | **2.80 A** |
| Pocket RMSD | 5.42 A | **0.80 A** |
| Peptide RMSD, pocket frame | 1.88 A | **2.32 A** |

The 6 A receptor deviation prompted a plausible but entirely invented
explanation - that the receptor was in a different activation state, with TM6
swung out. **That interpretation was an artefact of the misalignment and is
retracted.** The receptor agrees at 2.8 A, which is unremarkable.

**Rule going forward:** any script that maps residues between two structures
fits the offset by maximising sequence identity and **asserts that identity
exceeds 90%**, failing loudly otherwise. Numbering is never inferred from
provenance.

**Why this matters beyond one script:** a misaligned comparison produces
numbers in an entirely plausible range. There was no error, no exception, and
the result invited a mechanistic story. Only an independent check - amino-acid
identity at matched positions - exposed it.


<a name="d26"></a>
## D26 — Bilayer composition for intestinal muOR · FIRM

The project targets muOR in enteroendocrine cells of the intestine, so the
bilayer must suit that environment rather than a generic membrane. Intestinal
epithelial membranes are **strongly asymmetric between apical and basolateral
faces**, so "intestinal" alone does not specify a composition.

### Evidence

**Apical (brush border) membrane** is extreme: cholesterol : phospholipid :
glycolipid at roughly **1 : 1 : 1** by mole, i.e. about 50 mol% cholesterol,
with high glycosphingolipid content and correspondingly low fluidity.

**Basolateral membrane** is far more conventional: roughly **1 : 2.5 : 0.3**,
giving cholesterol at **~28-29 mol%** of cholesterol-plus-phospholipid.

### Which face, and which cells

Published localisation places intestinal muOR **mainly on enteric neurons** of
the myenteric and submucosal plexus, with epithelial expression reported in
goblet cells and associated with **basolateral** function - regulation of
chloride secretion, alongside the basolateral VIP and muscarinic receptors.
Beta-casomorphins have to cross the epithelium (shown in Caco-2 transport
studies) to reach these sites. Enteroendocrine cells apically **secrete**
endogenous opioids such as beta-endorphin and Met-enkephalin, which is a
separate matter from where the receptor sits.

### The decisive argument comes from our own benchmark data

**Every IC50 in the Stage 7 benchmark is from the GPI assay - guinea-pig ileum
longitudinal muscle / myenteric plexus** (decision D11).

**To be precise about what "neuronal" means here:** the myenteric plexus is
part of the **enteric nervous system**, roughly 500 million neurons embedded in
the wall of the intestine itself. These are **intestinal tissue, not brain
cells**. The GPI preparation is gut throughout. The distinction being drawn is
therefore not gut-versus-brain, but *which cell type within the gut wall* -
enteric neurons in the muscle layer, rather than the epithelial cells lining
the lumen.

So the benchmark affinities report on muOR in an **enteric neuron** plasma
membrane, not an enterocyte brush border.

Simulating an apical brush-border composition would mismatch both the
receptor's likely location and the data the pipeline is benchmarked against.
**A basolateral/neuronal-like membrane is the consistent choice.**

### Decision

**POPC : cholesterol at 7:3 (30 mol% cholesterol).**

- The build instructions specify "POPC with cholesterol". Moving to a
  multi-component bilayer (PE, PS, sphingomyelin, glycosphingolipids) would be
  a substitution of a specified component and is **not made silently**; it is
  recorded as a limitation below instead.
- 30 mol% is chosen to sit at the basolateral figure (~28-29 mol%) rather than
  the generic 25% used in the first build, which was not derived from anything.
- Cholesterol is not cosmetic for this receptor: muOR structures including
  6DDF and 8F7Q have cholesterol resolved, and GPCR function is known to depend
  on it.

### Limitations, to carry into the README

- **No leaflet asymmetry.** Real plasma membranes keep PS and PE inner-facing.
  Symmetric here.
- **No PE, PS, sphingomyelin or glycosphingolipid.** A basolateral membrane
  contains all of these. POPC stands in for the whole phospholipid fraction.
- **No glycocalyx**, which is a defining feature of the intestinal apical face
  but not of the basolateral one modelled here.
- These simplifications follow the build instructions. They mean the membrane
  is a *reasonable basolateral/neuronal mimic*, **not** a model of the
  enterocyte apical membrane, and results should not be read as the latter.

### A gap in the justification, stated plainly

The ~28-29 mol% figure is measured for the **intestinal epithelial basolateral
membrane**. The benchmark data, however, comes from **enteric neurons** - a
different cell type. **No enteric-neuron-specific membrane lipidomics was
found**, and none is claimed here.

30 mol% is therefore supported directly for the basolateral epithelial case and
is *plausible but not measured* for enteric neurons, where generic neuronal
plasma membranes fall in a similar range. It is a defensible choice for both,
not a measured value for the tissue the affinities actually come from.

**How to close this:** find enteric-neuron or peripheral-neuron plasma membrane
lipidomics and check 30 mol% against it. Until then the figure stands on the
epithelial measurement plus the general neuronal range, and that is what the
README should say.

<a name="d27"></a>
## D27 — Box size raised for charged-ligand ABFE · FIRM

The first build used `--dist 15 --dist_wat 17.5`, giving an 80 x 80 A box,
210 lipids and 69,656 atoms - **below the 100,000-150,000 the instructions
anticipate**. That was an accident of parameters rather than a decision.

Raised to `--dist 18 --dist_wat 22.5`. Two reasons, both specific:

1. **DAMGO carries a net charge of +1.** Absolute binding free energies for
   charged ligands have well-known finite-size artefacts under Ewald
   electrostatics, and the correction scales with box dimension. A larger box
   reduces the artefact at source rather than relying on a post-hoc correction.
   Stage 6 is the pilot's validation target, so this is the wrong place to
   economise.
2. **The complex spans 74 A in z against a 31 A bilayer**, so roughly 43 A of
   receptor projects beyond the membrane into solvent. A 17.5 A water layer
   leaves that region close to its periodic image.

**Cost:** more atoms means proportionally more MD time per nanosecond. Accepted;
the pilot's MD is short and Stage 6 accuracy matters more than Stage 4 speed.


<a name="d28"></a>
## D28 - QC reads the tool's own verdict, and is tested against known-bad input · FIRM

**Recorded because a failure was announced by the tool, in writing, and went
unread through an entire build, a topology conversion and two minimisation
attempts.**

### What happened

packmol did not converge. Its log ended:

```
packing problem with the desired distance tolerance.
...
contains the best solution found.
STOP 173
```

It had written its best attempt rather than a converged pack, leaving 551
inter-molecular pairs under 1.2 A, the worst at **0.090 A**, across 3,678
molecules - mostly water-water and lipid-lipid. GROMACS then reported
`Maximum force = inf` and steepest descent quit after 16 steps at
5.2e17 kJ/mol.

### Why the QC did not catch it

Two independent blind spots:

1. **It checked lipid-PROTEIN contacts only.** The actual overlaps were
   lipid-lipid and water-water, which it never looked at. The check was aimed
   at the failure mode I had imagined rather than the one packmol's tolerance
   actually governs.
2. **It never read packmol's status line.** The most reliable evidence
   available - the tool's own verdict on its own work - was ignored in favour
   of geometry I computed myself.

### Rules adopted

- **Read the tool's verdict first.** Before computing anything, parse whatever
  the upstream tool says about whether it succeeded. It knows more about its
  own convergence than any downstream geometric proxy.
- **Check the quantity the upstream tool is actually controlling.** packmol
  enforces a minimum distance between *molecules*; so the QC counts contacts
  between molecules, not between the two species I happened to be thinking of.
- **Test the QC against known-bad input.** QC v2 was run against the failed
  pack and required to reject it before being trusted on a new one. A check
  only ever validated against data believed to be good has never been shown to
  detect anything.

### A related self-inflicted error, recorded

The first attempt to rescue the bad pack moved individual atoms apart with no
restoring force on their bonded partners. Accumulated independent pushes
stretched C-H bonds to 3.05 A (r0 1.09) and tore water H1-H2 to 2.99 A
(r0 1.371). The structure was more broken after the repair than before it, and
minimisation still failed on the same atom. **Repairing coordinates atom-by-atom
without respecting molecular connectivity is not a valid remedy**; either move
whole molecules rigidly, or re-pack. Here the correct answer was to re-pack,
because 3,678 molecules were involved - a system-wide convergence failure, not
a few local defects.

### Earlier self-caught QC defects, kept for the same reason

- QC v1 **passed vacuously**: it searched for `POPC`/`CHL1` while Lipid21 writes
  the modular `PC` + `PA` + `OL` and `CHL`. It found zero lipids, so every
  geometric test silently had nothing to test, and it reported PASSED. A check
  that finds nothing to check now fails.
- Its first enclosure test was a raw neighbour count, which flagged 1159 atoms
  sitting in ordinary annular grooves of a 7-TM bundle. Replaced with a
  convex-hull enclosure test, which asks whether protein SURROUNDS the atom
  rather than whether protein is merely near it.


<a name="d29"></a>
## D29 - The membrane must be packed periodically, and contacts checked under minimum image · FIRM

**The single most consequential error in Stage 3, and it was invisible to every
check that had been written.**

### What was wrong

The system was never built periodically. packmol packs molecules into a
*region*; it knows nothing about periodic boundaries. Lipids at the edge of the
patch extend past it, so under PBC they overlap the lipids on the opposite
face.

Measured:

| Quantity | x | y | z |
|----------|---|---|---|
| Lipid extent | 91.1 A | 91.5 A | 48.6 A |
| Box | 86.06 A | 86.06 A | 119.0 A |
| **Overhang** | **+5.0** | **+5.4** | +1.7 |

14.4% of atoms lay outside a single box cell. Under minimum-image convention
the system carried **167 inter-molecular pairs under 0.5 A**, worst 0.118 A.

The protein was never implicated: it spans only 42.8 x 46.7 A.

### Why it took so long to find

**Every contact search used raw Cartesian distances.** All of them were
structurally blind to periodic images. The consequence was a system that looked
clean by every available measure - worst Cartesian contact 1.767 A, zero bonds
over 3 A, worst angle off by 17 degrees, rigidity verified to 1e-13 A - while
both GROMACS and sander reported astronomical forces on an atom whose nearest
Cartesian neighbour was 2.5 A away.

That contradiction was the clue, and it should have been treated as one much
sooner. An atom with no neighbour inside 2.5 A cannot carry an infinite force
under any correct force field, so either the force field was wrong or **the
distances being measured were not the distances the engine was using**. The
second is what was true.

The same data, measured correctly:

| Cutoff | Cartesian | Minimum image |
|--------|-----------|---------------|
| < 0.5 A | 0 | **90** |
| < 1.2 A | 55 | **1,315** |

### What ruled out the wrong hypotheses

Running the minimisation in **sander** was decisive. Both engines failed on the
same atom, which eliminated the Amber-to-GROMACS conversion - the leading
suspect at the time, since OPC is a 4-site model with a virtual site. Without
that test the search would have continued down the conversion path.

### Rules adopted

1. **Pack periodically.** `packmol-memgen --pbc` makes packmol respect the
   boundary and adapt its constraints, so the patch tiles instead of
   overhanging.
2. **Check contacts under minimum image, always.** The QC now reads the box
   from the packing regions, compares it against the coordinate extent, and
   runs a minimum-image contact census. A Cartesian-only check on a periodic
   system answers a question nobody asked.
3. **When a measurement contradicts an engine, suspect the measurement.** Two
   independent engines agreeing on an impossible force outweighed a geometric
   check that said everything was fine.

### Cost

Four packs (26 min, 5h39m, 1h09m, and the periodic rebuild) and five
minimisation attempts. Two of the three diagnoses along the way were real and
kept: the bilayer was over-packed (`--apl_offset`, an 11-fold reduction in
pathological contacts) and the residual overlaps needed rigid-body separation
rather than per-atom nudging. But the governing defect was the missing `--pbc`,
and it would have been found on the first attempt by checking the coordinate
extent against the box - a two-line comparison.


<a name="d30"></a>
## D30 - Benchmark restricted to agonists · FIRM

The receptor is **6DDF, the Gi-bound ACTIVE state**. Antagonists preferentially
bind the **inactive** conformation, so scoring them against this receptor tests
the wrong state. Three of the ten benchmark peptides were antagonists - casoxins
A, B and C - and they would have appeared in Stage 7 as failures for a reason
that has nothing to do with whether the pipeline works.

**Decision: exclude antagonists. The benchmark is n = 7 agonists**, spanning
0.2 to 59 uM, roughly a 300-fold range.

### Options considered

| Option | Effect |
|--------|--------|
| Restrict to agonists | Smaller n, cleanest claim. **Chosen.** |
| Add an inactive-state receptor (e.g. 4DKL) | Scientifically strongest; roughly doubles Stage 3-4 cost |
| Keep all ten, stratify by pharmacology | Makes the mismatch visible but still reports a number computed in the wrong state |

The second is the better science and the better story for a scale-up argument -
"the pipeline handles agonists and antagonists against their respective
receptor states" - and it is recorded here as the obvious extension rather than
as something ruled out.

### What this costs, stated in the README

n drops from 10 to 7. A rank correlation over seven points is weakly
determined, and Stage 7 must report n on the figure rather than only in the
caption. Combined with the inter-laboratory spread already recorded in D14,
this bounds what the benchmark can establish: it is a demonstration that the
comparison runs end to end, not a validation of the physics.

<a name="d31"></a>
## D31 - Bulk cation is Na+, not K+ · FIRM

packmol-memgen defaults to K+. The receptor's extracellular and interstitial
faces sit in **Na+-dominated fluid** - roughly 140 mM Na+ against 4 mM K+ -
while K+ dominance is the **intracellular** condition.

This is not cosmetic for this receptor family. Opioid receptors carry a
**conserved allosteric Na+ site at D2.50**, which in muOR is **Asp114** - the
residue PROPKA protonated to ASH in Stage 0. That protonation is *correct* for
the active state, where the sodium pocket collapses, so the model is internally
consistent. The point of this decision is that it should be a **stated choice
about a known allosteric site**, not an accident of a tool default.

**Revisit if** an inactive-state receptor is added ([D30](#d30)): there D2.50
is deprotonated and coordinates Na+, so the protonation must change with the
conformation.

<a name="d32"></a>
## D32 - The workflow is the deliverable · FIRM

Seven of nine Snakemake rules were stubs while every real result came from
hand-written sbatch scripts. For a repository whose stated purpose is
demonstrating that this **scales**, that inverted the deliverable: it showed
that an operator can hand-drive a membrane build, not that the pipeline runs
without one.

**What changed:** every stage now fans out over peptides as Snakemake
wildcards. The full pipeline is a **218-job DAG**; `stage4` alone is 172 jobs
across 10 peptides, and switching `md.peptides` from `benchmark` to `all` gives
325 jobs across 19. Scaling is a config change, not a workflow edit.

Paths are gone from the scripts - previously 14 of 28 hardcoded an absolute
scratch directory, so the repository could not be cloned and run anywhere else.

**Why this is the right emphasis:** a reviewer cannot verify a scaling claim
from a pile of shell scripts, which demonstrate only that a skilled person did
it once. A workflow file is a checkable assertion: run `snakemake -n`, see the
fan-out, change one line, see it grow. Stage 0 determinism is treated the same
way - a rule with an inspectable output rather than a sentence in a log.

**Retired to `scripts/contrib/`:** the coordinate-repair script whose constants
were tuned until one system worked, and six scripts superseded by the
parameterised Stage 3 entry points. Kept for history, off the main path, since
a reviewer finding a hand-rolled repair script with magic constants would
reasonably wonder what else was patched.
