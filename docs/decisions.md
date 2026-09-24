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
| [D19](#d19) | LICENSE copyright holder | OPEN | — |
| [D20](#d20) | OSU target architecture unconfirmed | OPEN | — |
| [D21](#d21) | ipSAE not produced by Chai-1 | OPEN | 2026-09-24 |
| [D22](#d22) | i-pLDDT taken as peptide-chain pLDDT | PROVISIONAL | 2026-09-24 |
| [D23](#d23) | Cross-seed agreement by receptor-superposed peptide RMSD | FIRM | 2026-09-24 |
| [D24](#d24) | Stage 2 test-peptide selection | FIRM | 2026-09-24 |

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
## D19 — LICENSE copyright holder · OPEN

`LICENSE` carries a placeholder. The user's full legal name is not known to me
and must not be guessed. **Blocks:** any distribution of the repository.

<a name="d20"></a>
## D20 — OSU target architecture · OPEN

ORCA (PSU, x86_64) is a **staging environment**; the stated target is a new
NVIDIA system at OSU. If Grace-based (GH200/GB200) its CPUs are **aarch64**,
which is why the instructions open with a `uname -m` check that looks
irrelevant here. **Not confirmed**, and drove [D2](#d2).
**Blocks:** a real porting risk assessment in `docs/arch_notes.md` — the
current table is x86_64 evidence and must not be read as aarch64 evidence.


<a name="d21"></a>
## D21 — ipSAE is not produced by Chai-1 · OPEN

The instructions require ipTM, **i-pLDDT, ipSAE**, aggregate score and PAE per
complex. Chai-1 provides `aggregate_score`, `ptm`, `iptm`, `per_chain_ptm`,
`per_chain_pair_iptm` and clash flags in `scores.model_idx_*.npz`; `plddt` and
`pae` come off the returned `StructureCandidates`. **ipSAE is not among them.**

ipSAE is a separate published interface metric derived from the PAE matrix, not
something Chai-1 emits. The column is written as `NOT_COMPUTED` rather than
omitted or filled with a substitute, so its absence is visible in the output
table rather than inferred from a missing column.

**Options:** implement it from the PAE matrix against the published definition
and validate that implementation; vendor the reference implementation; or
record it as out of scope with the reason stated. **Not yet decided — needs a
call before Stage 7**, since the instructions name it as a reported quantity.

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
