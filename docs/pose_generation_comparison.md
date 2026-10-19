# Four ways to place a peptide: comparing ADCP, FlexPepDock, Protenix and Chai-1

**Status: IN PROGRESS.** Purpose, design and decisions are recorded here as the
experiment is built. Results are filled in as each arm reports; sections marked
`PENDING` have not been run. Nothing in this file is a finding until it says so.

Started 2026-09-28. Prime directive from `instructions.md` applies: this is
**evidence of technical readiness, not a results campaign.** No arm is tuned
until it passes a control, and a negative result is a result.

## Why this experiment exists

Three scoring approaches have now failed the same control, documented in
`length_controls_result.md` and commit `56afbbe`:

| approach | control verdict |
|---|---|
| Chai-1 ipTM / ipSAE, single-sequence | fails - cannot discriminate at this peptide length |
| Chai-1 + template conditioning (arm B) | fails |
| MM/GBSA rescoring of MD poses | fails - size-biased; 14 kcal/mol GAFF2-vs-ff19SB spread on the *same* pose |

The common factor is not the scoring function. It is that **co-folding places
every peptide in the orthosteric pocket.** Every score therefore ranks
structures that are all already "bound", and the information that identifies a
non-binder has been discarded before scoring begins. None of the three tiers
beats the free baseline "does the sequence start with Tyr?" (AUROC 0.875 on the
real-vs-scramble comparison).

**The bottleneck is pose generation, not scoring.** This experiment tests that
claim by varying the pose generator while holding the evaluation fixed.

## The four arms

| arm | method | class | places non-binders outside the pocket? | status |
|---|---|---|---|---|
| **A** | Chai-1 0.6.1, single-sequence + ESM | co-folding (ML) | **No** - always builds a complex | done, fails |
| **B** | Protenix 2.0.0 `protenix_base_default_v1.0.0` | co-folding (ML) | **No** - same class | smoke passed, array pending |
| **C** | ADCP (AutoDock CrankPep) | physical docking, CPU, grid-based MC | **Yes in principle** - free search in a box | calibrating |
| **D** | Rosetta FlexPepDock in RosettaMP + franklin2019/2023 | physical refinement, membrane-aware | **No** - refinement needs a start pose | built, not run |

Arms A and B answer a narrow question: *is the failure specific to Chai-1, or
general to co-folding?* Neither addresses the pose-generation diagnosis, because
both place every peptide in the pocket by construction. **A Protenix pass would
not mean the pipeline works.**

Arms C and D are the ones that could in principle address the diagnosis, and
only C actually searches. This asymmetry is the point of the table and should
not be flattened in any summary of this experiment.

## Fixed evaluation, decided before any arm reported

Held constant so the four arms are comparable and so no arm can be selected on
the control:

- **Same 21 length controls**: beta-endorphin(1-10/15/20) as positives, 18 real
  milk-protein fragments with no N-terminal Tyr as negatives, 6 per length.
  Length varies in both arms, the opioid message in only one, so "too long to
  bind" and "lacks the pharmacophore" stay separable.
- **Same 5 seeds** (0-4) for the ML arms; fixed seed 42 for the docking arms.
- **`scripts/07b_length_controls_analysis.py` unmodified**, pointed at each
  arm's output. Arm B emits the identical `scores.tsv` schema so nothing
  downstream changes.
- **The baseline every arm must beat**: "does the sequence start with Tyr?",
  AUROC 0.875, free. An arm that does not beat it has not earned its compute.
- **Cross-seed / cross-replica agreement is reported for every arm.** In arm A
  this was the only signal that separated a correct pose from a failed one;
  ipTM and ipSAE could not.

## Arm B - Protenix: decisions

**D-P1. A separate conda environment (`mor-protenix`).** Protenix 2.0.0 pins
`numpy==2.4.1` and `torch==2.7.1`; `mor-chai` is `numpy<2`, `torch 2.6.0`. They
cannot share an environment, for the same reason `mor-chai` and `mor-pilot`
cannot - see `decisions.md` D5.

**D-P2. MSA-free single-sequence, matching arm A.** `infer_json_format.md`
states `pairedMsaPath`, `unpairedMsaPath` and `templatesPath` are all optional
and inference proceeds without them; `pred` also takes `--use_msa false`. So the
comparison is not confounded by one model receiving an MSA the other did not.

**D-P3. Matched inference budget is also Protenix's own default.** The run log
reports `Using default params for model protenix_base_default_v1.0.0: cycle=3,
step=200, use_msa=False`, and `--sample` defaults to 5. These are byte-identical
to arm A's recycles, timesteps and structures-per-run. Matching arm A therefore
does **not** handicap Protenix - a concern that was raised and dismissed on this
evidence.

**D-P4. The 2021-09-30-cutoff model, deliberately.** This decides what a pass
would mean:

| structure | ligand | released | in training? |
|---|---|---|---|
| 6DDF | DAMGO | 2018-06-13 | **yes** |
| 8F7Q | beta-endorphin | 2022-12-14 | **no - held out** |

So memorising the exact beta-endorphin complex is excluded, but the
YGGF-in-the-MOR-pocket motif is learned via 6DDF. A pass is informative with
that residual caveat. Using `protenix_base_20250630_v1.0.0` instead would put
8F7Q *into* training and destroy this property.

**D-P5. The fused kernel is pre-compiled once, before any array.** Protenix
JIT-compiles `fast_layer_norm_cuda_v2` at *import* time
(`protenix/model/layer_norm/layer_norm.py`, `try: import / except ImportError:
compile`), with no environment switch to skip it and no CPU fallback. It needs
`CUDA_HOME` (`module load cuda/12.9.0-none-none`), takes 188 s, and
`build_directory` resolves to **site-packages** - a shared path. 45 concurrent
array tasks would race `ninja` into one directory. Pre-compiling makes every
task hit the plain `import` branch.

**Measured cost:** 1201 s per (peptide, seed) on one L40S => the 105-task array
is ~35 GPU-hours, ~1 h wall at 45 concurrent, and each task fits `short`'s 4 h
limit.

## Arm C - ADCP: decisions and one discarded result

**D-C1. Fixed search effort across every peptide, not scaled with length.**
Scaling effort with length would search long peptides harder than short ones,
and length is the variable under test - any length trend could then be an
artefact of the sampling schedule. The spread across replicas is recorded
instead as a convergence check.

**D-C2. Box sized to the measured pocket, not to DAMGO's footprint.** Centred on
the 14 residues contacting DAMGO in >=50% of frames of the validated 50 ns run,
padded 11 A: 48.0 x 39.8 x 44.3 A, roughly double the measured pocket
(26.0 x 18.1 x 21.9 A) in every dimension. AGFR's default `-b ligand` would
have built the box around DAMGO's own 5.9 x 10.1 x 10.6 A footprint - under a
quarter of the pocket volume - forcing 15- and 20-mers to compress and
producing bad scores for long peptides as a pure box artefact.

**D-C3. ADCP requires a RELATIVE `-t` path.** It does
`os.path.join('./tmp_<job>/', target[:-4], 'translationPoints.npy')`, and
`os.path.join` discards everything before an absolute path, so an absolute `-t`
collapses the unzip directory away. This killed array `191241` - all 61 tasks -
while SLURM reported `COMPLETED`. Fix: symlink the `.trg` in and pass the bare
name.

**D-C4. `Constraining amino acids: 1` is NORMAL. Do not "fix" it.** ADCP writes
the `constrains` file itself, unconditionally, at `runADCP.py:182-184`, and
wires it into its own default option string
`-p Bias=NULL,external=5,constrains,1.0,1.0` at `:226`. It was initially
misdiagnosed here as a stray file left by the `-y` dry run. It is not.

### DISCARDED: the first production array (191302)

Array `191302` produced 61 peptides x 30 replicas = 1830 energies at
**`-n 100000`**. ADCP's own default is **2,500,000** (`runADCP.py:106`, argparse
default `:431`), with `-N 50` against our 30. **The assay ran at 1/25th of the
authors' sampling.** The symptoms were all one defect:

- replica spread reached 300-430 units at 15-20 aa vs ~8 at 4 aa (rho=+0.937
  with length); 34/61 peptides had a spread larger than their own best energy
- 268/1830 replicas returned **positive** energies, 157 above +50, max +417 -
  ADCP's 99999 / 2e6 failure sentinels surfacing as `totalE` because the search
  never escaped clashes
- replicas "finished" in a median of **5 s**

The discrimination statistics computed from it (length-control AUROC 0.463,
real-vs-scramble AUROC 0.488, medians -18.361 vs -18.415, plus an apparent
"longer peptides score worse" trend rho=+0.373) **measure sampling failure, not
binding, and are discarded.** They are recorded here only so they are never
mistaken for a fourth failed approach. The under-sampling penalty grows with
length, which manufactures precisely that spurious length trend.

`scripts/08_adcp_analysis.py` and `docs/adcp_energies.tsv` were written against
this run and are retained, because the analysis is correct and the rerun is a
one-command repeat.

### Calibration (array 191377): best energy has NOT converged at ADCP's default

3 peptides x 3 step budgets x 30 replicas. `nonneg` = replicas returning a
non-negative (sentinel) energy.

| peptide | budget | best | spread | nonneg/30 | wall |
|---|---|---|---|---|---|
| YGLF (4 aa) | 100k | -16.83 | 5.9 | 0 | 21 s |
| YGLF (4 aa) | 500k | -21.44 | 8.0 | 0 | 42 s |
| YGLF (4 aa) | **2.5M** | **-22.72** | 6.6 | 0 | 103 s |
| beta-end(1-10) | 100k | -15.18 | 122.6 | 3 | 32 s |
| beta-end(1-10) | 500k | -24.92 | 118.9 | 1 | 78 s |
| beta-end(1-10) | **2.5M** | **-29.16** | 121.2 | 1 | 198 s |
| beta-end(1-20) | 100k | -11.32 | 297.5 | **20** | 56 s |
| beta-end(1-20) | 500k | -24.19 | 250.5 | **17** | 173 s |
| beta-end(1-20) | **2.5M** | **-35.33** | 173.3 | **7** | 482 s |

Four things this establishes:

1. **The 100k run was badly wrong, not marginally wrong.** The 10-mer positive
   gains 14.0 units (92%) from 100k to 2.5M; the 20-mer gains 12.9 units by
   500k alone. At 100k the 10-mer positive scored -15.18 against a negatives
   median of -17.97, i.e. it *looked like a non-binder*. Properly sampled it is
   far better than that median. **Any ranking from the 100k run is void.**
2. **Best energy has not plateaued even at ADCP's default 2.5M** for 10 and
   20 aa. Gains per step increase: 4-mer 4.6 then 1.3 (converged); 10-mer 9.7
   then 4.2; 20-mer 12.9 then **11.1 - barely decelerating**. Sentinel replicas
   do fall with effort (20 -> 17 -> 7 for the 20-mer), so the search is
   genuinely improving rather than saturating.
3. **A monotonic size bias EMERGES as sampling improves, and the 100k run had
   it inverted.** At 2.5M the best energies order by length:
   20 aa (-35.33) < 10 aa (-29.16) < 4 aa (-22.72) - longer peptides score
   better simply because they have more atoms making more contacts. This is the
   same size bias measured for MM/GBSA (`mmgbsa_validation.md`). At 100k the
   apparent trend ran the *other* way (rho=+0.373, "longer scores worse")
   purely because the under-sampling penalty grows with length.

   **Consequence: raw ADCP `best target energy` is NOT comparable across
   lengths.** Only the per-length control structure - positive vs 6
   length-matched negatives - is a valid comparison, and the real-vs-scramble
   comparison is valid because scrambles are length-matched by construction.
   Any cross-length ranking of the real library on raw energy would be reading
   peptide size.
4. **The 20-mer still fails in 7/30 replicas at 2.5M** (a majority, 17/30, at
   500k). Consistent with ADCP's published degradation beyond ~12 residues.
   Reported as a capability limit, not tuned away.

**Consequence for the pre-registered criterion.** The criterion written into
`jobs/adcp_calibrate.sbatch` - spread a small fraction of |best| **and** zero
non-negative replicas - cannot be met at 20 aa at any budget tested. It is
therefore revised, and the revision is recorded rather than quietly applied:

- convergence is judged on the **spread among valid (negative-energy) replicas**
  and on **best-energy plateau**, not on raw spread, which is dominated by
  sentinels
- the **fraction of valid replicas is reported as a primary result per length**,
  because "ADCP does not converge for 20-mers at a practical budget" is a
  legitimate finding about cost and capability
- non-negative energies are **rejected explicitly**, never averaged in

## Arm D - FlexPepDock in RosettaMP: assessment before building

The proposal is FlexPepDock inside RosettaMP with `franklin2019` or
`franklin2023` - peptide-appropriate sampling with membrane-aware scoring. The
combination appears to be unbenchmarked, so it would need validating here.

**Constructible: yes.** `AddMembraneMover` and `FlexPepDockingProtocolCreator`
are both registered RosettaScripts movers (confirmed in
`init.MoverCreators.ihh`), and `franklin2018/2019/2021/2023.wts` all ship in the
database.

**Three obstacles:**

1. **`rosetta_scripts` was not built.** Job `191238` built only the
   `FlexPepDocking` app. The XML route needs the RosettaScripts binary.
2. **The FlexPepDocking app contains no membrane code at all** - no
   `AddMembraneMover`, no franklin reference, in neither the app nor
   `protocols/flexpep_docking/`. So `FlexPepDocking -score:weights franklin2019`
   loads the weights against a pose with **no `MembraneInfo`**, and
   `fa_water_to_bilayer` becomes a silent no-op or throws. **This is the most
   dangerous failure mode in the proposal: the output would look
   membrane-aware and would not be.** Any run must assert that a membrane was
   actually instantiated.
3. **A spanfile is required** for `AddMembraneMover` (receptor TM topology).

**Two scientific caveats:**

- **The membrane terms may not touch the interface.** franklin2019's only
  membrane term is `fa_water_to_bilayer` (weight 0.5) - a per-residue
  water-to-bilayer transfer energy; franklin2023 adds `f_elec_lipidlayer 0.128`
  and `fa_imm_elec 0.01`. The mu-OR orthosteric pocket sits in the extracellular
  vestibule, above the bilayer hydrophobic core. If the peptide binds above the
  phosphate plane, these terms contribute little *differential* signal between
  peptides, and their real value is keeping the TM bundle from distorting.
  That is receptor realism, not discrimination. **Cheap pre-check, to run before
  building anything: measure the bilayer-depth distribution of the peptide's
  contact residues in the existing validated 50 ns trajectory.**
- **FlexPepDock is refinement, not global search.** It needs a start within a
  few A, so it inherits the pose-generation bottleneck and cannot place a
  non-binder outside the pocket unless explicitly allowed. Fine for
  self-docking validation; not a fix for discrimination.

**D-D1. First validation target is 8F7Q / beta-endorphin, not DAMGO.** DAMGO
carries `DAL` (D-Ala), `MEA` (N-methyl-Phe) and `ETA` (C-terminal alcohol);
Rosetta handles D-amino acids natively but `MEA` and `ETA` need custom params -
see `damgo_notes.md`. beta-endorphin in 8F7Q is all-canonical **and already has a
Chai-1 pose-accuracy baseline in this repo**, giving a direct cross-method
comparison for free.

**D-D2. Unresolved contradiction, to be settled deliberately.**
`instructions.md:377` and `README.md:305` exclude ref2015-based tools **by
design**, for exactly the missing-membrane reason this arm proposes to fix. A
Rosetta build has already happened anyway. Either the exclusion is amended with
its reasoning, or this arm is dropped - it must not be resolved by silence.

## A note on exit codes

Three separate broken runs in this experiment reported SLURM `COMPLETED` /
exit `0:0`: ADCP array `191241` (dead on the absolute-path bug, zero output),
the first Protenix smoke test (`CUDA_HOME` unset, no inference, no output
directory), and ADCP array `191302` (ran at 1/25th sampling). The cause in each
case was a script that recorded a tool's return code and carried on, or piped a
tool through `tail` and let the pipeline's last command set the status.

This matters beyond tidiness: the arm A write-up leans on the sentence "all 105
jobs completed, so there is no computational failure at these lengths." That
sentence is only load-bearing if `COMPLETED` means something. Every job script
in this experiment now fails loudly, and `b48ce60` records the same false-report
pattern for PDBFixer.

## Results

| arm | length-control AUROC | real-vs-scramble AUROC | beats Tyr baseline (0.875)? | notes |
|---|---|---|---|---|
| A Chai-1 | 0.463 (ipTM, n=3 vs 18) | 0.488 | **no** | 21/21 DIVERGENT; `length_controls_result.md` |
| **B Protenix v2 (+MSA)** | 0.815 (p=0.101) | **0.290** (p=0.024) | **no** | prefers scrambles *significantly*; see below |
| B Protenix v1 | ~~0.278~~ **void** | - | - | receptor never folded; see below |
| C ADCP v2 | **0.463** (p=0.887) | **0.562** (p=0.508) | **no** | converged; see below |
| C ADCP v1 | ~~0.463~~ **void** | ~~0.488~~ **void** | - | 1/25th sampling; discarded |
| D FlexPepDock | `PENDING` | | | pipeline validated, funnel test pending |

### Arm B v2 result: the first arm to discriminate - with four live caveats

105/105 tasks, 0 FATALs, `N_msa` 8862-8876 on every task. Receptor accuracy holds
across the set, not just at the gate peptide: **pocket CA RMSD 1.02-1.14 A**
(all-CA 2.79-2.87) against arm A's 0.88-0.93 A pocket. The receptor is folded and
the orthosteric site is right.

| length | negatives median ipTM | POSITIVE | negatives >= positive |
|---|---|---|---|
| 10 aa | 0.8318 | **0.8619** | **0 / 6** |
| 15 aa | 0.7230 | **0.7868** | 1 / 6 |
| 20 aa | 0.6422 | **0.8274** | **0 / 6** |

**Length-control AUROC 0.815** (p=0.101, n=3 vs 18, descriptive). The positive
outscores every length-matched negative at 10 and 20 aa, and 5 of 6 at 15 aa.
ipTM now sits at 0.62-0.86, against 0.17-0.32 in the void v1 run and 0.21-0.35 in
arm A - consistent with a model that has an actual pocket to dock into.

For the first time in this pilot, a co-folding score separates the opioid message
from length-matched non-openers. That is a real change and it is worth stating
plainly. It is also not yet a working triage tier, for four reasons:

**1. It does not beat the free baseline, and on THIS control set the baseline is
unbeatable.** The length controls were built with positives that begin `YGGFM`
and negatives that begin with anything but Tyr - by construction. So "does the
sequence start with Tyr?" scores **AUROC 1.000** here, against ipTM's 0.815. The
baseline that mattered in `length_controls_result.md` was measured on
real-vs-scramble, where scrambles retain the tyrosine; arm B has **not** been run
on that set, so the comparison that could actually favour ipTM has not been made.
Running arm B against the 40-member real/scramble library is the next step, and
until then "beats the baseline" is unanswered, not answered negatively.

**2. ipTM is severely length-confounded.** rho = **-0.879** within negatives
alone (p=1.6e-06), -0.771 across all controls. Negative ipTM-vs-length means
longer peptides score worse regardless of sequence. The discrimination above
comes largely from the positive *failing to degrade* with length (0.862 -> 0.787
-> 0.827, flat) while negatives collapse (0.832 -> 0.723 -> 0.642). That is a
genuine signal, but it is a signal about length-robustness, not affinity, and the
per-length structure of the control set is the only reason it is visible at all.
Any cross-length ranking of a real library on ipTM would be reading length.

**3. High confidence sits on top of unreliable poses.** 18 of 21 peptides are
still DIVERGENT on cross-seed agreement (1 CONVERGED, 2 PARTIAL), and the worst
case is the strongest result: the 20-aa positive scores ipTM **0.827** with a
cross-seed peptide RMSD of **10.40 A**. `02_agreement.py` was written for exactly
this: "a high ipTM with poor cross-seed agreement is a stronger warning than
either number alone." The model is confidently producing different answers. Arm B
v2 improved agreement a great deal over v1 (which was 9.5-30 A, all DIVERGENT)
but has not fixed it.

**4. The memorisation caveat is now load-bearing.** It did not matter while every
arm failed; it matters now that one has succeeded. `protenix_base_default_v1.0.0`
has a 2021-09-30 cutoff, so **8F7Q (beta-endorphin, released 2022-12-14) is held
out - but 6DDF (DAMGO, 2018-06-13) is IN training.** The positives are
beta-endorphin N-terminal fragments beginning `YGGFM`, i.e. met-enkephalin, and
DAMGO is an enkephalin analogue bound to this receptor in a training structure.
So the model has seen a `YGGF`-motif peptide in the mu-OR orthosteric pocket.
A held-out test that separates learned pharmacophore recognition from recalled
complex geometry is now the single most important missing experiment in this arm.

### Arm B v2 on the fair baseline test: significantly WORSE than chance

200/200 tasks (40 peptides x 5 seeds), 0 FATALs. 20 real peptides vs 20
**exactly length-matched** scrambles - identical length histograms, so the
ipTM-vs-length confound cannot drive this comparison.

| metric | real median | scramble median | AUROC | p |
|---|---|---|---|---|
| ipTM | +0.8849 | **+0.9096** | **0.290** | 0.024 |
| ranking_score | +0.8912 | **+0.9122** | **0.285** | 0.021 |
| gpde interface | +1.2286 | **+1.0535** | **0.265** | 0.011 |

| | AUROC |
|---|---|
| "starts with Tyr" baseline | **0.875** |
| arm C ADCP, same 40 peptides | 0.562 |
| **arm B best metric** | **0.290** |

**All three of Protenix's own confidence metrics prefer the scrambles, and the
preference is statistically significant** (p=0.011-0.024). Seven of the top eight
peptides by ipTM are scrambles; the only real peptide in that group is
met-enkephalin. This is not a null result - the score is systematically
misled, in a consistent direction, across three independent metrics.

**The scramble caveat does not explain it.** If scrambles were weak binders
rather than non-binders (the ambiguity `length_controls_result.md` was built to
avoid) the expectation would be AUROC near 0.5, not 0.29. A systematic
*preference* for scrambles needs a different explanation. The most likely one:
ipTM rewards a confidently foldable, compact complex, and scrambling can produce
a sequence the model packs more happily than a real opioid peptide, which is
extended and conformationally constrained in the orthosteric site. On that
reading ipTM is measuring "can I build a tidy complex from this sequence", which
is not binding and is not even monotonically related to it.

**This reinterprets arm B's length-control result.** Taken alone, AUROC 0.815
looked like the first genuine signal in the pilot. Taken with 0.290 here, the
length-control number is almost certainly the **length-robustness artefact**
already visible in it: the positives happened not to degrade with length
(0.862 -> 0.787 -> 0.827) while negatives collapsed (0.832 -> 0.723 -> 0.642).
Remove the length difference, as this set does, and the signal does not merely
vanish - it inverts. The 0.815 should not be cited as discrimination.

**A prediction made before this run, which was wrong and is recorded as such:**
ipTM was expected to be "high and flat" at 4-10 aa, on the strength of arm A's
finding that ipTM is structurally unable to discriminate at short peptide
length. It is neither - it spans 0.247-0.958 with sd 0.168. ipTM at these
lengths is highly variable and confidently wrong, which is a worse failure mode
than being uninformative.

Cross-seed agreement: 5 CONVERGED, 18 PARTIAL, 17 DIVERGENT. Better than the
length controls, and it makes no difference to the outcome.

### Two defects in the control sets, found by inspection not by a run

**1. The `real_peptide` class violates D30, and both arm B's real40 run and arm C
inherit it.** `docs/reference_peptides.tsv` carries a `role` column and a
`benchmark_include` flag: 8 of its entries are **antagonists** (casoxin A/B/C/D,
casoxin_B_hum, lactoferroxin_A bov/hum, CTOP), all marked
`benchmark_include=NO`, and only 7 entries are `benchmark_include=YES` (all
agonists). But `data/templates/dock_set.tsv`'s `real_peptide` class contains
**6 of those antagonists**. Decision D30 states "Benchmark restricted to
agonists; receptor is active-state only" - antagonists stabilise the INACTIVE
conformation, so docking them into an active-state receptor is a conformational
mismatch. Three soymorphins in the same class are not in
`reference_peptides.tsv` at all, so they are positives without provenance, which
conflicts with the standing rule that every affinity be traceable.

**It does not explain the inverted AUROC.** Tested directly:

| subset | n | AUROC | p |
|---|---|---|---|
| all 20 real (as run) | 20v20 | 0.290 | 0.024 |
| minus the 6 antagonists (D30-compliant) | 14v14 | 0.286 | 0.057 |
| `benchmark_include=YES` agonists only | 7v7 | **0.265** | 0.165 |
| the 6 antagonists alone | 6v6 | 0.361 | 0.485 |

The effect size is stable (0.290 -> 0.286 -> 0.265; p weakens only because n
falls to 7), and antagonists do not score below agonists within the real class
(ipTM medians 0.8832 vs 0.8897, p=0.547). So the active-state receptor is not
penalising them detectably. **Restricted to the seven peptides the repo itself
designated benchmark-worthy agonists, Protenix still prefers the scrambles.**
The class must be cleaned for D30 compliance, but the conclusion stands.

**2. There is no composition-matched control at 10-20 aa.** The two sets are
disjoint by design:

| set | lengths | positives | negatives |
|---|---|---|---|
| length controls | 10-20 aa | 3 beta-endorphin fragments | 18 real non-Tyr milk fragments |
| real40 | 4-10 aa | 20 real food-derived peptides | 20 composition-matched scrambles |

So **scrambles and random sequences exist only at 4-10 aa**, and the 10-20 aa
negatives are real protein fragments selected for lacking an N-terminal Tyr -
not composition-matched, not random. This matters because 4-10 aa is precisely
the regime where arm A concluded ipTM is structurally unable to discriminate,
and because the whole length-robustness reading of arm B's 0.815 rests on the
10-20 aa set, where no scramble control exists to test it. The missing arm is
in-place scrambles of the 21 length controls: composition- and length-matched
negatives in the regime that set actually probes (~105 tasks, ~13 GPU-hours).

### Revised diagnosis: the bottleneck is not pose generation

Commit `56afbbe` concluded that co-folding's failure was a **pose-generation**
problem - every peptide is placed in the orthosteric pocket, so the information
identifying a non-binder is destroyed before scoring. Four arms later that
diagnosis does not survive:

- **ADCP supplied the missing capability and changed nothing.** It searches a
  48 x 40 x 44 A box freely and *can* place a peptide outside the pocket. With
  converged sampling (rank stability rho=+0.867) it scores AUROC 0.463 on the
  length controls and 0.562 on real-vs-scramble. Free global search did not
  produce discrimination.
- **Protenix with a correct pocket (1.03 A) is significantly anti-predictive.**
  The receptor is right, the poses are placed, and the score prefers scrambles.
- **Every method is beaten by a free one-character rule.** "Starts with Tyr"
  scores 0.875 on the set where it is beatable. No arm exceeds 0.562.

What the four arms share is not a pose generator. It is that **all of them score
a complex that has already been built, with functions that reward geometric and
packing plausibility** - learned confidence (arms A, B), a docking energy
dominated by peptide size (arm C, rho=-0.718), and a Rosetta total score
dominated by receptor repacking (arm D, `total_score` vs RMSD rho=+0.04).
Plausibility of the assembled complex is simply not the quantity that separates
a binder from a non-binder at this peptide size, and no amount of better
sampling or better folding supplies it.

That is a finding about the whole approach, and it is the pilot's most
defensible result. It also means the honest deliverable is a **negative
technical-readiness report with a quantified free baseline**, not a triage tier.

### VOID: arm B v1. Decision D-P2 was wrong, and the receptor never folded

All 105 tasks ran clean (array 191466, 0 FATALs, 525 structures, ~70 s each) and
the analysis produced numbers - length-control AUROC 0.278, 6/6 negatives above
the positive at 15 aa, 21/21 DIVERGENT with cross-seed peptide RMSD 9.5-30.0 A
against arm A's 2.7-23.6 A. **Those numbers are void.**

**The receptor was never correctly predicted.** Superposed on 8F7Q chain R:

| | receptor CA RMSD vs experiment |
|---|---|
| arm A (Chai-1) | **2.72 A** over 281 CA |
| arm B v1 (Protenix) | **16.16 / 16.37 A** over 281 CA |

At 16 A over 281 CA the 7TM bundle is not folded, so there was no orthosteric
pocket for any peptide to bind. The run measures a broken input, not Protenix's
ability to discriminate.

**Root cause: decision D-P2 is false.** It claimed MSA-free Protenix matched
arm A's single-sequence mode, so the comparison was unconfounded. Arm A actually
runs `use_msa_server=False` **with `use_esm_embeddings=True`** - Chai-1
substitutes an ESM language-model embedding for the MSA. Protenix with
`--use_msa false` has no MSA **and no such surrogate**. "Single-sequence" meant
two different things, and arm B was crippled rather than matched. Its low
confidences were visible from the pilot onward (ipTM ~0.2, plddt_complex 0.454)
and should have been treated as a red flag before the array was launched.

**This is the second time in this experiment that a "matched setting" was not
matched and produced a void negative result** - the first being ADCP at 1/25th
sampling. Both were caught by a control, not by the run failing. The pattern to
carry forward: a matched *flag name* is not a matched *condition*, and an arm
that reports implausibly low confidence should be validated against a known
structure before its discrimination statistics are computed at all.

**Arm B must be re-run with a real MSA** (`protenix msa` / `prep`, or template
conditioning against the experimental receptor) and gated on reproducing the
receptor to within a few A before any discrimination number is read from it.
Note that Protenix's MSA server modes send the receptor sequence to a
third-party service, which is a decision to take deliberately.

### Arm C result: ADCP fails the control at its authors' own sampling

61 peptides x 50 replicas at `-n 2500000` (array 191386), 3050 replicas, **4.0%
sentinel rejection against v1's 14.6%**, 37-50 valid replicas per peptide.

| length | negatives median | POSITIVE | negatives scoring better |
|---|---|---|---|
| 10 aa | -32.368 | -29.429 | **5 / 6** |
| 15 aa | -34.904 | -33.339 | 3 / 6 |
| 20 aa | -34.521 | **-36.429** | 2 / 6 |

- length controls: **AUROC 0.463, p=0.887** - no discrimination
- real vs scramble: **AUROC 0.562, p=0.508**, medians -26.425 vs -25.062
- the free Tyr baseline scores **0.875**, so ADCP does **not** earn its CPU time

**The ranking IS converged, so this is a real negative result and not another
sampling artefact.** The rank-stability check (array 191387) re-ran the 9-member
10-aa group at 10M steps, 4x production: **Spearman rho +0.867, p=0.0025**
against the 2.5M ranking, with a mean energy gain of only -3.85. Absolute
energies keep drifting; the ordering does not. That was the pre-registered
criterion, revised and met.

**Size bias is now the dominant signal: rho = -0.718, p=7.1e-11** (longer =
better energy), and it holds *within* homogeneous classes - real peptides
rho=-0.603 (p=0.005), scrambles rho=-0.658 (p=0.002). So ADCP's raw energy
ranks the real library by peptide length. The per-length controls exist to
strip that out, and with it stripped out there is no signal left.

**Why this matters more than the three earlier failures.** ADCP is the first arm
that *could* place a non-binder outside the pocket - it searches a 48 x 40 x 44 A
box freely rather than building a complex by construction. The pose-generation
diagnosis from `56afbbe` predicted that this was the missing capability. It was
supplied, and discrimination did not follow. So "co-folding puts everything in
the pocket" is **not a sufficient explanation** for the pilot's failure to rank
peptides, and the diagnosis in `56afbbe` needs revising rather than restating.

### Arm D interim: the ranking metric cannot be `total_score`

Pipeline validated on 8F7Q (job 191459): prepack 40 s, 10 decoys in 748 s
(~79 s/decoy), ref2015, `min_receptor_bb false`.

| metric vs RMSD | rho | p |
|---|---|---|
| `total_score` vs `rmsBB_if` | +0.224 | 0.533 |
| `total_score` vs `rmsBB` | +0.042 | 0.907 |
| **`I_sc` vs `rmsBB`** | **+0.794** | **0.006** |
| `I_sc` vs `rmsALL_if` | +0.685 | 0.029 |

`total_score` spans -390 to -521 REU, and that ~130 REU of variation is
receptor self-energy from repacking 287 residues; the peptide interface
contributes only -46 to -63 REU. The best-scoring decoy (-521.3, 2.00 A) beats
the native-most decoy (-469.9, 1.84 A) by 51 REU **in favour of the worse
pose**. Ranking on `total_score` would therefore measure receptor repacking, not
binding - the same class of error as MM/GBSA's size bias. **Arm D's primary
metric is `I_sc`.**

Caveat: n=10 and every decoy started at the native backbone
(`startRMSbb = 0`), so this is drift behaviour, not a funnel. A funnel test with
perturbed starts (`random_rot_start` / `random_trans_start`) is still required
before arm D reports.

## Met-enkephalin vs five non-opioid peptides, on all three pipelines

User-supplied set, 2026-09-28: `REKVLASS`, `FKCRRWQWRMK`, `QNALIVRYTR`,
`EGPKLVAS`, `DPYKLRP` (five, not six as stated). All runs clean - Protenix
25/25, Chai-1 25/25, ADCP 250 replicas, 0 FATALs. FlexPepDock cannot
participate: it is a refiner and has no starting pose for these.

**Provenance caveats, recorded in `data/templates/random5.tsv` itself:**
`FKCRRWQWRMK` is **lactoferricin B** (N-terminal 11 residues of the bovine
lactoferrin 17-41 fragment) - a characterised cationic antimicrobial peptide,
milk-derived, +5 charge, 2 Trp. It is not a random sequence. `QNALIVRYTR`
shares `IVRYT` with the receptor's own TM1 (`...LVMYVIVRYTKMK...`).

### The three pipelines contradict each other

| pipeline | met-enkephalin rank | 7 agonists vs 5 randoms | p |
|---|---|---|---|
| Protenix (+MSA) | **1 of 6** (ipTM 0.9545) | AUROC 0.657 | 0.432 |
| Chai-1 | **3 of 6** (ipTM 0.3090) | AUROC 0.486 | 1.000 |
| ADCP | **5 of 6** (E -24.34) | AUROC **0.114** | **0.030** |

Protenix calls met-enkephalin the best of the six. ADCP calls it second-worst,
and on the powered comparison ADCP is **significantly worse than chance**
(AUROC 0.114, p=0.030 - it prefers the non-opioid peptides). Chai-1 sees nothing.

### Why, and it is the same artefact in both directions

Score correlates with peptide LENGTH, not with opioid identity:

| pipeline | rho(score, length), n=6 | direction |
|---|---|---|
| Protenix ipTM | **-0.754** (p=0.084) | shorter -> higher confidence |
| ADCP energy | **-0.754** (p=0.084) | longer -> more negative energy |
| Chai-1 ipTM | -0.551 (p=0.257) | weak |

Within the five randoms alone the ordering is close to monotonic in length -
ADCP rho=-0.821, Protenix rho=-0.564. **Met-enkephalin is the shortest peptide
in the set at 5 aa**, and each pipeline's verdict is exactly what its own length
bias predicts: Protenix rewards short and ranks it first; ADCP's size bias
(rho=-0.718 across the full 61-peptide set) rewards long and ranks it fifth.

ADCP's top-scoring peptide is **lactoferricin B at -35.59** - the longest
peptide in the set, an antimicrobial with no opioid activity.

**Conclusion: this set cannot support the assertion.** Met-enkephalin's rank
flips from best to second-worst depending only on which length artefact the
method carries, so no pipeline is providing evidence about opioid binding here.

### Two design limits that would have to be fixed to ask this properly

1. **Length is not matched.** Met-enkephalin is 5 aa; the five comparators are
   7-11 aa. Given both non-null pipelines are length-driven, the comparison is
   confounded at the outset and a met-enkephalin "win" is uninterpretable. A
   length-extrapolation test on the 7 benchmark agonists was attempted and is
   inconclusive - the agonist-only fit has a positive slope (+0.163/aa, p=0.169,
   n=7, spanning only L=4-7), so it cannot be used to predict L=5 reliably. The
   fix is comparators at 5 aa, not a correction applied afterwards.
2. **None of the five begins with Tyr**, so the free "starts with Tyr" rule
   scores **AUROC 1.000** on this set. Even Protenix's rank-1 result only
   *matches* a rule that costs nothing. Testing whether any method beats the
   baseline requires non-opioid peptides that DO start with Tyr - which the
   existing scrambles provide, and where Protenix scored 0.290.

**Incidental observation worth keeping:** among the five supposedly equivalent
non-binders, Protenix's ipTM ranges from 0.2006 (`REKVLASS`) to 0.8823
(`DPYKLRP`) - a 0.68 spread across peptides that should all be equally
non-binding. Whatever ipTM is responding to, it varies enormously between
sequences with the same (absent) pharmacology.
