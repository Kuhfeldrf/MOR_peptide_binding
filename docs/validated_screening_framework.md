# A screening framework for MOR peptides that actually validates

**Question posed 2026-09-28:** find a scientifically valid, literature-documented
framework for screening opioid peptides against the mu-opioid receptor, after
five structure-based approaches failed to separate peptides with documented
opioid activity from likely non-binders.

**Answer: the framework that works is a sequence pharmacophore, not a
structure-based score.** The Schwyzer message-address model - N-terminal
tyrosine with a free alpha-amine plus an aromatic residue at position 3/4 -
reaches **AUROC 0.967** on a length-balanced benchmark of MOR binders versus
other GPCRs' cognate peptide ligands, and **0.971 under length stratification**.
Every structure-based method tested reaches at most 0.614 on the same set, and
**0.500 once peptide length is controlled**. The recommendation is therefore to
use the pharmacophore as the triage tier and to stop treating co-folding
confidence or docking energy as a ranking signal at this peptide size.

## 1. What the literature actually validates

Two papers reframed the problem. Neither supports the task as this pilot had
been posing it.

**Zhang et al., "Deep learning in GPCR drug discovery: benchmarking the path to
accurate peptide binding", *Briefings in Bioinformatics* 26(2) bbaf186, 2025.**
The published benchmark is **cognate ligand vs decoys drawn from OTHER
receptors** - 124 GPCR-peptide pairs, 1240 decoys at 10:1, decoys being either
<35% sequence identity or ligands of homologous receptors. Reported: AF2
AUC 0.86, AF3 0.82, **Chai-1 0.76**. Two details matter and are easy to miss:
the headline AUC is **aggregated over 124 receptors**, and per-receptor top-1
performance is only **58%**. A single-receptor campaign should expect far less
than 0.86.

**Kim, A.R. et al., "Enhanced Protein-Protein Interaction Discovery via
AlphaFold-Multimer", bioRxiv 2024.02.19.580970; https://github.com/flyark/AFM-LIS.**
Global confidence metrics dilute a small interface across a large complex. The
Local Interaction Score family fixes this: **LIA** = inter-chain residue pairs
with PAE <= 12 A; **LIS** = mean of those PAE values inverted to 0-1; **cLIS** =
same restricted to direct contacts (Cb-Cb <= 8 A); **iLIS = sqrt(LIS x cLIS)**.
The authors report it beats global metrics "in every pLDDT subgroup" and rescued
47 interactions ipTM missed, with gains concentrated in "flexible and small
interfaces". A 4-20mer peptide in a GPCR vestibule is exactly that, and ipTM
averaged over a 281-residue receptor is exactly the dilution being described.

Implemented in `scripts/02g_lis_rescore.py`. One honesty note: neither the paper
nor the repository gives the PAE->score arithmetic beyond "inverted to a scale of
0 to 1", so the standard reading `(12 - PAE)/12` is used. A different
normalisation would change absolute values but not the ranking.

## 2. The framework

Formalised from the message-address model (Schwyzer 1977; Portoghese). Rules
were taken from prior pharmacology, **not fitted to this data**:

| rule | definition |
|---|---|
| **R1** | `Tyr` at position 1 |
| **R2** | `Tyr1` **and** aromatic (F/W/Y) at position 3 **or** 4 |
| **R3** | `Tyr1` and aromatic at position 4 |
| **R4** | `YGGF...` or `YP.F...` |

R2 is recommended. All four beat every structure-based method on every
benchmark, so the conclusion does not depend on which rule is chosen - worth
stating because R2 was selected after seeing four rules on one set, which is
mild selection and would otherwise be a fair objection.

## 3. Validation

### Benchmark A - the documented task (`data/reference/cognate_benchmark.tsv`)

29 peptides: **14 MOR binders** (enkephalins, endomorphins, casomorphins,
dynorphins, beta-endorphin fragments) vs **15 cognate ligands of other GPCRs**
(angiotensin II, angiotensin-(1-7), bradykinin, substance P, alpha-MSH, CCK-8,
neurotensin-(8-13), galanin, apelin-12, kisspeptin-10, NPY, PYY, nociceptin
x2). Deliberate properties:

- **Length-balanced** (Mann-Whitney p=0.282, medians 9 vs 10). Built at p=0.081
  first, then beta-endorphin(1-15) and dynorphin A(1-17) were added to correct
  it - because every metric here favours short peptides, that imbalance would
  have manufactured a win.
- **Four decoys begin with Tyr.** NPY(1-10) and PYY(1-10) begin `YP`, exactly
  like the casomorphins and endomorphins. So R1 is beatable at 0.867 rather
  than trivially perfect.
- **Nociceptin is the sharpest control**: opioid-family, NOP-selective, and
  differs from met-enkephalin at position 1 only by Y->F.

| method | AUROC | p | length-controlled |
|---|---|---|---|
| **pharmacophore R2** | **0.967** | <1e-5 | **0.971** |
| pharmacophore R3 / R4 | 0.929 | <1e-5 | - |
| pharmacophore R1 (Tyr1) | 0.867 | 7e-5 | - |
| Chai-1 **iLIS** | 0.614 | 0.305 | **0.500** |
| Chai-1 aggregate_score | 0.576 | 0.499 | - |
| Chai-1 ipTM | 0.567 | 0.556 | - |
| Chai-1 interface PAE min | 0.538 | 0.743 | - |

R2 caught **14/14 binders with one false positive** (kisspeptin-10, `YNWNSFGLRF`
- a genuine Tyr1+aromatic-3 match that is not an opioid; a real limitation, not
a bug).

### Benchmark B - the hardest case: real peptides vs their own scrambles

Scrambles preserve composition and often the tyrosine, so this is where every
structure method failed worst.

| method | AUROC | p |
|---|---|---|
| pharmacophore R1 | 0.875 | <1e-5 |
| **pharmacophore R2, D30-compliant** | **0.893** | 3e-5 |
| pharmacophore R2, all 20 | 0.825 | 2e-5 |
| ADCP energy | 0.562 | 0.508 |
| Chai-1 iLIS | 0.412 | 0.389 |
| **Protenix ipTM (+MSA)** | **0.290** | **0.024** |

R2 produced **zero false positives** on the scrambles (0/20).

### Benchmark C - length controls

iLIS is a real improvement over ipTM here: **0.463 -> 0.741**, with perfect
separation (0/6 negatives above the positive) at 10 and 15 aa. It fails at
20 aa. Pharmacophore rules are trivially 1.000 on this set because the positives
are `YGGFM...` by construction, so it cannot discriminate the two approaches.

## 4. Why the structure-based methods fail: it is length, not chemistry

The single finding that explains five failures. Every structure-based score
correlates with peptide length far more strongly than with opioid identity:

| score | rho(score, length) |
|---|---|
| Chai-1 **iLIS** | **-0.856** (p=3e-9) |
| ADCP best energy | -0.718 (p=7e-11) |
| Protenix ipTM, length controls | -0.879 within negatives |
| Chai-1 ipTM | -0.332 |
| **pharmacophore R2** | **-0.207 (p=0.281, NS)** |

Controlling for length destroys the structure-based signal and leaves the
pharmacophore intact:

- iLIS: 0.614 -> **0.500** stratified; length-residualised **0.486**; logistic
  regression `class ~ iLIS + length` gives iLIS p=0.937.
- R2: 0.967 -> **0.971** stratified.

**After length residualisation, the two nociceptin decoys rank #1 and #2** - the
two peptides that must rank lowest rank highest. The scores are reading peptide
size and packing plausibility, not the opioid message.

This also predicts the contradiction seen earlier: on met-enkephalin vs five
non-opioid peptides, Protenix ranked it 1/6 and ADCP 5/6, because Protenix's
length bias rewards short and ADCP's rewards long, and met-enkephalin was the
shortest peptide in that set.

## 5. What this framework can and cannot do

**Can:** classify binder vs non-binder at high accuracy, for free, with no
receptor structure, no GPU and no sampling; give an interpretable reason; be
applied to any peptide library instantly.

**Cannot:**
- **Rank affinity.** It is binary. It cannot tell a 1 nM ligand from a 10 uM one,
  and nothing measured here can either.
- **Distinguish agonist from antagonist.** Casoxins and lactoferroxins carry the
  message and are antagonists.
- **Catch non-Tyr1 binders.** Nociceptin-family ligands and many synthetic
  opioids would be missed.
- **Avoid false positives that genuinely match.** Kisspeptin-10 is the example.

**The honest framing for the pilot:** the pharmacophore is the baseline that any
expensive method must beat, and no method tested beats it. That makes the
deliverable a negative technical-readiness report with a quantified free
baseline - which is what `instructions.md` asked for.

## 6. Recommended pipeline

1. **Tier 0 - pharmacophore filter (free).** R2 on sequence. This is the
   screen.
2. **Tier 1 - structure only for pose, never for ranking.** Co-folding with an
   MSA gives a usable pocket (Protenix +MSA reproduces the 8F7Q orthosteric site
   to **1.03 A**), and that is worth having for mechanism and for MD setup. Use
   `iLIS`, not ipTM, if a confidence number is needed, and report it alongside
   length. Never rank a library on it.
3. **Tier 2 - physics only on a validated pose, for a handful of compounds.**
   ABFE/FEP on 2-5 peptides with measured affinities. MM/GBSA is not adequate:
   14 kcal/mol spread between GAFF2 and ff19SB on the same pose.
4. **Gate everything on ground truth.** The blocking problem remains that
   affinities in this repo are unsourced IC50 in uM, n=6. Enrichment cannot be
   measured without labels. Sourcing Ki values from BindingDB / ChEMBL /
   Guide to Pharmacology is the prerequisite for any further method comparison.

## 7. Caveats

- **n is small.** 14 vs 15 on benchmark A. AUROC 0.967 has a wide interval.
- **19 of 29 sequences are `ASSERTED`**, not verified against UniProt - tracked
  in the `verification` column of the benchmark file. NPY, kisspeptin-10 and
  five already in this repo are confirmed. Verify before publication.
- **R2 was chosen among four rules on one set.** All four beat all
  structure-based methods, so the ordering is robust, but the 0.967 specifically
  should be re-measured on a held-out set.
- **Decoys are assumed non-binders** based on their being other receptors'
  cognate ligands. That is a stronger assumption than a measurement, though much
  stronger than the scramble assumption.
- **The pharmacophore encodes prior knowledge.** It is not a discovery method;
  it will not find a novel chemotype that lacks Tyr1.

---

## 8. Ground truth: what public data actually exists, and a trap it sets

Executed 2026-09-28 via `scripts/00c_fetch_affinities.py` against ChEMBL target
**CHEMBL233** (Homo sapiens, Mu-type opioid receptor). Every value carries its
ChEMBL document id, so each is traceable. Output: `docs/affinities_chembl.tsv`.

### Finding 1: the pilot's actual library has no curated MOR affinity data

9 of 24 queried peptides returned any activity. **Every milk-derived peptide is
absent** - beta-casomorphin-5, beta-casomorphin-7, alpha-lactorphin,
beta-lactorphin, casoxin A, casoxin C, soymorphin-5 - as are leu-enkephalin,
dynorphin A and beta-endorphin. So the library this pilot was built to screen
**cannot be validated against curated public data at all**; its ground truth has
to come from primary literature, peptide by peptide.

Two name-matching failures are recorded so they are not trusted later:
endomorphin-1 and endomorphin-2 both resolved to the SAME molecule
(`CHEMBL4177061`), so at most one is right; and "naloxone" resolved to
*naloxone methiodide*, a quaternary derivative. Name-based ChEMBL lookup is not
reliable for peptides and any curated table must be checked by structure or
sequence, not by name.

### Finding 2: the data that does exist is sound, and sets a noise floor

| peptide | pKi (median) | n | inter-study spread |
|---|---|---|---|
| morphine | 9.49 | 4 | 0.00 log (same value propagated) |
| endomorphin (ambiguous which) | 9.44 | 1 | - |
| met-enkephalin | 8.70 | 3 | **1.10 log (12.7x)** |
| alpha-neoendorphin | 8.03 | 1 | - |
| nociceptin | 6.36 | 3 | **1.05 log (11.2x)** |
| neuropeptide Y | 5.29 | 1 | - |

The ordering is biologically correct across four logs, and it independently
**confirms two of this benchmark's decoy assumptions by measurement rather than
assertion**: nociceptin is a weak MOR ligand (Ki 270-3030 nM, plus two censored
`>1000 nM` records) and NPY is essentially inactive (Ki 5100 nM).

It also establishes a **noise floor of roughly 1 log unit**: the same peptide's
Ki varies 11-13x between studies. No method can be expected to correlate with
affinity more tightly than the affinity data agrees with itself, and the repo's
existing single unsourced IC50 values were never adequate regardless of method.

### Finding 3: a worked example of the trap this whole exercise has been about

Five peptides have both a cited Ki and a score from the cognate benchmark:

| peptide | pKi | length | iLIS |
|---|---|---|---|
| endomorphin_1_free | 9.44 | 4 | 0.4585 |
| met_enkephalin | 8.70 | 5 | 0.3523 |
| alpha_neoendorphin | 8.03 | 10 | 0.2370 |
| nociceptin_1_17 | 6.36 | 17 | 0.2405 |
| NPY_1_20 | 5.29 | 20 | 0.0000 |

**rho(iLIS, pKi) = +0.900, p = 0.037.** Taken at face value that is a
significant, publishable-looking claim that iLIS predicts MOR affinity - and it
would be the first positive structure-based result in this pilot.

**It is an artefact. rho(length, pKi) = -1.000.** In this subset the potent
peptides are exactly the short ones: lengths in descending order of affinity are
4, 5, 10, 17, 20 - perfectly monotonic. A scoring function that reads nothing
but peptide length would score rho = 1.000 against measured affinity here. The
+0.900 therefore carries no information about whether iLIS sees binding.

This is the same confound that produced every earlier false signal, now
demonstrated against **real, cited affinity data** rather than class labels. It
is the clearest argument in this document for why the benchmark, not the method,
was the thing that needed fixing first.

### What an adequate affinity benchmark would require

1. **>= 15-20 peptides with measured Ki at MOR**, from primary literature, with
   assay conditions recorded (radioligand, cell/membrane source, species).
2. **Affinity decorrelated from length** - explicitly checked, not assumed. The
   design must include potent long peptides (beta-endorphin, dynorphin A) and
   weak short ones, or the test is void before it starts.
3. **Assay types kept separate.** Ki, IC50 and EC50 are not interchangeable;
   EC50 measures function, not binding. The repo's existing values are IC50 in
   uM and must not be pooled with nM Ki data.
4. **Censored values excluded**, not treated as points: five nociceptin records
   are `>1000 nM`, which is a bound, not a measurement.
5. **A power calculation.** At n=5 nothing is testable; at the ~1 log noise floor
   above, detecting a true rho of 0.6 at 80% power needs roughly n >= 20.

Until that exists, the defensible claim is the classification result in
sections 3-4: the pharmacophore separates MOR binders from other receptors'
peptide ligands at AUROC 0.967 and survives length control at 0.971, while every
structure-based score tested collapses to 0.500. Whether any of them can rank
*potency* remains **untested**, because no adequate dataset exists to test it.

---

## 9. The affinity test, properly powered: no score ranks potency

Section 8 said an adequate affinity benchmark did not exist. **That was wrong,
and the correction matters** - it came from searching ChEMBL by NAME, which
misses peptides whose `pref_name` is null. Enumerating all Ki records against
CHEMBL233 instead yields 142 peptides spanning six log units, of which 36 are
clean linear canonical sequences 4-20 aa (`scripts/00d_build_affinity_set.py`,
`data/reference/affinity_set.tsv`).

### The benchmark passes every check the earlier ones failed

| property | value |
|---|---|
| peptides | 36 |
| pKi range | 4.71 - 9.14 (**4.43 log units**) |
| **rho(length, pKi)** | **+0.021, p=0.902** |
| power for true rho=0.5 | **0.88** (0.68 at rho=0.4) |
| provenance | ChEMBL document ids per value |

**Affinity is decorrelated from length**, so the artefact that invalidated every
previous result cannot operate here. Compare the 5-peptide subset of section 8,
where rho(length, pKi) = -1.000 and iLIS scored a spurious rho = +0.900.

Amidation self-controls: 23/36 are C-terminally amidated and the model input is
the free acid, but same-sequence pairs differ by only 0.08-0.78 log
(`YGGFM` 7.92 amidated vs 8.70 free; `YGGFL` 7.94 vs 8.03), inside the ~1 log
inter-study noise floor, and pKi does not differ by amidation overall (p=0.576).

### Result: every metric fails

Chai-1, 180 tasks, 0 FATALs. n=36.

| score | rho vs pKi | p | rho vs length |
|---|---|---|---|
| iLIS | +0.196 | 0.251 | **-0.834** |
| cLIS | +0.191 | 0.264 | -0.814 |
| LIS | +0.161 | 0.347 | -0.855 |
| ipTM | +0.216 | 0.205 | -0.004 |
| aggregate_score | +0.206 | 0.227 | -0.098 |
| i-pLDDT peptide | +0.191 | 0.264 | -0.532 |
| interface PAE min | +0.247 | 0.147 | -0.655 |

Free-acid subset (n=13, score and label describe the same molecule): iLIS
rho = **+0.005**, ipTM rho = +0.093.

Two distinct failure modes are now separable, which earlier sets could not do:
**iLIS is length-driven** (rho = -0.834 even here), while **ipTM is not
length-confounded on this set at all** (rho = -0.004) and simply carries **no
affinity signal**. The problem is therefore not only the length artefact - the
metrics genuinely do not encode binding strength.

### The internal control settles it

The set contains a structure-activity series that was not designed in: the
`YGGFL...` peptides (pKi 7.39-8.71) against the `YGPFL...` peptides
(pKi 4.71-5.51). A single **Gly3 -> Pro** substitution costs up to **1600x** in
measured affinity at matched length and near-identical composition.

| | measured affinity | iLIS |
|---|---|---|
| separates YGGFL* from YGPFL* | **AUROC 1.000** | **AUROC 0.409** (p=0.591) |

iLIS does not merely miss the substitution - it ranks the **inactive** series
slightly higher. And the extremes make the point plainly: `YAWF` (pKi 5.50,
weak) ranks 2nd by iLIS, while dynorphin A(1-17) `YGGFLRRIRPKLKWDNQ` (pKi 8.03,
potent) ranks **last of 36**. The top six by iLIS are all 4-7 aa; the bottom six
are 10-17 aa.

### Conclusion

**No structure-based score tested in this pilot ranks MOR peptide affinity**, on
a benchmark with 4.4 logs of cited dynamic range, length decorrelated from the
label, and 88% power to detect rho=0.5. This is a powered negative result, not
an absence of evidence, and it is the strongest claim this project can make.

It also completes the framework recommendation in section 5 with positive
evidence rather than inference: use the **pharmacophore** for classification
(AUROC 0.967, length-stratified 0.971), use **structure for pose only** - a
1.03 A orthosteric pocket is worth having for mechanism and MD - and do not
expect any of it to prioritise a library by potency. Ranking potency at this
peptide size requires physics (ABFE/FEP) on a handful of compounds, or a
ligand-based model trained on the affinity data now assembled here.
