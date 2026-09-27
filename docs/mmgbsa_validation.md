# MM/GBSA does not discriminate either, and we know exactly why

The physics tier was meant to replace the co-folding confidence scores after
those failed ([`triage_validation.md`](triage_validation.md)). It fails too -
but unlike the ipTM result, the mechanism here is identified and quantitative.

**MM/GBSA's total is a proxy for ligand size.**

## The measurements

16 real peptides, 18 scrambled decoys, scored from minimised complexes.

### Discrimination: chance

| metric | real median | decoy median | AUROC | p |
|---|---|---|---|---|
| delta_total | -40.31 | **-43.14** | 0.483 | 0.88 |
| vdwaals | -56.47 | -58.62 | 0.479 | 0.85 |
| eel | -81.90 | -75.33 | 0.517 | 0.88 |
| egb | 104.97 | 108.87 | 0.497 | 0.99 |
| esurf | -8.26 | -8.12 | 0.493 | 0.96 |

The decoys score **more favourably** than the real peptides.

### Calibration: backwards

| metric | rho vs -log10(IC50) | p |
|---|---|---|
| delta_total | **-0.643** | 0.119 |
| vdwaals | -0.643 | 0.119 |

Negative rho means the score ranks these peptides **the wrong way round**: the
weakest binders get the best scores.

## The cause, measured

| term | vs ligand heavy-atom count |
|---|---|
| **vdwaals** | **Spearman rho = -0.888, p = 2.6e-12; Pearson r = -0.930** |
| delta_total | rho = -0.557, p = 6.3e-4; r = -0.752 |
| eel | rho = -0.019, p = 0.92 |

The van der Waals term is very nearly a **linear function of how many atoms the
ligand has**, and it dominates the total. Electrostatics, correctly, does not
care about size.

This is the textbook MM/GBSA pathology: with no configurational-entropy term,
every additional atom adds favourable contacts and nothing pays for the
flexibility it costs. Bigger always looks better.

It explains both failures at once:

| observation | explanation |
|---|---|
| anti-correlation with potency | the shorter casomorphins happen to be the more potent ones, so a size-driven score ranks them backwards |
| AUROC 0.483 against decoys | a scramble has **identical composition and length**, so identical size, so the same score |

Seen directly in the benchmark set - the score tracks atom count, not IC50:

| peptide | IC50 µM | ligand atoms | delta_total |
|---|---|---|---|
| met-enkephalin | 0.2 | 75 | -40.6 |
| beta-casomorphin-5 bov | 6.5 | 79 | -31.5 |
| beta-casomorphin-5 hum | 13.5 | 89 | -35.9 |
| beta-casomorphin-4 bov | 22.0 | 72 | -36.1 |
| beta-casomorphin-7 hum | 29.0 | **122** | **-52.4** |
| beta-casomorphin-7 bov | 57.0 | 110 | -48.3 |
| neocasomorphin-6 | 59.0 | 103 | -46.5 |

The two 7-mers - among the weakest binders in the set - get the two best scores.

## Size normalisation removes the artefact but finds no signal

| | AUROC vs decoys | rho vs IC50 |
|---|---|---|
| raw delta_total | 0.483 (p = 0.88) | **-0.643** |
| per ligand heavy atom | **0.601** (p = 0.33) | -0.036 |

Dividing by atom count does what it should: the spurious anti-correlation
collapses to zero. AUROC rises to 0.601, which is the most encouraging number
in the project so far - and it is **not significant** at 16 vs 18 (p = 0.33).
It is a lead, not a result.

## Both tiers now fail, and both fail on size

| tier | AUROC | direction of size bias |
|---|---|---|
| ipTM, sequence only | 0.376 | **penalises** small ligands (fewer interface pairs) |
| ipTM, 6DDF template | 0.415 | same |
| MM/GBSA raw | 0.483 | **rewards** large ligands (no entropy term) |
| MM/GBSA per-atom | 0.601 | corrected; no signal left |

Neither tier measures **specificity**. They measure size, in opposite
directions. Matched-permutation decoys are what makes this visible: they hold
composition and length fixed, so a score that survives them has to be
responding to sequence order.

## The common factor is the pose, not the scoring function

Both tiers score a peptide that **co-folding placed in the orthosteric pocket**,
and co-folding does that just as confidently for a scramble as for a real
peptide. Any scoring function applied afterwards is being asked to rank
structures that are all "bound", and the thing that would distinguish a
non-binder - that it does not adopt a favourable pose at all - has already been
removed.

**The bottleneck is pose generation, not scoring.** A better scoring function on
these poses cannot fix it.

This matters directly for Stage 6: **ABFE inherits the same defect.** It would
compute a careful free energy for a pocket-bound scramble, at roughly 36 lambda
windows per peptide. Worth knowing before spending that.

## What follows

1. **Any MM/GBSA ranking across peptides of different lengths must be
   size-normalised.** Raw totals are a length assay.
2. **Score poses that were not assumed to bind.** Docking into the experimental
   receptor lets a non-binder fail to find a favourable pose - the discriminating
   information that co-folding discards.
3. **Stage 6 stays scoped to method validation on DAMGO**, not screening. Its
   pose is experimental, so it does not inherit this problem; nothing else
   currently qualifies.
4. **Report this as a finding.** Two independent tiers, matched controls, a
   measured mechanism (r = -0.930), and a specific diagnosis. That is a
   stronger result than an unvalidated correlation would have been.

## Coverage

34 of 37 peptides scored. Three failed and are excluded, not worked around:

| peptide | failure |
|---|---|
| `scr__beta_casomorphin_7_hum` | charge assertion: summed to -1.0110, expected -1 |
| `alphas1_casein_exorphin_7` | antechamber failed |
| `scr__lactoferroxin_A_hum` | antechamber failed |

The charge case is a tolerance question - 0.011 e is physically negligible, but
the assertion that rejects it is the same one that caught the whole-unit DAMGO
charge error, so it is not being loosened without thought. The two antechamber
failures are unexplained and are the longest peptides involved, which suggests
sqm convergence rather than anything specific to them.
