# The anchor set confirms it: the Protenix score is a length detector

19 anchor peptides, 5 seeds, 5 samples each — 475 predicted complexes on a
correctly folded receptor with a real MSA — and the score still reads peptide
length rather than affinity. This is the same conclusion the scrambled-decoy
and positive-control experiments reached, now with **cited affinities from four
primary documents** instead of class labels.

## What was different this time

Every previous arm B attempt could be dismissed on technical grounds. This one
cannot:

| Prior objection | Status here |
|---|---|
| Receptor never folded (arm B v1: 16.2 Å) | Gate passed at **1.74 Å** CA RMSD vs 8F7Q chain R; orthosteric pocket **0.70 Å** |
| Ran MSA-free, unlike Chai-1 | Receptor MSA with **N_msa = 8862** (4801 paired + 4178 unpaired) |
| Only binary labels, no real affinities | Cited Ki / KD / IC50 for 19 of 20 anchors |
| Affinities pooled across incomparable assays | **No cross-block statistic is computed** |

So the negative result is not a failure of the setup. The machinery worked, and
the score still does not measure binding.

## The headline number

```
score vs peptide LENGTH, all 19 peptides:   rho = -0.878   p < 0.001
```

Within block B, where one document and one assay cover eight peptides:

```
score vs LENGTH:    rho = -0.927   p = 0.001   n = 8
score vs AFFINITY:  rho = +0.571   p = 0.139   n = 8
```

Length explains block B's scores; affinity does not reach significance. The two
Koch length ladders are perfectly monotone in length — bovine
`rho = -1.000`, human `rho = -1.000` — in a block where length and affinity are
themselves uncorrelated (`rho = -0.293`, p = 0.48). That is the cleanest
possible demonstration: where length and affinity come apart, the score follows
length.

## The comparison that needs no statistics

| peptide | length | measured affinity | score |
|---|---|---|---|
| `endomorphin_1` (YPWF) | 4 | **Ki 1.86 nM** | 0.922 |
| `kappa_casein_f61_63` (YPY) | 3 | **> 1000 µM** | 0.921 |

Five-plus orders of magnitude apart in affinity; **0.001 apart in score.** The
near-inactive 3-mer outscores β-endorphin (Ki 77.6 nM, score 0.853) and
dynorphin A 1-13 (Ki 955 nM, score 0.821). A score that ranks a non-binder
above two real agonists cannot prioritise a synthesis list.

These two sit in different blocks (A and C), so no correlation is computed
across them — but the point does not rest on a correlation.

## Aim 2's pair is not resolved

The casomorphin-7 pair is the design's cleanest affinity test: **same length,
same block, same assay, same document**, 3× apart in KD. Length is matched, so
the artefact is held constant.

| peptide | KD | score | seed spread |
|---|---|---|---|
| `bov_casomorphin_7` (YPFPGPI) | 2 µM | 0.892 | 0.024 |
| `hum_casomorphin_7` (YPFVEPI) | 6 µM | 0.878 | 0.010 |

The gap is **+0.014** — ordered correctly, but **smaller than the run-to-run
spread of the tighter binder (0.024)**. The ordering is not reproducible enough
to count. Aim 2 cannot rest on this score.

## Per-block detail

| block | source / assay | n scored | score vs affinity | score vs length |
|---|---|---|---|---|
| A | Wang 2023 Cell Table S2, pKi | 5 (+1 N.A.) | +0.700 (p = 0.19) | −0.564 (p = 0.32) |
| B | Koch 1985 Table 2, KD | 8 | +0.571 (p = 0.14) | **−0.927 (p = 0.001)** |
| C | Chiba 1989 Table 1, IC50 | 4 | −1.000 | −0.800 |
| D | Tani 1990 Table 1, IC50 | 1 | — | — |

Three cautions on this table, all of which cut against over-reading it:

- **Block A cannot settle anything.** Its apparent affinity correlation
  (+0.700) is not significant at n = 5, and within block A length and affinity
  are themselves confounded (`rho = -0.872`, p = 0.054). The block is too small
  and too collinear to separate the two hypotheses.
- **Block C's perfect inversion is fragile.** `rho = -1.000` on n = 4 is one
  rank swap away from collapsing, and one of its four values
  (`kappa_casein_f61_63`, > 1000 µM) is a **censored lower bound**, not a
  measurement. Read it as "certainly no positive signal," not as "perfectly
  anti-correlated."
- **Nociceptin FQ ranks last in block A, and that proves nothing.** It is a
  measured non-binder at µ, so ranking it 6/6 looks like a success — but at 17
  residues it is the second-longest peptide in the block, which is exactly what
  the length artefact predicts. The right conclusion is that this test cannot
  distinguish the two explanations, not that the score got it right.

## Reproducibility caveat on the gate

The gate ran twice on the same node (orcaga07) with the same seed 0. Four of
five samples reproduced to within 0.03 Å; **sample 3 moved from 2.59 Å to
1.74 Å.** The headline 1.74 Å is therefore the least reproducible of the five
numbers, because "best of 5" selects precisely the outlier.

The **verdict** is robust — both runs pass the 4.0 Å threshold on the corrected
metric, and the pocket is sub-ångström in both. The **number** should be quoted
as ~2.6 Å best-of-5 typical, with 1.74 Å as one observed draw, not as the
method's accuracy.

## A measurement bug was found and fixed en route

The gate initially FAILED at 4.47 Å. The cause was in the gate, not the model:
its RMSD paired CA atoms **by list index**, and Protenix renumbers its output
1–281 while 8F7Q chain R runs 66–352. Residue identity under the pairing the
code used was **11/281 (3.9%)**; at the correct one-residue offset it is
**276/281 (98.2%)**. Corrected, the same five structures are 2.59–2.81 Å.

The threshold was **not** changed. The calculation now lives in
`scripts/gate_rmsd.py`, which establishes the residue correspondence first and
**refuses to print an RMSD** unless some offset reaches 90% identity — a
misalignment must never again be reportable as a folding failure. It is 90% and
not 100% because the 6DDF-derived receptor and 8F7Q chain R genuinely differ at
5 of 281 positions (`V2→I68, N73→T139, F92→W158, V123→I189, I242→V308`).

**Consequence for the arm A comparison:** arm A's 2.72 Å baseline was computed
with the same index-paired code and is therefore also suspect. Until both are
recomputed with `gate_rmsd.py`, no arm A vs arm B structural comparison should
be quoted — including in `docs/pose_generation_comparison.md`.

## What this means for the grant

This strengthens rather than weakens the case built in
`docs/validated_screening_framework.md`. The pharmacophore filter remains the
only method in this project shown to discriminate, and the reason is now
documented against cited affinities from four independent papers rather than
against assumed labels: **co-folding confidence scores rank peptides by length,
and length is not affinity.** The hero run is justified by what cheap scoring
provably cannot do, and that claim now rests on a correctly folded receptor
with a real MSA, which removes the last technical escape route.

## Reproducing

```bash
sbatch jobs/anchor_pqr_env.sbatch      # pdb2pqr 3.6.1 / propka 3.5.1 (pinned)
sbatch jobs/anchor_stage0.sbatch       # receptor from 6DDF
sbatch jobs/anchor_ptx_msa.sbatch      # ONE receptor MSA, shared
sbatch jobs/anchor_ptx_gate.sbatch     # hard gate, <= 4.0 A, before any GPU spend
sbatch jobs/anchor_ptx_array.sbatch    # 19 peptides x 5 seeds = 95 tasks

python scripts/02e_aggregate.py --cofold results/02_cofold_ptx_anchor \
    --out docs/anchor_ptx_scores.tsv
python scripts/08b_anchor_analysis.py --scores docs/anchor_ptx_scores.tsv \
    --anchors data/reference/anchor_set.tsv \
    --library data/templates/anchor20.tsv --out docs/anchor_ptx_result.txt
```

Jobs: env 242493, Stage 0 242494, MSA 242495, gate 242507 (242496 failed on the
metric bug), array 242508 (95/95 COMPLETED). Raw per-model scores:
`docs/anchor_ptx_scores.tsv` (475 rows). Full analysis output:
`docs/anchor_ptx_result.txt`.
