# OpenFold3 vs. published IC50 — MOR agonists, 8F7Q templated

- models measured: **420** across 28 peptides
- receptor frame verified per model (offset auto-detected: [65])
- anchor = peptide anchor-Tyr alpha-amino N to Asp149 OD, gate 5.0 A

## 1. Per-peptide results (real agonists)

| peptide | seq | len | IC50 (uM) | anchor best (A) | anchor med | pocket med | ipTM med | anchor formed |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| met_enkephalin | YGGFM | 5 | 0.2 | 2.03 | 5.49 | 10.0 | 0.667 | 7/15 |
| beta_casomorphin_5_bov | YPFPG | 5 | 6.5 | 3.26 | 5.46 | 8.0 | 0.657 | 7/15 |
| beta_casomorphin_5_hum | YPFVE | 5 | 14.0 | 2.18 | 4.59 | 9.0 | 0.705 | 8/15 |
| beta_casomorphin_4_bov | YPFP | 4 | 22.0 | 3.37 | 5.18 | 9.0 | 0.683 | 7/15 |
| beta_casomorphin_7_hum | YPFVEPI | 7 | 25.0 | 2.62 | 5.66 | 9.0 | 0.647 | 5/15 |
| alpha_lactorphin | YGLF | 4 | 50.0 | 3.06 | 6.94 | 9.0 | 0.675 | 5/15 |
| beta_casomorphin_7_bov | YPFPGPI | 7 | 57.0 | 3.99 | 8.23 | 9.0 | 0.614 | 1/15 |
| neocasomorphin_6 | YPVEPF | 6 | 59.0 | 3.52 | 12.10 | 10.0 | 0.673 | 2/15 |
| beta_lactorphin | YLLF | 4 | 160.0 | 2.54 | 5.01 | 10.0 | 0.680 | 7/15 |
| alphas1_casein_exorphin_7 | RYLGYLE | 7 | - | 3.71 | 6.31 | 10.0 | 0.659 | 4/15 |
| alphas1_casein_exorphin_6 | RYLGYL | 6 | - | 4.20 | 7.14 | 10.0 | 0.708 | 1/15 |
| soymorphin_5 | YPFVV | 5 | - | 2.30 | 6.05 | 10.0 | 0.711 | 6/15 |
| soymorphin_6 | YPFVVN | 6 | - | 3.46 | 7.40 | 9.0 | 0.700 | 5/15 |
| soymorphin_7 | YPFVVNA | 7 | - | 3.45 | 6.53 | 10.0 | 0.661 | 3/15 |

## 2. Correlation with IC50 (n = 9)

| quantity | agreement direction | Spearman rho | p | partial rho (length controlled) | p |
|---|---|---:|---:|---:|---:|
| anchor distance, best | positive | **0.433** | 0.2440 | 0.447 | 0.2276 |
| anchor distance, median | positive | **0.350** | 0.3558 | 0.404 | 0.2809 |
| pocket contacts, median | negative | **0.391** | 0.2977 | 0.394 | 0.2942 |
| ipTM, median | negative | **0.000** | 1.0000 | 0.019 | 0.9616 |
| depth, median | positive | **0.333** | 0.3807 | 0.539 | 0.1339 |
| anchor reproducibility (frac of samples <=5A) | negative | **-0.568** | 0.1103 | -0.665 | 0.0508 |

*Agreement direction* is what a working predictor would show: a shorter anchor distance, more pocket contacts or higher confidence should accompany a LOWER IC50.

- Spearman rho(length, log IC50) = **0.017** (p = 0.9647) — the confound, measured on this set.
- Spearman rho(ipTM, length) = **-0.736** (p = 0.0236)

## 3. Real vs. scrambled controls

Anchor distance, best sample. Positive = the real peptide anchors closer than its composition-matched scramble.

- pairs: **14**
- real anchors closer in **5 of 14**
- mean difference: **0.62 A**
- Wilcoxon signed-rank p = **0.6257**

| peptide | real anchor | scramble anchor | diff |
|---|---:|---:|---:|
| met_enkephalin | 2.03 | 2.11 | 0.08 |
| beta_casomorphin_5_hum | 2.18 | 2.59 | 0.41 |
| soymorphin_5 | 2.30 | 2.00 | -0.29 |
| beta_lactorphin | 2.54 | 1.91 | -0.63 |
| beta_casomorphin_7_hum | 2.62 | 11.92 | 9.30 |
| alpha_lactorphin | 3.06 | 2.68 | -0.38 |
| beta_casomorphin_5_bov | 3.26 | 6.44 | 3.18 |
| beta_casomorphin_4_bov | 3.37 | 2.62 | -0.75 |
| soymorphin_7 | 3.45 | 1.67 | -1.78 |
| soymorphin_6 | 3.46 | 8.66 | 5.20 |
| neocasomorphin_6 | 3.52 | 2.24 | -1.28 |
| alphas1_casein_exorphin_7 | 3.71 | 2.01 | -1.70 |
| beta_casomorphin_7_bov | 3.99 | 2.93 | -1.06 |
| alphas1_casein_exorphin_6 | 4.20 | 2.65 | -1.55 |

