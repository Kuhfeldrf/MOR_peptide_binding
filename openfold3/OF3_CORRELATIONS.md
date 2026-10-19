# OpenFold3 quantities vs. IC50 — exhaustive search

OpenFold3 emits **no binding affinity**; it is a structure predictor. Every confidence quantity it does emit, plus every geometric observable of the predicted pose, is tested below.

- peptides with a published IC50: **9**
- quantities tested: **13**
- Bonferroni threshold for 13 tests: alpha = **0.0038**

## 1. All peptides with an IC50

n = 9 (met_enkephalin, beta_casomorphin_5_bov, beta_casomorphin_5_hum, beta_casomorphin_4_bov, beta_casomorphin_7_hum, alpha_lactorphin, beta_casomorphin_7_bov, neocasomorphin_6, beta_lactorphin)

| quantity | agrees if | Spearman rho | p | direction |
|---|---|---:|---:|---|
| sample_ranking_score (OF3 composite) | rho < 0 | **0.033** | 0.9322 | **wrong sign** |
| ipTM | rho < 0 | **0.000** | 1.0000 | none |
| pTM | rho < 0 | **-0.017** | 0.9661 | OK |
| chain pTM, peptide chain | rho < 0 | **-0.050** | 0.8984 | OK |
| avg pLDDT | rho < 0 | **-0.150** | 0.7001 | OK |
| gPDE | rho > 0 | **0.350** | 0.3558 | OK |
| disorder | rho > 0 | **0.017** | 0.9647 | OK |
| interface PAE (computed) | rho > 0 | **0.000** | 1.0000 | none |
| anchor distance, best sample | rho > 0 | **0.433** | 0.2440 | OK |
| anchor distance, median | rho > 0 | **0.350** | 0.3558 | OK |
| anchor reproducibility (frac <=5A) | rho < 0 | **-0.568** | 0.1103 | OK |
| pocket contacts, median | rho < 0 | **0.391** | 0.2977 | **wrong sign** |
| depth, median | rho > 0 | **0.333** | 0.3807 | OK |

- correct-direction quantities with raw p < 0.05: **0**
- surviving the Bonferroni threshold: **0**

## 2. Sensitivity analysis — excluding met_enkephalin

**Stated rationale:** met_enkephalin has an IC50 of 0.2 uM against a range of 6.5-160.0 uM for the remaining peptides — 32x below the next most potent. It is also the endogenous reference ligand rather than a milk-derived peptide.

**This is a sensitivity analysis, not the headline result.** Removing a point from an already small sample can only be read as a robustness check. Section 3 shows what happens when every other peptide is dropped in turn, which is the guard against reporting whichever exclusion happens to help.

## 2a. Excluding met_enkephalin

n = 8 (beta_casomorphin_5_bov, beta_casomorphin_5_hum, beta_casomorphin_4_bov, beta_casomorphin_7_hum, alpha_lactorphin, beta_casomorphin_7_bov, neocasomorphin_6, beta_lactorphin)

| quantity | agrees if | Spearman rho | p | direction |
|---|---|---:|---:|---|
| sample_ranking_score (OF3 composite) | rho < 0 | **-0.119** | 0.7789 | OK |
| ipTM | rho < 0 | **-0.167** | 0.6932 | OK |
| pTM | rho < 0 | **-0.452** | 0.2604 | OK |
| chain pTM, peptide chain | rho < 0 | **-0.095** | 0.8225 | OK |
| avg pLDDT | rho < 0 | **-0.643** | 0.0856 | OK |
| gPDE | rho > 0 | **0.786** | 0.0208 | OK |
| disorder | rho > 0 | **0.012** | 0.9768 | OK |
| interface PAE (computed) | rho > 0 | **0.286** | 0.4927 | OK |
| anchor distance, best sample | rho > 0 | **0.190** | 0.6514 | OK |
| anchor distance, median | rho > 0 | **0.405** | 0.3199 | OK |
| anchor reproducibility (frac <=5A) | rho < 0 | **-0.552** | 0.1557 | OK |
| pocket contacts, median | rho < 0 | **0.866** | 0.0054 | **wrong sign** |
| depth, median | rho > 0 | **0.048** | 0.9108 | OK |

- correct-direction quantities with raw p < 0.05: **1** — gPDE
- surviving the Bonferroni threshold: **0**

## 3. Full leave-one-out

Each peptide dropped in turn, for the quantities that showed any correct-direction signal. If one row stands out, that result depends on a single point.

**sample_ranking_score (OF3 composite)** (agrees if rho < 0)

| dropped | n | rho | p |
|---|---:|---:|---:|
| *(none)* | 9 | **0.033** | 0.9322 |
| met_enkephalin | 8 | -0.119 | 0.7789 |
| beta_casomorphin_5_bov | 8 | -0.119 | 0.7789 |
| beta_casomorphin_5_hum | 8 | 0.214 | 0.6103 |
| beta_casomorphin_4_bov | 8 | 0.214 | 0.6103 |
| beta_casomorphin_7_hum | 8 | 0.143 | 0.7358 |
| alpha_lactorphin | 8 | -0.024 | 0.9554 |
| beta_casomorphin_7_bov | 8 | 0.214 | 0.6103 |
| neocasomorphin_6 | 8 | -0.119 | 0.7789 |
| beta_lactorphin | 8 | -0.119 | 0.7789 |

**avg pLDDT** (agrees if rho < 0)

| dropped | n | rho | p |
|---|---:|---:|---:|
| *(none)* | 9 | **-0.150** | 0.7001 |
| met_enkephalin | 8 | -0.643 | 0.0856 |
| beta_casomorphin_5_bov | 8 | 0.024 | 0.9554 |
| beta_casomorphin_5_hum | 8 | 0.024 | 0.9554 |
| beta_casomorphin_4_bov | 8 | -0.071 | 0.8665 |
| beta_casomorphin_7_hum | 8 | -0.071 | 0.8665 |
| alpha_lactorphin | 8 | -0.143 | 0.7358 |
| beta_casomorphin_7_bov | 8 | -0.143 | 0.7358 |
| neocasomorphin_6 | 8 | -0.262 | 0.5309 |
| beta_lactorphin | 8 | -0.095 | 0.8225 |

**gPDE** (agrees if rho > 0)

| dropped | n | rho | p |
|---|---:|---:|---:|
| *(none)* | 9 | **0.350** | 0.3558 |
| met_enkephalin | 8 | 0.786 | 0.0208 |
| beta_casomorphin_5_bov | 8 | 0.262 | 0.5309 |
| beta_casomorphin_5_hum | 8 | 0.262 | 0.5309 |
| beta_casomorphin_4_bov | 8 | 0.238 | 0.5702 |
| beta_casomorphin_7_hum | 8 | 0.238 | 0.5702 |
| alpha_lactorphin | 8 | 0.357 | 0.3851 |
| beta_casomorphin_7_bov | 8 | 0.238 | 0.5702 |
| neocasomorphin_6 | 8 | 0.381 | 0.3518 |
| beta_lactorphin | 8 | 0.381 | 0.3518 |

**anchor reproducibility (frac <=5A)** (agrees if rho < 0)

| dropped | n | rho | p |
|---|---:|---:|---:|
| *(none)* | 9 | **-0.568** | 0.1103 |
| met_enkephalin | 8 | -0.552 | 0.1557 |
| beta_casomorphin_5_bov | 8 | -0.552 | 0.1557 |
| beta_casomorphin_5_hum | 8 | -0.498 | 0.2091 |
| beta_casomorphin_4_bov | 8 | -0.503 | 0.2036 |
| beta_casomorphin_7_hum | 8 | -0.520 | 0.1865 |
| alpha_lactorphin | 8 | -0.520 | 0.1865 |
| beta_casomorphin_7_bov | 8 | -0.536 | 0.1706 |
| neocasomorphin_6 | 8 | -0.536 | 0.1706 |
| beta_lactorphin | 8 | -0.847 | 0.0080 |

**pocket contacts, median** (agrees if rho < 0)

| dropped | n | rho | p |
|---|---:|---:|---:|
| *(none)* | 9 | **0.391** | 0.2977 |
| met_enkephalin | 8 | 0.866 | 0.0054 |
| beta_casomorphin_5_bov | 8 | 0.282 | 0.4991 |
| beta_casomorphin_5_hum | 8 | 0.391 | 0.3379 |
| beta_casomorphin_4_bov | 8 | 0.391 | 0.3379 |
| beta_casomorphin_7_hum | 8 | 0.391 | 0.3379 |
| alpha_lactorphin | 8 | 0.391 | 0.3379 |
| beta_casomorphin_7_bov | 8 | 0.391 | 0.3379 |
| neocasomorphin_6 | 8 | 0.206 | 0.6242 |
| beta_lactorphin | 8 | 0.206 | 0.6242 |

