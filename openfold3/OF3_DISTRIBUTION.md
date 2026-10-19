# OpenFold3: known agonists vs. diverse background

Same 8F7Q-templated protocol for both sets (3 seeds x 5 samples), so the groups are directly comparable.

| group | n | description |
|---|---:|---|
| known_agonist | 14 | published MOR agonists, N-terminal Tyr pharmacophore |
| background_tyr | 14 | background peptides containing a Tyr |
| background_noY | 36 | background peptides with NO Tyr -- cannot form the message motif |
| agonist_scramble | 14 | scrambles of the agonists |
| background_scramble | 50 | scrambles of the background |

**The no-Tyr group is the sharpest negative control available**: those peptides lack the opioid pharmacophore entirely. Their anchor distance is measured from residue 1 and is NOT comparable to the Tyr groups, so it is reported separately and never pooled.

## 1. Distributions

**sample_ranking_score** (higher is better)

| group | n | median | sd | min | max |
|---|---:|---:|---:|---:|---:|
| known_agonist | 14 | 0.717 | 0.021 | 0.670 | 0.748 |
| background_tyr | 14 | 0.608 | 0.102 | 0.296 | 0.700 |
| background_noY | 36 | 0.611 | 0.110 | 0.295 | 0.794 |
| agonist_scramble | 14 | 0.741 | 0.037 | 0.687 | 0.821 |
| background_scramble | 50 | 0.619 | 0.080 | 0.423 | 0.802 |

**ipTM** (higher is better)

| group | n | median | sd | min | max |
|---|---:|---:|---:|---:|---:|
| known_agonist | 14 | 0.674 | 0.027 | 0.614 | 0.711 |
| background_tyr | 14 | 0.526 | 0.131 | 0.136 | 0.652 |
| background_noY | 36 | 0.539 | 0.136 | 0.178 | 0.766 |
| agonist_scramble | 14 | 0.706 | 0.046 | 0.632 | 0.801 |
| background_scramble | 50 | 0.545 | 0.104 | 0.298 | 0.775 |

**pTM** (higher is better)

| group | n | median | sd | min | max |
|---|---:|---:|---:|---:|---:|
| known_agonist | 14 | 0.850 | 0.009 | 0.827 | 0.861 |
| background_tyr | 14 | 0.823 | 0.036 | 0.724 | 0.860 |
| background_noY | 36 | 0.825 | 0.035 | 0.728 | 0.868 |
| agonist_scramble | 14 | 0.852 | 0.009 | 0.832 | 0.863 |
| background_scramble | 50 | 0.827 | 0.028 | 0.738 | 0.867 |

**chain pTM (peptide)** (higher is better)

| group | n | median | sd | min | max |
|---|---:|---:|---:|---:|---:|
| known_agonist | 14 | 0.105 | 0.017 | 0.084 | 0.136 |
| background_tyr | 14 | 0.060 | 0.109 | 0.041 | 0.463 |
| background_noY | 36 | 0.073 | 0.055 | 0.028 | 0.322 |
| agonist_scramble | 14 | 0.105 | 0.019 | 0.084 | 0.141 |
| background_scramble | 50 | 0.065 | 0.046 | 0.029 | 0.235 |

**avg pLDDT** (higher is better)

| group | n | median | sd | min | max |
|---|---:|---:|---:|---:|---:|
| known_agonist | 14 | 81.130 | 0.726 | 79.242 | 81.990 |
| background_tyr | 14 | 78.426 | 2.935 | 70.957 | 81.762 |
| background_noY | 36 | 78.486 | 2.856 | 72.044 | 82.689 |
| agonist_scramble | 14 | 81.180 | 0.718 | 79.649 | 82.313 |
| background_scramble | 50 | 78.927 | 2.410 | 72.080 | 82.417 |

**gPDE** (lower is better)

| group | n | median | sd | min | max |
|---|---:|---:|---:|---:|---:|
| known_agonist | 14 | 0.506 | 0.022 | 0.483 | 0.565 |
| background_tyr | 14 | 0.626 | 0.146 | 0.491 | 0.974 |
| background_noY | 36 | 0.635 | 0.147 | 0.463 | 1.041 |
| agonist_scramble | 14 | 0.511 | 0.023 | 0.472 | 0.552 |
| background_scramble | 50 | 0.615 | 0.111 | 0.478 | 0.913 |

**anchor best (A)** (lower is better)

| group | n | median | sd | min | max |
|---|---:|---:|---:|---:|---:|
| known_agonist | 14 | 3.314 | 0.684 | 2.026 | 4.197 |
| background_tyr | 14 | 5.073 | 3.073 | 2.499 | 13.152 |
| background_noY | 36 | 5.068 | 6.505 | 1.399 | 26.900 |
| agonist_scramble | 14 | 2.605 | 3.070 | 1.671 | 11.921 |
| background_scramble | 50 | 5.061 | 6.142 | 1.510 | 28.308 |

**anchor reproducibility** (higher is better)

| group | n | median | sd | min | max |
|---|---:|---:|---:|---:|---:|
| known_agonist | 14 | 0.333 | 0.157 | 0.067 | 0.533 |
| background_tyr | 14 | 0.000 | 0.137 | 0.000 | 0.400 |
| background_noY | 36 | 0.000 | 0.332 | 0.000 | 1.000 |
| agonist_scramble | 14 | 0.600 | 0.390 | 0.000 | 1.000 |
| background_scramble | 50 | 0.000 | 0.323 | 0.000 | 1.000 |

**pocket contacts** (higher is better)

| group | n | median | sd | min | max |
|---|---:|---:|---:|---:|---:|
| known_agonist | 14 | 9.500 | 0.646 | 8.000 | 10.000 |
| background_tyr | 14 | 9.500 | 2.598 | 0.000 | 10.000 |
| background_noY | 36 | 9.000 | 1.496 | 2.000 | 10.000 |
| agonist_scramble | 14 | 10.000 | 0.514 | 9.000 | 10.000 |
| background_scramble | 50 | 9.000 | 0.694 | 8.000 | 10.000 |

**depth (A)** (lower is better)

| group | n | median | sd | min | max |
|---|---:|---:|---:|---:|---:|
| known_agonist | 14 | 7.356 | 0.967 | 6.124 | 9.202 |
| background_tyr | 14 | 10.587 | 4.687 | 7.635 | 25.609 |
| background_noY | 36 | 10.503 | 3.937 | 7.065 | 20.904 |
| agonist_scramble | 14 | 6.942 | 1.344 | 5.230 | 9.337 |
| background_scramble | 50 | 10.592 | 2.832 | 6.483 | 16.962 |

## 2. Can any quantity separate known agonists from background?

AUROC of known agonists over the background (Tyr-containing and no-Tyr pooled as `background_real`). 0.5 = no discrimination.

| quantity | AUROC | 95% CI | Mann-Whitney p |
|---|---:|---|---:|
| sample_ranking_score | **0.881** ** | 0.793 - 0.953 | 0.0000 |
| ipTM | **0.894** ** | 0.809 - 0.961 | 0.0000 |
| pTM | **0.831** ** | 0.727 - 0.921 | 0.0002 |
| chain pTM (peptide) | **0.820** ** | 0.716 - 0.909 | 0.0003 |
| avg pLDDT | **0.844** ** | 0.741 - 0.931 | 0.0001 |
| gPDE | **0.871** ** | 0.779 - 0.949 | 0.0000 |
| anchor best (A) | **0.711** ** | 0.586 - 0.827 | 0.0166 |
| anchor reproducibility | **0.790** ** | 0.675 - 0.891 | 0.0006 |
| pocket contacts | **0.594** | 0.440 - 0.741 | 0.2443 |
| depth (A) | **0.906** ** | 0.817 - 0.967 | 0.0000 |

`**` marks a CI excluding 0.5. Any quantity whose CI includes 0.5 cannot distinguish a known MOR agonist from an arbitrary milk peptide.

### 2a. Length-matched AUROC

Known agonists span **4-7 residues**; the background spans **5-29**. Section 4 shows most quantities track length at |rho| > 0.9, so the AUROCs above may be measuring peptide size rather than MOR affinity. Restricting the background to the agonists' own length range is the control that separates the two.

- background peptides within 4-7 residues: **9** of 50

| quantity | AUROC (all) | AUROC (length-matched) | 95% CI | p |
|---|---:|---:|---|---:|
| sample_ranking_score | 0.881 | **0.595** | 0.325 - 0.849 | 0.4688 |
| ipTM | 0.894 | **0.627** | 0.357 - 0.873 | 0.3289 |
| pTM | 0.831 | **0.365** | 0.119 - 0.635 | 0.2986 |
| chain pTM (peptide) | 0.820 | **0.659** | 0.405 - 0.889 | 0.2193 |
| avg pLDDT | 0.844 | **0.349** | 0.103 - 0.619 | 0.2439 |
| gPDE | 0.871 | **0.421** | 0.159 - 0.683 | 0.5495 |
| anchor best (A) | 0.711 | **0.556** | 0.270 - 0.817 | 0.6822 |
| anchor reproducibility | 0.790 | **0.627** | 0.353 - 0.881 | 0.3253 |
| pocket contacts | 0.594 | **0.560** | 0.357 - 0.766 | 0.6169 |
| depth (A) | 0.906 | **0.595** | 0.349 - 0.825 | 0.4688 |

A quantity whose length-matched CI includes 0.5 was discriminating on peptide size, not on MOR binding.

### 2b. Agonists vs. their own scrambles (length identical by construction)

| quantity | agonist median | scramble median | AUROC | p |
|---|---:|---:|---:|---:|
| sample_ranking_score | 0.717 | 0.741 | **0.321** | 0.1129 |
| ipTM | 0.674 | 0.706 | **0.321** | 0.1129 |
| pTM | 0.850 | 0.852 | **0.434** | 0.5657 |
| chain pTM (peptide) | 0.105 | 0.105 | **0.500** | 1.0000 |
| avg pLDDT | 81.130 | 81.180 | **0.469** | 0.8005 |
| gPDE | 0.506 | 0.511 | **0.505** | 0.9817 |
| anchor best (A) | 3.314 | 2.605 | **0.352** | 0.1904 |
| anchor reproducibility | 0.333 | 0.600 | **0.344** | 0.1659 |
| pocket contacts | 9.500 | 10.000 | **0.449** | 0.6185 |
| depth (A) | 7.356 | 6.942 | **0.388** | 0.3232 |

Scrambles are length-matched to their parents by construction, so this comparison is free of the length confound entirely. AUROC below 0.5 means the scramble scores BETTER than the real peptide.

## 3. Does the pharmacophore matter?

Peptides with no Tyr cannot form the Tyr->D3.32 message-motif contact. If they nonetheless reach the pocket as readily as the agonists, pose formation is not reporting on opioid pharmacology.

- **pocket contacts**: agonists median 9.500, no-Tyr background median 9.000, AUROC 0.629 (p = 0.1237)
- **depth (A)**: agonists median 7.356, no-Tyr background median 10.503, AUROC 0.897 (p = 0.0000)
- **ipTM**: agonists median 0.674, no-Tyr background median 0.539, AUROC 0.861 (p = 0.0001)
- **sample_ranking_score**: agonists median 0.717, no-Tyr background median 0.611, AUROC 0.843 (p = 0.0002)

## 4. Length confound

- rho(sample_ranking_score, length) = **-0.916** (p = 0.0000, n = 64)
- rho(ipTM, length) = **-0.925** (p = 0.0000, n = 64)
- rho(pTM, length) = **-0.891** (p = 0.0000, n = 64)
- rho(chain pTM (peptide), length) = **-0.407** (p = 0.0009, n = 64)
- rho(avg pLDDT, length) = **-0.920** (p = 0.0000, n = 64)
- rho(gPDE, length) = **0.937** (p = 0.0000, n = 64)
- rho(anchor best (A), length) = **0.342** (p = 0.0057, n = 64)
- rho(anchor reproducibility, length) = **-0.400** (p = 0.0011, n = 64)
- rho(depth (A), length) = **0.941** (p = 0.0000, n = 64)

