# Assumptions and unverified values

Everything in this repository that is assumed, unverified, or carried forward
on trust. Checked items record what was actually verified and when.

Last updated: 2026-09-22.

---

## 1. Verified against the live cluster (2026-09-22)

These were checked directly on ORCA and are **not** assumptions.

| Item | Verified value |
|---|---|
| Partitions | `short` (4 h), `normal` (1 day, default), `long` (7 days), `osg` (3 days) |
| GPU nodes | 25 nodes, 4 GPUs each; `orcaga01-06` L40S (24 GPUs), `orcaga07-25` A30 (76 GPUs) |
| GPU memory | **L40S 46068 MiB**, A30 24 GB. The difference is material for complex prediction. |
| Node spec | 64 cores, 570 GB RAM, 440 GB usable `/tmp` |
| Apptainer | **Module only, not on PATH by default**: `apptainer/1.4.1-gcc-13.4.0` |
| `module` command | **Compute nodes only.** It is *not* defined on the login node. |
| Python | `python/3.12.12-gcc-13.4.0`; system python3 on a compute node is 3.9 |
| CUDA | `cuda/12.9.0-none-none` (loaded by default via spack `StdEnv`) |
| Plotting | `intel-python/24.0.0` exists, as the prior work's convention assumed |
| Outbound network | Compute nodes reach pypi, github and huggingface (HTTP 200) |
| Workspace tools | `ws_list`, `ws_allocate`, `ws_extend`, `ws_find` all present in `/usr/local/bin` |
| Workspace | `Kuhfeld_temp` did not exist; created 2026-09-22, expires 2026-10-22, 5 extensions left |
| Home quota | 11.5 GB used of a 92 GB soft / 100 GB hard quota |
| `module` in sbatch | **Not defined in a non-login shell.** Every sbatch here uses `#!/bin/bash -l`. |

### 1.1 Parallel scratch was full — storage policy deviates from the brief

`df -h /scratch` reported **57 T used of 57 T, 6.9 GB free (100%)**. This is a
shared-filesystem condition across all ORCA users, not this project's
footprint; the prior run recorded the same mount at ~92% in July 2026.

The build brief specifies that run data lives on the workspace and that home
holds scripts and configs only. With 6.9 GB free, the workspace cannot hold
model weights plus outputs, and an array submitted against it would fail one
task at a time as tasks tried to write.

**Deviation, made deliberately and at the user's direction:**
`storage.root_mode` is set to **`home`**, with `home_subdir: bionemo_mor_run`.
Home had ~80 GB of headroom.

Mitigations, since home is slow, quota-limited and backed up nightly:
  * hot intermediates stay on node-local `/tmp` (440 GB) and only the final
    `.cif` and `.json` per task are copied back, so nightly backup churn is small;
  * `storage.min_free_gb` (25 GB) is checked before anything is submitted —
    the check that would have caught the scratch condition up front;
  * the workspace path is still resolved at runtime, never hardcoded, and the
    full `ws_list` / `ws_allocate` / `ws_extend` logic is implemented and
    tested in `src/setup_workspace.py`.

**Set `root_mode` back to `workspace` once scratch has capacity.** No code
change is needed; it is one config value.

---

## 2. Receptor structures

### 2.1 Reconciliation of the local and cluster copies — they are different entries

The brief asks whether the cluster copy differs from the local copy. **It does,
and not only in numbering: they are different PDB entries.**

| | Local (this repo uses these) | Cluster `preliminary_agonist_binding_data_v2/` |
|---|---|---|
| Active state | `5c1m_MOR.pdb` — **5C1M**, 295 res, auth 52–347 | `receptor_clean.pdb` — **8F7Q**, 287 res, auth 66–352 |
| Inactive state | `4dkl_MOR.pdb` — **4DKL**, 282 res, auth 65–352 | (none; the inactive arm never completed) |
| D3.32 | **Asp147** (canonical) | **Asp149** (+2 offset) |
| Pocket residues | 147 / 148 / 293 / 297 | 149 / 150 / 295 / 299 |

Verified by direct residue-identity check on all four files:
`5c1m_MOR.pdb` and `4dkl_MOR.pdb` both have 147=ASP, 148=TYR, 293=TRP,
297=HIS; `receptor_clean.pdb` has 147=**SER**, 148=**ILE**, 293=**VAL**,
297=**PRO** and 149=ASP, 150=TYR, 295=TRP, 299=HIS.

**Which this repo uses and why:** the local 5C1M/4DKL pair, as the brief
specifies. They are the canonical active/inactive pair, they are both
prepared in canonical numbering (so the literature's Asp147 is literally
residue 147 in the file), and the cluster has no prepared inactive structure
at all. 8F7Q is retained in `config/receptors.yaml` as an optional positive
control only.

**Nothing hardcodes 147.** Each receptor declares its own `d332` and
`binding_site_residues`, and `run_msa.verify_residues()` aborts the run if a
declared residue is not the expected type. Three numbering schemes are in play
across this project and they do not transfer — this has already cost the
project two silent errors (see `PRIOR_WORK.md`).

### 2.2 5C1M is sterically sealed to a peptide — scope of that finding

Prior work measured rigid 5C1M with BU72 removed as **sealed at 1.8, 2.6 and
3.0 Å probe radii** (open only to water at 1.4 Å), and concluded that no
docking program can place a pentapeptide in that pocket using rigid 5C1M.

**That finding scoped a rigid-receptor docking protocol.** A structure
predictor does not dock into a fixed receptor; it predicts the complex, and
can model an induced-fit pocket opening. The finding therefore does not
transfer directly — but it does not vanish either:

* **Assumption:** that OpenFold3/Boltz-2 will produce a physically sensible
  peptide-bound 5C1M-like active state. **Unverified.**
* A predicted complex is **not** evidence that the pocket opens. If predicted
  active-state complexes place the peptide in the pocket, that is a model
  output, not a measurement.
* `rigid_docking_competent: false` is recorded on the active receptor in
  `config/receptors.yaml` so the caveat travels with the file.

### 2.3 4DKL T4-lysozyme chimera

4DKL fuses T4 lysozyme into **the same chain A**, residues 1002–1161. A naive
"receptor = chain A" yields 442 residues of MOR+T4L, and an active-vs-inactive
comparison then partly measures T4L. The local file was prepared with
`prepare_receptor.py --max-resi 999`; **verified** — 282 residues, maximum
residue number 352, no T4L present.

---

## 3. Peptides

* Sequences, names and the agonist/antagonist class come **verbatim** from the
  prior run (`peptides.py`, `scrambles.tsv`). None were invented.
* Scrambles are **reused**, not regenerated, so the two stacks are comparable.
  Composition was verified to match the parent sequence for all 20 pairs.
* **IC50 values: only 6 of 20 real peptides have one**, taken from the prior
  analysis (`anganost.py`). The remaining 14 are blank. No value was invented.
  * **Unverified:** the primary literature source for those 6 values. They are
    recorded as `prior_analysis:anganost.py`, which is a provenance
    breadcrumb, not a citation. **Before these numbers appear in a proposal or
    manuscript, trace each to its original paper.**
* **Unverified:** the agonist/antagonist class labels carry no citation in the
  prior work either. Same caveat.

### 3.1 The five-residue minimum — what it actually is

The brief asks that the ≥5-residue filter be applied and stated. It is
implemented (`screen.min_peptide_length: 5`) and it excludes 8 of 40 peptides
(4 real, 4 scrambles): **YPFP** (β-casomorphin-4), **YPYY** (casoxin B),
**YGLF** (α-lactorphin), **YLLF** (β-lactorphin).

**Its provenance must not be overstated.** Prior work traced this threshold to
a single input filter in the earlier pipeline —
`if len(sequence) < 5: return False, "Sequence too short"` — an ESMFold API
guard. It is **not** a simulation finding, and the four 4-mers it excludes all
have published activity. Earlier project documents reported it as a simulation
result; that is a known error being corrected here, not repeated.

**Consequence for this demo, and it matters:** casoxin B is one of only six
antagonists. After filtering, **5 antagonists remain against 11 agonists.**
Any class-level conclusion is underpowered before a single GPU is used. This
is reported in the `--dry-run` output of `validate_against_known.py` so it
cannot be missed.

---

## 4. Models

### 4.0 Verified on ORCA, 2026-09-22

A Python 3.12.12 virtualenv was built on a compute node and the following were
checked directly:

| Item | Status |
|---|---|
| **Boltz-2** | **Verified installed: `boltz` 2.2.1**, on `torch` 2.14.0+cu130 |
| Boltz-2 CLI flags | **Verified present**: `--seed`, `--use_msa_server`, `--msa_server_url`, `--cache`, `--diffusion_samples`, `--recycling_steps`, `--sampling_steps`, `--output_format`, `--accelerator`, `--affinity_mw_correction`, `--diffusion_samples_affinity`, `--sampling_steps_affinity`, `--method`, `--override` |
| MMseqs2 client | **Verified present**: `boltz.data.msa.mmseqs2.run_mmseqs2` — stage 1 uses it, so no separate ColabFold install is required |
| **OpenFold3** | **Published on PyPI as `openfold3` (0.5.0).** Package exists; whether its weights are fetchable without gated access is **NOT verified**. |
| `bionemo-recipes` | **Verified ABSENT from PyPI.** "No matching distribution found." |
| `bionemo-framework` | **Verified ABSENT from PyPI.** Same. |
| ColabFold | Published on PyPI (1.6.3). Not installed; the Boltz client is used instead. |
| Outbound network | pypi / github / huggingface reachable from compute nodes (HTTP 200) |

**Consequence for the "BioNeMo framework" framing.** The BioNeMo *framework*
packages are not pip-installable; BioNeMo recipe code lives in its source
repository and on NGC. What this pipeline actually runs is the **open-source
model stack that BioNeMo packages** — Boltz-2 (MIT), OpenFold3, and an MMseqs2
MSA search — installed directly from open sources. This is consistent with the
brief's requirement to avoid NIM containers and AI-Enterprise-licensed
components, but the repository should not be described as *running on the
BioNeMo framework* when it installs none of it. Described accurately in the
README as the open-source model stack.

### 4.1 Still unverified

| Item | Status |
|---|---|
| OpenFold3 open weights fetchable | **Unverified.** |
| OpenFold2 weights available under an open licence | **Unverified.** Stage disabled by default. |
| OpenFold3 CLI interface | **Unverified.** `src/predict_complex.py` is written to a plausible interface and will report a missing/failed engine rather than a stack trace. |
| MSA endpoint exercised end to end | In progress at time of writing. |

### 4.2 Seeds

The brief asks that the three-seed replication carry forward *if the models
expose a seed*. **Boltz-2's `--seed` flag is verified present** ("Seed to use
for random number generator. Default is None (no seeding."), so
`917 / 424242 / 20260715` are carried forward unchanged.

**Unverified:** that the seed actually changes the output — i.e. that three
seeds produce three *different* predictions rather than three identical ones.
If seeding turns out to be a no-op, the three "replicates" are one measurement
repeated, the seed-to-seed SD collapses to zero, and the noise floor in
`rank.py` becomes meaningless. **This must be checked on the first few
completed tasks before the full array is trusted** — it is the single cheapest
way to invalidate the whole replicate structure.

### 4.3 Boltz-2 cannot score a peptide binder — MEASURED, not assumed

This was an open question; it is now answered. Running Boltz-2 2.2.1 on ORCA
with MOR (5C1M, 295 aa) as chain A and Met-enkephalin (YGGFM) as chain B, with
`properties: affinity: binder: B`:

```
ValueError: Chain B is not a ligand! Affinity is currently only supported for ligands.
Failed to process probe.yaml. Skipping.
```

**Two consequences, and the second is the dangerous one:**

1. Boltz-2's affinity head accepts **ligand chains only**. It will not score a
   peptide, whatever the peptide is. This is a hard interface restriction, not
   a distribution caveat.
2. When the request is rejected, Boltz **skips the entire record** and exits
   **0**. `processed/manifest.json` came back `{"records": []}` — no structure,
   no confidence, no error exit code. A naive array would have produced 192
   silent no-ops and looked like it succeeded.

**What this means for the proposed pipeline.** Stage 4 as specified in the
brief — "Boltz-2: complex structure plus binding affinity estimate; supplies
the ranking" — **does not work for peptide ligands.** The complex-structure
half works; the binding-affinity half is unavailable for this class of binder.

**How the pipeline responds:**

* `models.affinity.request_affinity` defaults to **`never`**. The affinity is
  not requested, so the record is not dropped and the complex is produced.
* The ranking signal is the **interface confidence (ipTM)**. This is a
  legitimate ranking signal and an entirely different quantity from a binding
  affinity. It is never presented as one:
  * `scores.csv` records `boltz2:iptm` in the `model` column;
  * `rank.py` derives the direction of merit per row (ipTM: higher is better)
    and **refuses to rank** a file mixing quantities with opposite directions;
  * `predict_affinity.py` now **aborts** a task that produced no confidence
    output rather than recording an empty row, so the silent-skip failure mode
    cannot recur.
* `affinity_source` is recorded per task as `iptm_not_affinity`.

**This must be stated plainly in any report built on this run.** The headline
correlation is between *published IC50* and *predicted interface confidence*,
not between IC50 and a predicted affinity. Calling it an affinity correlation
would misrepresent what was computed.

---

## 5. Ranking

`ranking.differential_margin: 0.5` is a **placeholder**, not a physical
constant. It cannot be set correctly before the stack's noise floor is
measured.

`rank.py` computes the median within-peptide seed-to-seed SD and **warns when
the margin sits below it**, flagging every affected peptide
`below_noise_floor`. This is a direct response to the prior stack's failure
mode, where every class contrast measured was an order of magnitude below its
own seed-to-seed noise (score SD ≈ 29) and no call was meaningful.

**Set this from the measured noise floor after the demo run.**

---

## 6. Deviations from the build brief, collected

1. **`root_mode: home` rather than workspace** — parallel scratch was 100%
   full. §1.1. At the user's direction.
2. **Apptainer is a module, not a PATH binary** — the brief's "matching the
   existing container convention on ORCA" holds, but every sbatch must
   `module load apptainer/1.4.1-gcc-13.4.0` first, and `module` does not exist
   on the login node.
3. **`src/_common.py` and `src/_engines.py`** are shared helpers not listed in
   the brief's layout. `slurm/make_tasklist.py` and `slurm/submit_all.sh`
   likewise. They exist to keep `src/` free of duplicated path and engine logic.
4. **A venv fallback is allowed** (`containers.allow_venv_fallback`) so the
   pipeline is testable before the `.sif` is built. Still open-source only.
