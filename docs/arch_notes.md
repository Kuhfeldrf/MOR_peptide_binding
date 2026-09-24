# Architecture notes

Per-tool build status, reusable as a porting risk assessment.

## Headline

`uname -m` returns **`x86_64`** (AMD Zen4), recorded 2026-09-24.

The build instructions anticipated `aarch64`, where Python wheel availability
is patchy. **That risk does not apply on this cluster.** Every tool in the
stack is expected to install from a prebuilt manylinux wheel or a conda-forge
binary; no source builds are anticipated.

This table is retained because it remains the correct artefact for assessing a
future port to aarch64 hardware, but its current contents describe an
x86_64 build and should not be read as evidence about aarch64.

## Per-tool status

Populated in step 3 as the environment is built.

| Tool | Arch | Wheel or source | Patches required | Time spent |
|------|------|-----------------|------------------|------------|
| _pending_ | | | | |

## Stage 4 engine choice: GROMACS

Decided 2026-09-24. The instructions permit either GROMACS or OpenMM and
require the reason to be recorded here.

### The deciding criterion is the target machine, not this one

ORCA is a **staging environment**. The intended target is a **new NVIDIA system
at OSU**. If that machine is Grace-based (GH200 / GB200) its CPUs are
**aarch64**, which is why the build instructions open with a `uname -m` check
that looks irrelevant on ORCA's x86_64 login node.

**The OSU system's architecture has not been confirmed.** GROMACS was chosen
partly *because* of that uncertainty: it is the option that degrades gracefully
if the target turns out to be ARM.

### Reasons

1. **Architecture-agnostic by construction.** GROMACS is C++/CMake and builds
   from source on any architecture. OpenMM's CUDA support arrives through the
   conda-forge binary ecosystem, whose aarch64 + CUDA coverage could not be
   verified. A source-buildable package does not have that failure mode.
2. **NVIDIA maintains optimised ARM/GH200 GROMACS builds via NGC containers**,
   and the prior work on this account already uses NVIDIA BioNeMo containers,
   so this fits the existing stack.
3. **The source build is itself a deliverable.** The instructions require build
   flags in `provenance.md` and per-tool build status here. Building GROMACS
   against CUDA on ORCA produces exactly that record and rehearses the OSU port
   rather than deferring it.
4. **Faster per ns for membrane systems** of this size, which matters for the
   scaling estimates in the README.

### Correction to an earlier assessment

An earlier draft of this analysis claimed GROMACS would require hand-assembled
lambda schedules for Stage 6, and counted that against it. **That was wrong and
is recorded here rather than quietly removed.** GROMACS free-energy
perturbation is a first-class feature: `couple-moltype`, `init-lambda-state`,
and Boresch restraints via `intermolecular_interactions`, emitting `dhdl.xvg`
that `alchemlyb` / `pymbar` consume directly for MBAR. Once corrected, OpenMM's
principal advantage largely disappears.

### Costs accepted

- A CUDA 12.9 source build on Zen4, which may take a few hours of wall time.
- **Topology conversion at the Stage 3 boundary.** `packmol-memgen` emits Amber
  topologies; GROMACS needs them converted via ParmEd. A routine but real extra
  step, and a genuine place for silent error. **Stage 3 must verify atom count
  and total system charge across the conversion and fail loudly on mismatch.**
  (This risk was briefly removed by a CHARMM36m decision on 2026-09-24 and
  returned when that decision was reversed - see `docs/forcefield_decision.md`.)

### To revisit

If the OSU system is confirmed x86_64 (e.g. a DGX platform), OpenMM becomes the
better choice - no source build, and `openmmtools` supplies the alchemical
machinery directly. Re-open this decision if that is established.

## Dependency conflict: numpy, Chai-1, and the conda-forge stack

**Found 2026-09-24 on the first environment build. This is a porting finding,
not a local accident, and it should be expected to recur on the OSU system.**

A single environment containing both `chai_lab` and the conda-forge scientific
stack **does not work**. `pip install chai_lab` pulls `torch 2.6.0+cu124`,
which downgrades numpy to **1.26.4**. Every conda-forge package in the
environment is built against numpy 2.x, so they then fail at import with:

```
AttributeError: module 'numpy' has no attribute 'long'
```

Observed breakage: `MDAnalysis` (Stages 4, 5), `pymbar` (Stage 6 MBAR),
`packmol-memgen` (Stage 3), plus unsatisfied numpy constraints reported for
`scipy`, `mdtraj`, `jax` and `snakemake-executor-plugin-slurm`. SLURM reported
the build job as `COMPLETED`; the environment was nonetheless unusable. **Job
exit status was not a reliable signal here.**

`pdb2pqr` was additionally missing: `packmol-memgen` requires it but it is not
pulled in transitively.

### Resolution: two environments

| Environment | File | Used by | numpy |
|-------------|------|---------|-------|
| `mor-pilot` | `environment.yml` | Stages 0, 1, 2b, 3-7 | >= 2 |
| `mor-chai`  | `envs/chai.yml`    | Stage 2 only | < 2 |

Snakemake per-rule `conda:` directives select between them, with a shared
`conda-prefix` so the torch environment is built once. The split is honest
rather than tidy: Stage 2 genuinely is an isolated GPU job with an incompatible
dependency set, and pinning numpy < 2 globally would mean fighting conda-forge
across the entire stack.

**Consequence for reproducibility:** there are two environment files to pin,
not one. Both must be recorded in `provenance.md`.

## Per-tool build status

Recorded on x86_64. See the headline note above: this is not aarch64 evidence.

| Tool | Arch | Wheel or source | Patches | Time | Notes |
|------|------|-----------------|---------|------|-------|
| GROMACS 2025.3 | x86_64 | **spack source build** | none | ~5 min / 16 cores | CUDA 12.9, cuda_arch 80,89. No binary cache; compiled clean first attempt. |
| AmberTools (packmol-memgen) | x86_64 | conda-forge binary | none | - | Needs `pdb2pqr`, which is **not** pulled in transitively. |
| chai_lab 0.6.1 | x86_64 | pip wheel | none | - | Drags `torch 2.6.0+cu124`; forces numpy<2. Isolated env. |
| pymbar / alchemlyb / MDAnalysis | x86_64 | conda-forge binary | none | - | Broke under numpy<2; fixed by the env split. |

### SIMD flag note (performance, not correctness)

Spack targeted `linux-zen4` but emitted `-march=skylake-avx512`. GROMACS
selected `AVX_512`, which Zen4 does support natively, so the build is correct
and runs. It is, however, **not tuned for Zen4** - `znver4` would be the
matching target, and GROMACS on Zen4 is sometimes faster with `AVX2_256` than
`AVX_512` depending on kernel.

This is left as-is for the pilot: it is a performance question, not a
correctness one, and the pilot is not a benchmark. **It should be measured
before any scaling decision**, since MD throughput drives the GPU-hour estimate
in the README scaling section.

### Build failure worth recording

The first GROMACS job was reported `FAILED` by SLURM although the build had
**succeeded**. The failure was in the script's own recording step: `set -u`
tripped over `GMXRC`, which references unbound shell variables (`shell`,
`GMXLDLIB`). Source `GMXRC` with `set +u`.

Together with the numpy incident above, where a broken environment was reported
`COMPLETED`, the lesson is the same in both directions: **on this cluster, job
exit status alone establishes nothing.** Verify the artefact.
