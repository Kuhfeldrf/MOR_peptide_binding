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
  topologies; GROMACS needs them converted via ParmEd or acpype. This is a
  routine but real extra step, and a genuine place for silent error. Stage 3
  must verify atom counts and total charge across the conversion.

### To revisit

If the OSU system is confirmed x86_64 (e.g. a DGX platform), OpenMM becomes the
better choice - no source build, and `openmmtools` supplies the alchemical
machinery directly. Re-open this decision if that is established.

## Per-tool build status

Populated in step 3 as the environment is built.

| Tool | Arch | Wheel or source | Patches required | Time spent |
|------|------|-----------------|------------------|------------|
| _pending_ | | | | |
