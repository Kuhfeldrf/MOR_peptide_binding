# Provenance

Every tool, version, module load, build flag, and date obtained.

## Cluster

| Item | Value | Date recorded |
|------|-------|---------------|
| Cluster | ORCA, Portland State University ARC (`login.orca.pdx.edu`) | 2026-09-24 |
| Login node OS | Rocky Linux 9, kernel 5.14.0-687.29.1.el9_8 | 2026-09-24 |
| Architecture | `x86_64` (AMD Zen4) | 2026-09-24 |
| Scheduler | SLURM | 2026-09-24 |
| Scratch root | `/scratch/kuhfeldr-Kuhfeld_temp` (admin-provisioned) | 2026-09-24 |

**Note on institution:** the build instructions describe ORCA as an Oregon
State University cluster. The host reached is `login.orca.pdx.edu` with
`arcnfs.rc.pdx.edu` filesystems, i.e. Portland State University. The
discrepancy is recorded rather than resolved.

## Scratch

`/scratch` is not writable at the top level; per-user project directories
follow a `<user>-<project>` convention and are provisioned by administrators.
`/scratch/kuhfeldr-Kuhfeld_temp` already existed (created 2026-09-22, empty)
and is used as `PILOT_ROOT`. Filesystem: 57 T total, 8.4 T free at 86% use.
Home directory quota is 92 G (31.7 G used) and is too small for the pilot.

## SLURM partitions

| Partition | GRES | Walltime | Nodes |
|-----------|------|----------|-------|
| `short`   | `gpu:l40s:4` (orcaga01-06), `gpu:a30:4` (orcaga07-19) | 4:00:00 | 19 |
| `normal`* | `gpu:l40s:4` (orcaga01-05), `gpu:a30:4` (orcaga10-25) | 1-00:00:00 | 21 |
| `long`    | `gpu:a30:4` (orcaga11-19) | 7-00:00:00 | 9 |
| `osg`     | `gpu:a30:4` (orcaga20-25) | 3-00:00:00 | 6 |

Compute nodes: 64 cores, ~570 GB RAM. No A100/H100 present. L40S (Ada) and
A30 (Ampere) both support bfloat16, satisfying the Chai-1 requirement.

## Modules available

Lmod must be initialised explicitly in non-interactive shells and job scripts:

```bash
source /etc/profile.d/z00_lmod.sh
```

| Module | Version |
|--------|---------|
| `apptainer` | 1.4.1-gcc-13.4.0 |
| `cuda` | 11.8.0, 12.9.0 |
| `python` | 3.12.12, 3.14.0 |
| `gcc` | 13.4.0 |
| `cmake` | 3.31.9 |
| `openmpi` | 4.1.8, 5.0.8 |
| `uv` | 0.9.29 |
| `spack` | v1.0.2, v1.1.1 |

## Not available as modules

- conda / mamba / micromamba - Miniforge installed into `PILOT_ROOT` instead.
- GROMACS - absent; contributed to the OpenMM choice for Stage 4.
- AmberTools - absent; obtained from conda-forge for `packmol-memgen`.

## Tool versions

| Tool | Version | Source | Date |
|------|---------|--------|------|
| GROMACS | 2025.3-spack | spack source build, CUDA | 2026-09-24 |
| CUDA (build) | 12.9.0 | spack (`cuda@12.9.0`) | 2026-09-24 |
| gcc | 13.4.0 | spack | 2026-09-24 |
| Miniforge / conda | 26.7.2 | installer into PILOT_ROOT | 2026-09-24 |
| mamba | 2.9.0 | Miniforge | 2026-09-24 |
| chai_lab | 0.6.1 | pip (`mor-chai` env) | 2026-09-24 |
| torch | 2.6.0+cu124 | pip, via chai_lab | 2026-09-24 |

Full package lists: `docs/env_solved_main.txt`, `docs/env_solved_chai.txt`.

## GROMACS build

Built with spack rather than a hand-rolled CMake invocation so that the full
dependency graph and hash are recorded rather than described.

```bash
source /etc/profile.d/z00_lmod.sh
module load spack/v1.1.1
spack env create -d $PILOT_ROOT/spack-env
spack -e $PILOT_ROOT/spack-env config add "config:install_tree:root:$PILOT_ROOT/spack-install"
spack -e $PILOT_ROOT/spack-env add "gromacs@2025.3 +cuda cuda_arch=80,89 ~mpi +openmp"
spack -e $PILOT_ROOT/spack-env install --fail-fast -j16
```

- **Spec hash:** `45hxd2evh2mxcjqopqa56cr3atcbwlrv`
- **Prefix:** `spack-install/linux-zen4/gromacs-2025.3-45hxd2evh2mxcjqopqa56cr3atcbwlrv`
- **cuda_arch 80,89** covers both GPU types on ORCA: A30 (sm_80) and L40S (sm_89).
- Build time: ~5 minutes on 16 cores, from source (no binary cache available).

Reported configuration:

| Property | Value |
|----------|-------|
| GROMACS version | 2025.3-spack |
| Precision | mixed |
| GPU support | CUDA |
| SIMD instructions | AVX_512 |
| CPU FFT library | fftw-3.3.10 (sse2/avx/avx2/avx512) |
| GPU FFT library | cuFFT |
| C/C++ compiler | GNU 13.4.0 |
| CUDA compiler | nvcc 12.9.41 |

### GPU verification (job 182455, node orcaga01)

| Property | Value |
|----------|-------|
| GPU | NVIDIA L40S |
| Driver | 610.57.04 |
| Compute capability | 8.9 (matches `cuda_arch=89`) |
| CUDA driver / runtime | 13.30 / 12.90 |

`gmx mdrun -version` on a GPU node reports CUDA support, cuFFT, and NBNxM GPU
setup active. The driver (13.30) is newer than the runtime (12.90), which is
the supported direction for CUDA compatibility.
