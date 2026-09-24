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

Populated in step 3 once the environment solves.

| Tool | Version | Source | Date |
|------|---------|--------|------|
| _pending_ | | | |
