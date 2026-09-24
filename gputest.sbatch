#!/bin/bash
#SBATCH --job-name=mor_gputest
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/gputest_%j.log
#SBATCH --error=logs/gputest_%j.log
echo "=== node: $(hostname) ==="
nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv
P=/scratch/kuhfeldr-Kuhfeld_temp/spack-install/linux-zen4/gromacs-2025.3-45hxd2evh2mxcjqopqa56cr3atcbwlrv
$P/bin/gmx --version 2>&1 | grep -iE "GROMACS version|GPU support|CUDA driver|CUDA runtime"
echo "=== gmx GPU detection ==="
$P/bin/gmx -quiet mdrun -version 2>&1 | grep -iA4 "GPU" | head -12 || true
