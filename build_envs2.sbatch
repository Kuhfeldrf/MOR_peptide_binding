#!/bin/bash
#SBATCH --job-name=mor_env2
#SBATCH --partition=short
#SBATCH --time=03:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=logs/env_build2_%j.log
#SBATCH --error=logs/env_build2_%j.log
set -x
bash /scratch/kuhfeldr-Kuhfeld_temp/build_envs2.sh
