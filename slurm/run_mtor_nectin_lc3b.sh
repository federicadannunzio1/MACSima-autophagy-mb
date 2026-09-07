#!/bin/bash
# =============================================================================
# run_mtor_nectin_lc3b.sh
# SLURM job script: partial Spearman correlation pMTOR ~ Nectin2 | LC3B
#
# Usage (submit from the project root directory):
#   sbatch slurm/run_mtor_nectin_lc3b.sh
#
# Adjust conda path and #SBATCH directives to match your cluster.
# =============================================================================

#SBATCH --job-name=partial_cor
#SBATCH --output=logs/partial_cor_%j.log
#SBATCH --error=logs/partial_cor_%j.log
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=BEGIN,END,FAIL
# #SBATCH --mail-user=your.email@institution.edu   # uncomment and set your email

PROJECT_DIR="${SLURM_SUBMIT_DIR}"

source /lustre/software/anaconda/2022.10_all/etc/profile.d/conda.sh
conda activate seurat_env

export MACSIMA_DIR="${PROJECT_DIR}"
export MACSIMA_SCRIPTS_DIR="${PROJECT_DIR}/scripts"

mkdir -p "${PROJECT_DIR}/logs"

Rscript "${PROJECT_DIR}/scripts/11_mtor_nectin_lc3b_analysis.R"
