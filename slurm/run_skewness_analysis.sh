#!/bin/bash
# =============================================================================
# run_skewness_analysis.sh
# SLURM job script: LC3B/P62 cytoplasmic skewness analysis and statistics
#
# Runs scripts 08 and 09 in sequence:
#   08_skewness_analysis.R   — per-cell skewness distributions and plots
#   09_skewness_statistics.R — Wilcoxon tests and punctate proportions
#
# Usage (submit from the project root directory):
#   sbatch slurm/run_skewness_analysis.sh
#
# Adjust conda path and #SBATCH directives to match your cluster.
# =============================================================================

#SBATCH --job-name=skewness
#SBATCH --output=logs/skewness_%j.log
#SBATCH --error=logs/skewness_%j.log
#SBATCH --time=00:30:00
#SBATCH --mem=16G
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

Rscript "${PROJECT_DIR}/scripts/08_skewness_analysis.R"
Rscript "${PROJECT_DIR}/scripts/09_skewness_statistics.R"
