#!/bin/bash
# =============================================================================
# run_mtor_skewness.sh
# SLURM job script: pMTOR vs LC3B/P62 skewness correlation per patient
#
# Usage (submit from the project root directory):
#   sbatch slurm/run_mtor_skewness.sh
#
# Adjust conda path and #SBATCH directives to match your cluster.
# =============================================================================

#SBATCH --job-name=mtor_skew
#SBATCH --output=logs/mtor_skew_%j.log
#SBATCH --error=logs/mtor_skew_%j.log
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

Rscript "${PROJECT_DIR}/scripts/10_mtor_skewness_correlation.R"
