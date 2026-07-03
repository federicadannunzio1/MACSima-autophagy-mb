#!/bin/bash
#SBATCH --job-name=mtor_skew
#SBATCH --output=/lustre/home/gfiscon/projects/MACSima_pipeline/logs/mtor_skew_%j.log
#SBATCH --error=/lustre/home/gfiscon/projects/MACSima_pipeline/logs/mtor_skew_%j.log
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=federica.dannunzio@uniroma1.it

PROJECT_DIR="/lustre/home/gfiscon/projects/MACSima_pipeline"

source /lustre/software/anaconda/2022.10_all/etc/profile.d/conda.sh
conda activate seurat_env

export MACSIMA_DIR="${PROJECT_DIR}"
export MACSIMA_SCRIPTS_DIR="${PROJECT_DIR}/scripts"

Rscript "${PROJECT_DIR}/scripts/mtor_skewness_correlation.R"
