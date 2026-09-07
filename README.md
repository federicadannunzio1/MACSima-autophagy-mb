# MACSima Pipeline — Medulloblastoma Autophagy Study

Single-cell spatial proteomics analysis supporting:

> Giansanti et al., *Autophagy inhibition restores DNAM-1-mediated NK cell immunity in Group 3 medulloblastoma* (under revision)

---

## Overview

This repository contains the full bioinformatics pipeline for the MACSima (Miltenyi Biotec) multiplexed cyclic immunofluorescence analysis of human medulloblastoma tumour sections. The dataset comprises 9 patients (5 Group 3, 4 SHH), 23 regions of interest (ROIs), and approximately 4.7 million segmented cells.

**Six-marker panel:** LC3B · p62 · phospho-mTOR · Nectin-2 · PVR · Ki67

**Key biological question:** Do cells with high autophagic activity (high LC3B, low pMTOR) show reduced DNAM-1 ligand expression (Nectin-2, PVR) in primary tumour tissue?

---

## Repository structure

```
.
├── scripts/               # R analysis scripts (run in numbered order)
│   ├── config.R           # paths, sample metadata, palette — edit before running
│   ├── run_pipeline.R     # master runner for steps 01–07
│   │
│   ├── 01_load_data.R     # load raw MACS iQ View CSV exports
│   ├── 02_qc_normalization.R   # QC filters, arcsinh normalisation (cofactor = 282)
│   ├── 03_integration.R   # Harmony batch correction across patients
│   ├── 04_clustering.R    # Seurat clustering (Leiden)
│   ├── 05_comparison_G3_SHH.R  # differential expression G3 vs SHH
│   ├── 06_correlations.R  # Spearman correlations: autophagy markers vs DNAM-1 ligands
│   ├── 07_visualization.R # UMAP and summary plots
│   │
│   ├── 08_skewness_analysis.R      # LC3B/P62 cytoplasmic skewness distributions
│   ├── 09_skewness_statistics.R    # Wilcoxon tests, punctate cell proportions
│   ├── 10_mtor_skewness_correlation.R  # pMTOR vs skewness correlation per patient
│   ├── 11_mtor_nectin_lc3b_analysis.R  # partial Spearman: pMTOR ~ Nectin2 | LC3B
│   ├── 12_lc3b_mfi_vs_skewness.R   # MFI vs skewness scatter (reviewer response)
│   ├── 13_ki67_sensitivity_analysis.R  # sensitivity: exclude Ki67-high cells
│   │
│   └── figure_panel.R     # composite figure assembly
│
├── slurm/                 # SLURM job scripts for HPC execution
│   ├── run_pipeline.sh         # full pipeline (steps 01–07)
│   ├── run_skewness_analysis.sh
│   ├── run_mtor_skewness.sh
│   └── run_mtor_nectin_lc3b.sh
│
├── setup_packages.R       # install all required R packages
├── data/                  # raw MACS iQ View exports — not tracked (see .gitignore)
├── output/                # pipeline outputs (data + plots) — not tracked
└── logs/                  # SLURM and pipeline logs — not tracked
```

---

## Quick start

### 1. Install dependencies

```r
Rscript setup_packages.R
```

### 2. Configure paths and samples

Edit `scripts/config.R`:
- Set `BASE_DIR` to the project root
- Set `DATA_DIR` to the folder containing patient subdirectories
- Verify `SAMPLES` list matches your data

### 3. Run the core pipeline

```bash
# Run all steps locally
Rscript scripts/run_pipeline.R

# Run specific steps
Rscript scripts/run_pipeline.R --steps 1,2,3

# Restart from a step
Rscript scripts/run_pipeline.R --from 5
```

### 4. Run downstream analyses

Each script in `08_`–`13_` can be run independently after the core pipeline:

```bash
Rscript scripts/08_skewness_analysis.R
Rscript scripts/13_ki67_sensitivity_analysis.R
# etc.
```

### 5. HPC (SLURM)

Submit from the project root:

```bash
sbatch slurm/run_pipeline.sh
sbatch slurm/run_skewness_analysis.sh
```

---

## Input data format

Raw data: per-ROI CSV exports from MACS iQ View. Two files per ROI, row-aligned (one row = one segmented cell):

| File pattern | Content |
|---|---|
| `*_WO.csv` | Mean Fluorescence Intensity (MFI) per protein per cell |
| `*Skewness*.csv` | Cytoplasmic intensity skewness per protein per cell |

Expected directory layout:
```
data/
└── <patient_id>/
    └── <ROI_name>/
        ├── *_WO.csv
        └── *Skewness*.csv
```

---

## Key outputs

| File | Description |
|---|---|
| `output/data/06_spearman_correlations.csv` | Spearman rho per patient per marker pair |
| `output/data/skewness_stats_per_patient.csv` | Median skewness and punctate proportions per patient |
| `output/data/ki67_sensitivity_correlations.csv` | Correlations on Ki67-low cell subset |
| `output/data/mtor_nectin_lc3b_partial_cor.csv` | Partial Spearman results |
| `output/plots/` | All figures, organised by analysis stage |

---

## Methods note — skewness calculation

For each segmented cell, MACS iQ View computes the third standardised moment of the pixel intensity distribution within the cytoplasmic mask:

```
skewness = [(1/N) * sum_i (x_i - mu)^3] / sigma^3
```

where `N` is the number of pixels, `x_i` the raw fluorescence intensity of pixel `i`, and `mu`, `sigma` the within-cell mean and standard deviation. A positive skewness indicates signal concentrated in a minority of bright pixels, consistent with punctate localisation (autophagosomes).

---

## Dependencies

R >= 4.3. Main packages: `Seurat`, `harmony`, `data.table`, `ggplot2`, `dplyr`, `tidyr`, `patchwork`, `ggridges` (optional).

Run `setup_packages.R` to install all dependencies.

---

## Authors

Federica D'Annunzio · Ignazio Caruana · Giansanti Lab

---

## Data availability

Processed single-cell data and source data for figures will be deposited at Zenodo upon acceptance. Analysis scripts are available in this repository. Raw MACSima imaging data are available from the corresponding author upon reasonable request.
