# MACSima Pipeline — Medulloblastoma Autophagy Study

Single-cell spatial proteomics analysis supporting:

> Giansanti et al., *Autophagy inhibition restores DNAM-1-mediated NK cell immunity and reprograms tumor microenvironment in Group 3 medulloblastoma* (under revision, 2026)

---

## Biological background

This pipeline analyses MACSima multiplexed cyclic immunofluorescence data from human medulloblastoma tumour sections. The central question is whether cells with high autophagic activity show reduced expression of DNAM-1 ligands (Nectin-2, PVR) in primary tumour tissue, and whether this relationship differs between Group 3 (G3) and SHH medulloblastoma subtypes.

**Dataset:**
- 9 patients: 5 Group 3, 4 SHH
- 23 regions of interest (ROIs)
- ~4.7 million segmented cells
- 6-marker panel: LC3B · p62 · phospho-mTOR · Nectin-2 · PVR · Ki67

---

## Repository structure

```
.
├── README.md
├── scripts/
│   ├── config.R                # paths, sample list, parameters — edit this first
│   ├── run_pipeline.R          # master runner for the core pipeline (steps 01-07)
│   ├── 00_setup_packages.R     # install all R dependencies — run once before anything else
│   │
│   │   ── Core pipeline (run in order via run_pipeline.R) ──────────────
│   ├── 01_load_data.R          # load and validate raw MACS iQ View CSV exports
│   ├── 02_qc_normalization.R   # QC, arcsinh normalisation (cofactor = 282)
│   ├── 03_integration.R        # Harmony batch correction across patients
│   ├── 04_clustering.R         # Seurat clustering (Leiden, resolution = 0.1)
│   ├── 05_comparison_G3_SHH.R  # differential expression G3 vs SHH
│   ├── 06_correlations.R       # Spearman correlations: autophagy vs DNAM-1 ligands
│   ├── 07_visualization.R      # UMAP and summary figures
│   │
│   │   ── Downstream analyses (run independently after step 06) ─────────
│   ├── 08_skewness_analysis.R       # LC3B/P62 cytoplasmic skewness distributions
│   ├── 09_skewness_statistics.R     # Wilcoxon tests, punctate cell proportions
│   ├── 10_mtor_skewness_correlation.R  # pMTOR vs skewness per patient
│   ├── 11_mtor_nectin_lc3b_analysis.R  # partial Spearman: pMTOR ~ Nectin2 | LC3B
│   ├── 12_lc3b_mfi_vs_skewness.R   # MFI vs skewness scatter (4M cells)
│   ├── 13_ki67_sensitivity_analysis.R  # re-run correlations excluding Ki67-high cells
│   │
│   └── 14_figure_panel.R      # composite figure assembly for the manuscript
│
└── slurm/                     # SLURM job scripts for HPC — submit from project root
    ├── run_pipeline.sh             # full pipeline (steps 01-07)
    ├── run_skewness_analysis.sh    # runs scripts 08 and 09
    ├── run_mtor_skewness.sh        # runs script 10
    └── run_mtor_nectin_lc3b.sh    # runs script 11
```

---

## Requirements

- **R >= 4.3**
- All package dependencies installed via `scripts/00_setup_packages.R` (see below)
- Raw MACS iQ View CSV exports in the `data/` directory (see *Data* section)

Main packages: `Seurat`, `harmony`, `data.table`, `ggplot2`, `dplyr`, `tidyr`, `patchwork`, `ggridges`, `clustree`.

---

## How to reproduce the analyses

### Step 0 — Clone the repository

```bash
git clone https://github.com/federicadannunzio1/MACSima-autophagy-mb.git
cd MACSima-autophagy-mb
```

### Step 1 — Install R dependencies

```r
Rscript scripts/00_setup_packages.R
```

This checks which packages are installed and installs any that are missing from CRAN and Bioconductor. Run once before the first execution.

### Step 2 — Place the raw data

The `data/` directory is not included in this repository (see *Data availability*). Once you have access to the raw data, organise it as follows:

```
data/
├── Gr3_23-3017/
│   ├── ROI1/
│   │   ├── <export>_WO.csv           # MFI per cell per protein
│   │   └── <export>_Skewness.csv     # cytoplasmic skewness per cell per protein
│   └── ROI2/
│       └── ...
├── Gr3_23-3106/
│   └── ...
└── SHH_24-8477/
    └── ...
```

Each ROI folder contains exactly two CSV files exported from MACS iQ View:
- `*_WO.csv`: mean fluorescence intensity (MFI) per protein per segmented cell
- `*Skewness*.csv`: cytoplasmic intensity skewness per protein per segmented cell

The two files are row-aligned: row *i* in both files corresponds to the same segmented cell.

### Step 3 — Configure paths

Open `scripts/config.R` and verify:

```r
BASE_DIR <- Sys.getenv("MACSIMA_DIR", unset = getwd())
```

By default `BASE_DIR` is the directory from which you launch R. If you always run from the project root (recommended), no changes are needed. Alternatively, set the environment variable before running:

```bash
export MACSIMA_DIR=/absolute/path/to/MACSima-autophagy-mb
```

The `SAMPLES` list in `config.R` defines which patients and ROIs are included. The current configuration matches the 9 patients and 23 ROIs in the paper. Do not modify unless you are adapting the pipeline to a different dataset.

### Step 4 — Run the core pipeline (steps 01-07)

**Locally:**

```bash
# Full pipeline
Rscript scripts/run_pipeline.R

# Selected steps only
Rscript scripts/run_pipeline.R --steps 1,2,3

# Restart from a specific step
Rscript scripts/run_pipeline.R --from 5
```

**On an HPC cluster with SLURM:**

```bash
# Always submit from the project root
cd /path/to/MACSima-autophagy-mb
sbatch slurm/run_pipeline.sh
```

Edit `slurm/run_pipeline.sh` to adjust the conda environment path (`conda activate seurat_env`) and SLURM resource directives (`--mem`, `--time`, `--cpus-per-task`) to match your cluster.

The pipeline writes intermediate outputs to `output/data/` and plots to `output/plots/`. Log files go to `logs/`.

> **Note:** Steps 01-07 require significant RAM. Step 03 (Harmony integration on ~4.7M cells) needs at least 64 GB. Step 04 (clustering) needs at least 128 GB. On a laptop, run on a subset of patients first.

### Step 5 — Run the downstream analyses (steps 08-13)

These scripts are independent of each other but all depend on outputs from the core pipeline (in particular `output/data/` files produced by steps 01-06). Run them in any order after the core pipeline is complete.

```bash
# Skewness distributions and statistics (run together)
Rscript scripts/08_skewness_analysis.R
Rscript scripts/09_skewness_statistics.R

# pMTOR vs skewness correlation
Rscript scripts/10_mtor_skewness_correlation.R

# Partial Spearman: pMTOR ~ Nectin2 controlling for LC3B
Rscript scripts/11_mtor_nectin_lc3b_analysis.R

# LC3B MFI vs skewness scatter (reads raw data directly, ~4M cells)
Rscript scripts/12_lc3b_mfi_vs_skewness.R

# Sensitivity analysis: exclude Ki67-high cells, re-run correlations
Rscript scripts/13_ki67_sensitivity_analysis.R
```

**On HPC:**

```bash
sbatch slurm/run_skewness_analysis.sh    # runs 08 + 09
sbatch slurm/run_mtor_skewness.sh        # runs 10
sbatch slurm/run_mtor_nectin_lc3b.sh     # runs 11
```

Scripts 12 and 13 do not have dedicated SLURM scripts; they can be run locally or with a simple `srun`/`sbatch` one-liner.

### Step 6 — Assemble the manuscript figure panel

```bash
Rscript scripts/14_figure_panel.R
```

Requires outputs from `06_correlations.R` and `08_skewness_analysis.R`. Produces `output/plots/figures/Figure_Panel_combined.pdf` and the individual panels A and B.

---

## Script dependency chain

```
raw data (data/)
    │
    ▼
01_load_data.R  ──────────────────────────────────────────────► 01_cells_raw.rds
02_qc_normalization.R  ───────────────────────────────────────► 02_seurat_normalised.rds
                                                                  02_cofactor_quantiles.csv
03_integration.R  ────────────────────────────────────────────► 03_seurat_integrated.rds
04_clustering.R  ─────────────────────────────────────────────► 04_seurat_clustered.rds
05_comparison_G3_SHH.R  ──────────────────────────────────────► 05_*.csv
06_correlations.R  ───────────────────────────────────────────► 06_spearman_correlations.csv
07_visualization.R  (final summary plots)

    ├── 08_skewness_analysis.R  ◄── raw data
    │       └── 09_skewness_statistics.R  ◄── skewness_cells.csv
    │
    ├── 10_mtor_skewness_correlation.R  ◄── 04_seurat_clustered.rds
    ├── 11_mtor_nectin_lc3b_analysis.R  ◄── 04_seurat_clustered.rds
    │                                        06_spearman_correlations.csv
    │
    ├── 12_lc3b_mfi_vs_skewness.R  ◄── raw data
    │
    ├── 13_ki67_sensitivity_analysis.R  ◄── 07_integrated_matrix_arcsinh.tsv
    │                                        06_spearman_correlations.csv
    │
    └── 14_figure_panel.R  ◄── 06_spearman_correlations.csv
                                skewness_cells.csv (from 08)
                                ► Figure_Panel_combined.pdf
```

---

## Key outputs

All outputs are written to `output/` (not tracked by git).

| File | Produced by | Description |
|---|---|---|
| `output/data/06_spearman_correlations.csv` | `06_correlations.R` | Spearman rho per patient per marker pair |
| `output/data/skewness_stats_per_patient.csv` | `09_skewness_statistics.R` | Median skewness, punctate cell proportions per patient |
| `output/data/ki67_sensitivity_correlations.csv` | `13_ki67_sensitivity_analysis.R` | Correlations on Ki67-low cell subset |
| `output/data/mtor_nectin_lc3b_partial_cor.csv` | `11_mtor_nectin_lc3b_analysis.R` | Partial Spearman results |
| `output/plots/correlations/` | `06`, `11`, `13` | Correlation plots and heatmaps |
| `output/plots/skewness/` | `08`, `09`, `12` | Skewness distributions and MFI vs skewness scatters |

---

## Methods note — normalisation and skewness

**arcsinh normalisation:** Raw MFI values are transformed as `arcsinh(MFI / cofactor)` with `cofactor = 282` (median of per-protein, per-patient 50th-percentile values across all 54 protein–patient combinations). This compresses the dynamic range while preserving differences at low expression levels.

**Cytoplasmic skewness:** For each segmented cell, MACS iQ View computes the third standardised moment of the pixel intensity distribution within the cytoplasmic mask:

```
skewness = [(1/N) * sum_i (x_i - mu)^3] / sigma^3
```

where `N` is the number of pixels in the cytoplasmic mask, `x_i` is the raw fluorescence intensity of pixel `i`, and `mu`, `sigma` are the within-cell mean and standard deviation. A positive skewness indicates that most pixels have low intensity while a minority have high intensity — consistent with signal concentrated in discrete puncta (autophagosomes).

---

## Data availability

Raw MACSima imaging data and processed single-cell CSV files will be deposited upon acceptance. Analysis scripts are available in this repository.

To request early access to the data, contact the corresponding author.

---

## Citation

If you use this code, please cite the associated paper (citation will be updated upon acceptance).

