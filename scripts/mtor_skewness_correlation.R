# =============================================================================
# mtor_skewness_correlation.R
# Correlazione tra pMTOR (MFI arcsinh) e skewness di LC3B / P62
#
# Ipotesi biologica:
#   pMTOR alto -> mTOR attivo -> autofagia soppressa -> LC3B diffuso
#   => correlazione NEGATIVA pMTOR ~ LC3B skewness
#   pMTOR alto -> P62 non degradato -> P62 che si accumula
#   => correlazione POSITIVA pMTOR ~ P62 skewness
#
# I file *Skewness*.csv di MACS iQ View contengono nella stessa riga
# sia la skewness di LC3B/P62 che l'MFI grezzo di pMTOR, permettendo
# l'analisi senza merge tra file distinti.
#
# Output:
#   output/data/mtor_skewness_correlations.csv
#   output/plots/mTOR_Skew_01_scatter_LC3B.pdf
#   output/plots/mTOR_Skew_02_scatter_P62.pdf
#   output/plots/mTOR_Skew_03_rho_summary.pdf
#   output/plots/mTOR_Skew_04_rho_by_group.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

.sd <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

message("=== pMTOR ~ LC3B/P62 SKEWNESS CORRELATION ===")

# --------------------------------------------------------------------------
# Configurazione
# --------------------------------------------------------------------------

# Cofactor usato nella pipeline principale (02_qc_normalization.R)
# Riletto dal CSV per garantire coerenza
cofactor_csv <- file.path(OUT_DATA, "02_cofactor_quantiles.csv")
if (file.exists(cofactor_csv)) {
  q_df     <- read.csv(cofactor_csv)
  COFACTOR <- median(q_df$p50)
  message(sprintf("Cofactor riletto da 02_cofactor_quantiles.csv: %.0f", COFACTOR))
} else {
  COFACTOR <- 282
  message(sprintf("02_cofactor_quantiles.csv non trovato — uso cofactor = %d", COFACTOR))
}

# Colonne nei file skewness di MACS iQ View
COL_MTOR_RAW  <- "Phospo mTor 49F9 Cyto Exp"
COL_LC3B_SKEW <- "LC3B Cytoplasm Intensity Skewness"
COL_P62_SKEW  <- "P62 Cytoplasm Intensity Skewness"
REQUIRED_COLS <- c(COL_MTOR_RAW, COL_LC3B_SKEW, COL_P62_SKEW)

# --------------------------------------------------------------------------
# Funzione: trova il file skewness con tutte le colonne necessarie
# --------------------------------------------------------------------------

find_skewness_file <- function(roi_dir) {
  if (!dir.exists(roi_dir)) return(NULL)
  files <- list.files(roi_dir, pattern = "Sk?ewness", full.names = TRUE,
                      ignore.case = TRUE)
  if (length(files) == 0) return(NULL)
  for (f in files) {
    header <- tryCatch(colnames(fread(f, nrows = 0)),
                       error = function(e) character(0))
    if (all(REQUIRED_COLS %in% header)) return(f)
  }
  return(NULL)
}

# --------------------------------------------------------------------------
# Caricamento: scansiona tutte le ROI di ogni paziente in DATA_DIR
# (come skewness_analysis.R: include ROI non presenti in SAMPLES se
#  il file skewness esiste e ha tutte le colonne necessarie)
# --------------------------------------------------------------------------

patient_dirs <- list.dirs(DATA_DIR, recursive = FALSE, full.names = TRUE)
cell_list    <- list()
skipped_rois <- character(0)

for (pat_dir in sort(patient_dirs)) {
  pid <- basename(pat_dir)

  # Ricava il gruppo da SAMPLES; se non trovato, salta (e.g. SHH_25-6667)
  grp <- NULL
  for (s in SAMPLES) {
    if (s$patient_id == pid) { grp <- s$group; break }
  }
  if (is.null(grp)) {
    message(sprintf("\n[SKIP] %s — non in SAMPLES (paziente escluso)", pid))
    next
  }

  message(sprintf("\n%s (%s):", pid, grp))

  roi_dirs <- list.dirs(pat_dir, recursive = FALSE, full.names = TRUE)
  for (roi_dir in sort(roi_dirs)) {
    roi_name <- basename(roi_dir)
    skew_file <- find_skewness_file(roi_dir)

    if (is.null(skew_file)) {
      message(sprintf("  [SKIP] %s: nessun file skewness con pMTOR+LC3B+P62", roi_name))
      skipped_rois <- c(skipped_rois, paste0(pid, "_", roi_name))
      next
    }

    dt <- tryCatch(
      fread(skew_file, data.table = FALSE),
      error = function(e) {
        message(sprintf("  [ERROR] %s: %s", roi_name, e$message))
        NULL
      }
    )
    if (is.null(dt)) next

    # Verifica colonne richieste (dopo lettura completa, evita problemi di encoding con select=)
    missing_cols <- setdiff(REQUIRED_COLS, colnames(dt))
    if (length(missing_cols) > 0) {
      message(sprintf("  [SKIP] %s: colonne mancanti: %s", roi_name,
                      paste(missing_cols, collapse = ", ")))
      skipped_rois <- c(skipped_rois, paste0(pid, "_", roi_name))
      next
    }
    dt <- dt[, REQUIRED_COLS, drop = FALSE]

    # Converti a numerico e rimuovi NA
    dt <- as.data.frame(lapply(dt, as.numeric))
    na_mask <- rowSums(is.na(dt)) > 0
    if (any(na_mask)) {
      message(sprintf("  [WARN] %s: %d righe con NA rimosse", roi_name, sum(na_mask)))
      dt <- dt[!na_mask, , drop = FALSE]
    }

    if (nrow(dt) == 0) {
      message(sprintf("  [SKIP] %s: nessuna cella valida", roi_name))
      next
    }

    # Trasformazione arcsinh su pMTOR
    dt$pMTOR_arcsinh  <- asinh(dt[[COL_MTOR_RAW]] / COFACTOR)
    dt$LC3B_skewness  <- dt[[COL_LC3B_SKEW]]
    dt$P62_skewness   <- dt[[COL_P62_SKEW]]
    dt$patient_id     <- pid
    dt$group          <- grp
    dt$roi_id         <- paste0(pid, "_", roi_name)

    cell_list[[paste0(pid, "_", roi_name)]] <- dt[, c("patient_id", "group", "roi_id",
                                                       "pMTOR_arcsinh",
                                                       "LC3B_skewness",
                                                       "P62_skewness")]
    message(sprintf("  %s -> %s (%d cellule)", roi_name, basename(skew_file), nrow(dt)))
  }
}

if (length(cell_list) == 0) stop("Nessun dato caricato.")

cells <- do.call(rbind, cell_list)
rownames(cells) <- NULL

message(sprintf("\nTotale cellule: %d da %d ROI", nrow(cells), length(cell_list)))
if (length(skipped_rois) > 0)
  message(sprintf("ROI senza dati (%d): %s", length(skipped_rois),
                  paste(skipped_rois, collapse = ", ")))

# --------------------------------------------------------------------------
# Correlazioni Spearman per paziente
# --------------------------------------------------------------------------

message("\nCalcolo correlazioni Spearman per paziente...")

cor_rows <- list()
for (pid in unique(cells$patient_id)) {
  sub <- cells[cells$patient_id == pid, ]
  grp <- sub$group[1]
  n   <- nrow(sub)

  r_lc3b <- cor.test(sub$pMTOR_arcsinh, sub$LC3B_skewness,
                     method = "spearman", exact = FALSE)
  r_p62  <- cor.test(sub$pMTOR_arcsinh, sub$P62_skewness,
                     method = "spearman", exact = FALSE)

  cor_rows[[pid]] <- data.frame(
    patient_id        = pid,
    group             = grp,
    n_cells           = n,
    rho_pMTOR_LC3Bskew = round(r_lc3b$estimate, 4),
    p_pMTOR_LC3Bskew   = signif(r_lc3b$p.value, 3),
    rho_pMTOR_P62skew  = round(r_p62$estimate, 4),
    p_pMTOR_P62skew    = signif(r_p62$p.value, 3),
    stringsAsFactors = FALSE
  )
}

cor_df <- do.call(rbind, cor_rows)
rownames(cor_df) <- NULL

message("\nCorrelazioni Spearman pMTOR ~ Skewness per paziente:")
print(cor_df[, c("patient_id", "group", "n_cells",
                 "rho_pMTOR_LC3Bskew", "rho_pMTOR_P62skew")])

# Riepilogo per gruppo
for (grp in c("G3", "SHH")) {
  sub <- cor_df[cor_df$group == grp, ]
  message(sprintf("\n  %s (n=%d pazienti):", grp, nrow(sub)))
  message(sprintf("    pMTOR ~ LC3B skewness: mean rho = %.3f  (range %.3f - %.3f)",
                  mean(sub$rho_pMTOR_LC3Bskew),
                  min(sub$rho_pMTOR_LC3Bskew),
                  max(sub$rho_pMTOR_LC3Bskew)))
  message(sprintf("    pMTOR ~ P62 skewness:  mean rho = %.3f  (range %.3f - %.3f)",
                  mean(sub$rho_pMTOR_P62skew),
                  min(sub$rho_pMTOR_P62skew),
                  max(sub$rho_pMTOR_P62skew)))
}

write.csv(cor_df, file.path(OUT_DATA, "mtor_skewness_correlations.csv"),
          row.names = FALSE)
message("\nSalvato: mtor_skewness_correlations.csv")

# --------------------------------------------------------------------------
# Wilcoxon G3 vs SHH (esplorativo, n piccolo)
# --------------------------------------------------------------------------

message("\nWilcoxon G3 vs SHH (esplorativo):")
for (metric in c("rho_pMTOR_LC3Bskew", "rho_pMTOR_P62skew")) {
  g3  <- cor_df[cor_df$group == "G3",  metric]
  shh <- cor_df[cor_df$group == "SHH", metric]
  w   <- wilcox.test(g3, shh, exact = FALSE)
  message(sprintf("  %s: W=%.0f  p=%.4f", metric, w$statistic, w$p.value))
}

# --------------------------------------------------------------------------
# PLOT 1 — Scatter pMTOR vs LC3B skewness per paziente
# --------------------------------------------------------------------------

message("\nGenerazione plot...")

# Subsample per rendere i scatter leggibili (max 5000 celle per paziente)
set.seed(SEED)
cells_sub <- cells %>%
  group_by(patient_id) %>%
  slice_sample(n = min(5000, n())) %>%
  ungroup() %>%
  left_join(cor_df[, c("patient_id", "rho_pMTOR_LC3Bskew")], by = "patient_id") %>%
  mutate(label = sprintf("%s\nrho=%.3f", patient_id, rho_pMTOR_LC3Bskew),
         group = factor(group, levels = c("G3", "SHH")))

p1 <- ggplot(cells_sub, aes(x = pMTOR_arcsinh, y = LC3B_skewness, colour = group)) +
  geom_point(alpha = 0.15, size = 0.3) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.6) +
  facet_wrap(~ label, ncol = 4) +
  scale_colour_manual(values = PALETTE_GROUP) +
  labs(title = "pMTOR (arcsinh) vs LC3B skewness",
       subtitle = "Spearman rho shown per patient | subsample 5,000 cells",
       x = "pMTOR arcsinh(MFI / 282)",
       y = "LC3B cytoplasmic skewness",
       colour = "Group") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 7),
        legend.position = "bottom")

ggsave(file.path(OUT_PLOTS, "mTOR_Skew_01_scatter_LC3B.pdf"),
       p1, width = 14, height = 8)
message("Salvato: mTOR_Skew_01_scatter_LC3B.pdf")

# --------------------------------------------------------------------------
# PLOT 2 — Scatter pMTOR vs P62 skewness per paziente
# --------------------------------------------------------------------------

cells_sub2 <- cells %>%
  group_by(patient_id) %>%
  slice_sample(n = min(5000, n())) %>%
  ungroup() %>%
  left_join(cor_df[, c("patient_id", "rho_pMTOR_P62skew")], by = "patient_id") %>%
  mutate(label = sprintf("%s\nrho=%.3f", patient_id, rho_pMTOR_P62skew),
         group = factor(group, levels = c("G3", "SHH")))

p2 <- ggplot(cells_sub2, aes(x = pMTOR_arcsinh, y = P62_skewness, colour = group)) +
  geom_point(alpha = 0.15, size = 0.3) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.6) +
  facet_wrap(~ label, ncol = 4) +
  scale_colour_manual(values = PALETTE_GROUP) +
  labs(title = "pMTOR (arcsinh) vs P62 skewness",
       subtitle = "Spearman rho shown per patient | subsample 5,000 cells",
       x = "pMTOR arcsinh(MFI / 282)",
       y = "P62 cytoplasmic skewness",
       colour = "Group") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 7),
        legend.position = "bottom")

ggsave(file.path(OUT_PLOTS, "mTOR_Skew_02_scatter_P62.pdf"),
       p2, width = 14, height = 8)
message("Salvato: mTOR_Skew_02_scatter_P62.pdf")

# --------------------------------------------------------------------------
# PLOT 3 — rho per paziente: LC3B vs P62, ordinato per gruppo
# --------------------------------------------------------------------------

cor_long <- cor_df %>%
  select(patient_id, group, rho_pMTOR_LC3Bskew, rho_pMTOR_P62skew) %>%
  pivot_longer(cols = c(rho_pMTOR_LC3Bskew, rho_pMTOR_P62skew),
               names_to = "marker",
               values_to = "rho") %>%
  mutate(
    marker = recode(marker,
                    rho_pMTOR_LC3Bskew = "pMTOR ~ LC3B skewness",
                    rho_pMTOR_P62skew  = "pMTOR ~ P62 skewness"),
    group      = factor(group, levels = c("G3", "SHH")),
    patient_id = reorder(patient_id, ifelse(group == "G3", 1, 2))
  )

p3 <- ggplot(cor_long, aes(x = rho, y = patient_id, fill = group, shape = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_point(aes(colour = group), size = 3.5) +
  facet_wrap(~ marker, ncol = 2) +
  scale_colour_manual(values = PALETTE_GROUP) +
  scale_fill_manual(values = PALETTE_GROUP) +
  labs(title = "Spearman rho: pMTOR vs LC3B/P62 skewness per patient",
       subtitle = "Dashed line = rho 0 (no correlation)",
       x = "Spearman rho",
       y = NULL,
       colour = "Group") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_PLOTS, "mTOR_Skew_03_rho_per_patient.pdf"),
       p3, width = 10, height = 5)
message("Salvato: mTOR_Skew_03_rho_per_patient.pdf")

# --------------------------------------------------------------------------
# PLOT 4 — rho per gruppo: boxplot G3 vs SHH
# --------------------------------------------------------------------------

p4 <- ggplot(cor_long, aes(x = group, y = rho, colour = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_boxplot(width = 0.4, outlier.shape = NA, fill = NA, linewidth = 0.7) +
  geom_jitter(width = 0.1, size = 3) +
  facet_wrap(~ marker, ncol = 2) +
  scale_colour_manual(values = PALETTE_GROUP) +
  labs(title = "pMTOR ~ skewness: G3 vs SHH",
       subtitle = "Each dot = one patient (exploratory, n=4 per group)",
       x = NULL,
       y = "Spearman rho",
       colour = "Group") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(OUT_PLOTS, "mTOR_Skew_04_rho_by_group.pdf"),
       p4, width = 8, height = 5)
message("Salvato: mTOR_Skew_04_rho_by_group.pdf")

message("\n=== mTOR SKEWNESS CORRELATION completata ===")
message(sprintf("Plot in: %s", OUT_PLOTS))
