# =============================================================================
# ki67_sensitivity_analysis.R
# Analisi di sensitività: correlazioni autofagia vs ligandi DNAM-1
# escludendo cellule ad alta proliferazione (Ki67 > 75° percentile per paziente)
#
# MOTIVAZIONE (Risposta Reviewer 1A):
#   Il revisore chiede se la correlazione pMTOR ~ Nectin2/PVR potrebbe essere
#   un artefatto della co-espressione in cellule proliferanti (Ki67-high),
#   dove sia pMTOR che i ligandi DNAM-1 sono elevati per ragioni legate al
#   ciclo cellulare piuttosto che all'autofagia.
#
# METODO:
#   1. Per ciascun paziente, calcola il 75° percentile di Ki67 (arcsinh)
#   2. Esclude le cellule sopra quella soglia (Ki67-high)
#   3. Ricalcola le correlazioni Spearman sulle cellule Ki67-low/mid
#   4. Confronta con le correlazioni originali (da 06_spearman_correlations.csv)
#
# Input:  output/data/07_integrated_matrix_arcsinh.tsv
#         output/data/06_spearman_correlations.csv
# Output: output/data/ki67_sensitivity_correlations.csv
#         output/plots/correlations/Ki67_01_sensitivity_comparison.pdf
#         output/plots/correlations/Ki67_02_rho_heatmap_ki67low.pdf
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

message("=== Ki67 SENSITIVITY ANALYSIS — Correlazioni su cellule non-proliferanti ===")

# --------------------------------------------------------------------------
# Carica dati
# Usa 07_integrated_matrix_arcsinh.tsv (già arcsinh-normalizzato, per-cell)
# --------------------------------------------------------------------------

mat_path <- file.path(OUT_DATA, "07_integrated_matrix_arcsinh.tsv")
if (!file.exists(mat_path)) stop("File non trovato: ", mat_path)

message("Caricamento 07_integrated_matrix_arcsinh.tsv...")
df_mat <- fread(mat_path, data.table = FALSE)
message(sprintf("  %s cellule, colonne: %s",
                format(nrow(df_mat), big.mark = ","),
                paste(colnames(df_mat), collapse = ", ")))

# Carica correlazioni originali per confronto
orig_corr <- read.csv(file.path(OUT_DATA, "06_spearman_correlations.csv"),
                      stringsAsFactors = FALSE)

AUTOPHAGY_MARKERS <- c("LC3B", "P62", "pMTOR")
DNAM1_LIGANDS     <- c("PVR", "Nectin2")

PAIRS <- expand.grid(autophagy = AUTOPHAGY_MARKERS,
                     ligand    = DNAM1_LIGANDS,
                     stringsAsFactors = FALSE)

# Colonna Ki67
ki67_col <- colnames(df_mat)[grepl("Ki.?67|Ki_67", colnames(df_mat), ignore.case = TRUE)]
if (length(ki67_col) == 0) stop("Colonna Ki67 non trovata nel file.")
ki67_col <- ki67_col[1]
message(sprintf("  Colonna Ki67 identificata: '%s'", ki67_col))

# --------------------------------------------------------------------------
# Per ogni paziente: identifica Ki67-high (> 75° percentile) ed escludi
# --------------------------------------------------------------------------

message("\nFiltro Ki67 > 75° percentile per paziente...")

ki67_filter_stats <- list()
corr_ki67low      <- list()

for (pid in sort(unique(df_mat$patient_id))) {
  sub   <- df_mat[df_mat$patient_id == pid, ]
  grp   <- unique(sub$group)
  n_all <- nrow(sub)

  ki67_vals <- sub[[ki67_col]]
  ki67_q75  <- quantile(ki67_vals, 0.75, na.rm = TRUE)

  sub_low <- sub[ki67_vals <= ki67_q75, ]
  n_low <- nrow(sub_low)
  pct_excluded <- 100 * (1 - n_low / n_all)

  message(sprintf("  %s (%s): %d cellule totali → %d Ki67-low (%.0f%% escluse)",
                  pid, grp, n_all, n_low, pct_excluded))

  ki67_filter_stats[[pid]] <- data.frame(
    patient_id   = pid,
    group        = grp,
    n_total      = n_all,
    ki67_q75     = round(ki67_q75, 4),
    n_ki67low    = n_low,
    pct_excluded = round(pct_excluded, 1),
    stringsAsFactors = FALSE
  )

  if (n_low < 100) {
    message(sprintf("  [SKIP] %s: troppo poche cellule Ki67-low (%d)", pid, n_low))
    next
  }

  for (i in seq_len(nrow(PAIRS))) {
    x_name <- PAIRS$autophagy[i]
    y_name <- PAIRS$ligand[i]
    x_vals <- sub_low[[x_name]]
    y_vals <- sub_low[[y_name]]

    ct <- cor.test(x_vals, y_vals, method = "spearman", exact = FALSE)

    corr_ki67low[[length(corr_ki67low) + 1]] <- data.frame(
      patient_id = pid,
      group      = grp,
      autophagy  = x_name,
      ligand     = y_name,
      pair       = paste0(x_name, " ~ ", y_name),
      rho        = round(ct$estimate, 4),
      p_value    = ct$p.value,
      n_cells    = n_low,
      stringsAsFactors = FALSE
    )
  }
}

corr_ki67low_df   <- do.call(rbind, corr_ki67low)
ki67_stats_df     <- do.call(rbind, ki67_filter_stats)
rownames(corr_ki67low_df) <- NULL
rownames(ki67_stats_df)   <- NULL

message("\n--- Statistiche filtro Ki67 per paziente ---")
print(as.data.frame(ki67_stats_df))

write.csv(corr_ki67low_df, file.path(OUT_DATA, "ki67_sensitivity_correlations.csv"),
          row.names = FALSE)
message("\nSalvato: ki67_sensitivity_correlations.csv")

# --------------------------------------------------------------------------
# Confronto: rho originale vs rho Ki67-low
# --------------------------------------------------------------------------

orig_sub <- orig_corr %>%
  select(patient_id, group, pair, rho) %>%
  rename(rho_all = rho)

ki67_sub <- corr_ki67low_df %>%
  select(patient_id, pair, rho, n_cells) %>%
  rename(rho_ki67low = rho, n_ki67low = n_cells)

comparison_df <- orig_sub %>%
  left_join(ki67_sub, by = c("patient_id", "pair")) %>%
  mutate(
    delta_rho = round(rho_ki67low - rho_all, 4),
    group     = factor(group, levels = c("G3", "SHH"))
  )

message("\n--- Confronto rho_all vs rho_ki67low (media per gruppo) ---")
summary_comp <- comparison_df %>%
  group_by(pair, group) %>%
  summarise(
    rho_all_mean     = round(mean(rho_all,     na.rm = TRUE), 3),
    rho_ki67low_mean = round(mean(rho_ki67low, na.rm = TRUE), 3),
    delta_mean       = round(mean(delta_rho,   na.rm = TRUE), 3),
    .groups = "drop"
  )
print(as.data.frame(summary_comp))

# --------------------------------------------------------------------------
# PLOT 1 — Confronto rho: all cells vs Ki67-low
# --------------------------------------------------------------------------

message("\nGenerazione plot...")

plot_comp <- comparison_df %>%
  select(patient_id, group, pair, rho_all, rho_ki67low) %>%
  pivot_longer(cols = c(rho_all, rho_ki67low),
               names_to = "subset", values_to = "rho") %>%
  mutate(
    subset = recode(subset,
                    rho_all     = "All cells",
                    rho_ki67low = "Ki67-low (excl. proliferating)"),
    subset = factor(subset, levels = c("All cells", "Ki67-low (excl. proliferating)"))
  )

p1 <- ggplot(plot_comp, aes(x = rho, y = patient_id, colour = group, shape = subset)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_line(aes(group = paste(patient_id, pair)), colour = "grey70", linewidth = 0.5) +
  geom_point(size = 3, alpha = 0.9) +
  facet_wrap(~ pair, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = PALETTE_GROUP) +
  scale_shape_manual(values = c("All cells" = 16, "Ki67-low (excl. proliferating)" = 15)) +
  labs(
    title    = "Sensitivity analysis: correlations with and without Ki67-high cells",
    subtitle = "Lines connect the same patient; triangles = Ki67-low subset (excl. top 25% proliferating cells)",
    x        = "Spearman rho",
    y        = NULL,
    colour   = "Group",
    shape    = "Cell subset"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(size = 9))

ggsave(file.path(OUT_PLOTS_CORR, "Ki67_01_sensitivity_comparison.pdf"),
       p1, width = 12, height = 10)
message("  Salvato: Ki67_01_sensitivity_comparison.pdf")

# --------------------------------------------------------------------------
# PLOT 2 — Heatmap rho Ki67-low per gruppo
# --------------------------------------------------------------------------

corr_mean_ki67low <- corr_ki67low_df %>%
  group_by(group, autophagy, ligand) %>%
  summarise(rho_mean = round(mean(rho), 3), .groups = "drop") %>%
  mutate(group = factor(group, levels = c("G3", "SHH")))

p2 <- ggplot(corr_mean_ki67low,
             aes(x = ligand, y = autophagy, fill = rho_mean)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = rho_mean), size = 4) +
  facet_wrap(~ group) +
  scale_fill_gradient2(low = "#1A5276", mid = "white", high = "#C0392B",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman rho") +
  labs(
    title    = "Spearman rho — Ki67-low cells only (excl. top 25% proliferating)",
    subtitle = "Mean rho across patients per group",
    x        = "DNAM-1 ligand",
    y        = "Autophagy marker"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_PLOTS_CORR, "Ki67_02_rho_heatmap_ki67low.pdf"),
       p2, width = 8, height = 5)
message("  Salvato: Ki67_02_rho_heatmap_ki67low.pdf")

# --------------------------------------------------------------------------
# STAMPA FINALE: Valutazione robustezza
# --------------------------------------------------------------------------

message("\n=== VALUTAZIONE ROBUSTEZZA ===")

# Per la coppia principale: pMTOR ~ Nectin2 e pMTOR ~ PVR
for (pr in c("pMTOR ~ Nectin2", "pMTOR ~ PVR")) {
  sub <- comparison_df[comparison_df$pair == pr, ]
  message(sprintf("\n%s:", pr))
  message(sprintf("  rho ALL:      G3 mean=%.3f, SHH mean=%.3f",
                  mean(sub$rho_all[sub$group == "G3"], na.rm = TRUE),
                  mean(sub$rho_all[sub$group == "SHH"], na.rm = TRUE)))
  message(sprintf("  rho Ki67-low: G3 mean=%.3f, SHH mean=%.3f",
                  mean(sub$rho_ki67low[sub$group == "G3"], na.rm = TRUE),
                  mean(sub$rho_ki67low[sub$group == "SHH"], na.rm = TRUE)))
  message(sprintf("  delta medio:  G3=%.3f, SHH=%.3f",
                  mean(sub$delta_rho[sub$group == "G3"], na.rm = TRUE),
                  mean(sub$delta_rho[sub$group == "SHH"], na.rm = TRUE)))

  robust <- all(abs(sub$delta_rho) < 0.1, na.rm = TRUE)
  message(sprintf("  Robusto (|delta| < 0.1 per tutti i pazienti): %s", robust))
}

message("\n=== Ki67 SENSITIVITY ANALYSIS completata ===")
