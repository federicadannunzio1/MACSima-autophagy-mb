# =============================================================================
# skewness_statistics.R
# Statistiche riassuntive per la distribuzione di skewness LC3B e P62
#
# Approccio:
#   - Unità statistica = paziente (n=4 G3, n=4 SHH)
#   - Test di confronto G3 vs SHH: Wilcoxon rank-sum (Mann-Whitney) — non
#     parametrico, appropriato per n piccolo
#   - Effect size: rank-biserial correlation (r_rb = 1 - 2W / (n1 * n2)),
#     interpretato come r < 0.3 piccolo, 0.3-0.5 medio, > 0.5 grande
#   - Per le proporzioni punctate: bootstrap 95% CI gia' calcolato in
#     skewness_analysis.R; qui si riassumono a livello di gruppo
#
# Input:  output/data/skewness_patient_medians.csv
#         output/data/skewness_punctate_proportion.csv
# Output: output/data/skewness_stats_per_patient.csv
#         output/data/skewness_stats_group_summary.csv
#         output/data/skewness_stats_wilcoxon.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

.sd <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

message("=== SKEWNESS STATISTICS ===")

# --------------------------------------------------------------------------
# Carica dati
# --------------------------------------------------------------------------

medians_df <- read.csv(file.path(OUT_DATA, "skewness_patient_medians.csv"),
                       stringsAsFactors = FALSE)
punctate_df <- read.csv(file.path(OUT_DATA, "skewness_punctate_proportion.csv"),
                        stringsAsFactors = FALSE)

# --------------------------------------------------------------------------
# TABELLA 1: Riepilogo per paziente
# Unisce mediane + proporzioni punctate a tre soglie + CI bootstrap
# --------------------------------------------------------------------------

prop_wide <- punctate_df %>%
  select(patient_id, threshold_lab, prop_punctate, ci_lower, ci_upper) %>%
  pivot_wider(
    names_from  = threshold_lab,
    values_from = c(prop_punctate, ci_lower, ci_upper),
    names_sep   = "_"
  )

# Rinomina colonne per chiarezza
rename_map <- c(
  "prop_punctate_skewness > 0"   = "prop_LC3B_pos",
  "prop_punctate_skewness > 0.5" = "prop_LC3B_gt0.5",
  "prop_punctate_skewness > 1"   = "prop_LC3B_gt1",
  "ci_lower_skewness > 0"        = "ci_lower_LC3B_pos",
  "ci_lower_skewness > 0.5"      = "ci_lower_LC3B_gt0.5",
  "ci_lower_skewness > 1"        = "ci_lower_LC3B_gt1",
  "ci_upper_skewness > 0"        = "ci_upper_LC3B_pos",
  "ci_upper_skewness > 0.5"      = "ci_upper_LC3B_gt0.5",
  "ci_upper_skewness > 1"        = "ci_upper_LC3B_gt1"
)
colnames(prop_wide) <- ifelse(colnames(prop_wide) %in% names(rename_map),
                               rename_map[colnames(prop_wide)],
                               colnames(prop_wide))

per_patient <- medians_df %>%
  select(patient_id, group,
         LC3B_skewness_median, P62_skewness_median, n_cells.x) %>%
  rename(n_cells = n_cells.x) %>%
  left_join(prop_wide, by = "patient_id") %>%
  arrange(group, patient_id) %>%
  mutate(across(where(is.numeric) & !n_cells, ~ round(.x, 4)))

write.csv(per_patient, file.path(OUT_DATA, "skewness_stats_per_patient.csv"),
          row.names = FALSE)
message("Salvato: skewness_stats_per_patient.csv")
message("\n--- Riepilogo per paziente ---")
print(as.data.frame(per_patient[, c("patient_id", "group", "n_cells",
                                     "LC3B_skewness_median", "P62_skewness_median",
                                     "prop_LC3B_gt0.5")]))

# --------------------------------------------------------------------------
# TABELLA 2: Riepilogo per gruppo (media ± SD dei valori per paziente)
# --------------------------------------------------------------------------

metrics <- c("LC3B_skewness_median", "P62_skewness_median",
             "prop_LC3B_pos", "prop_LC3B_gt0.5", "prop_LC3B_gt1")

group_summary_rows <- lapply(metrics, function(m) {
  for_g3  <- per_patient[[m]][per_patient$group == "G3"]
  for_shh <- per_patient[[m]][per_patient$group == "SHH"]

  data.frame(
    metric    = m,
    G3_mean   = round(mean(for_g3,  na.rm = TRUE), 4),
    G3_sd     = round(sd(for_g3,    na.rm = TRUE), 4),
    G3_min    = round(min(for_g3,   na.rm = TRUE), 4),
    G3_max    = round(max(for_g3,   na.rm = TRUE), 4),
    SHH_mean  = round(mean(for_shh, na.rm = TRUE), 4),
    SHH_sd    = round(sd(for_shh,   na.rm = TRUE), 4),
    SHH_min   = round(min(for_shh,  na.rm = TRUE), 4),
    SHH_max   = round(max(for_shh,  na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
})
group_summary <- do.call(rbind, group_summary_rows)

write.csv(group_summary, file.path(OUT_DATA, "skewness_stats_group_summary.csv"),
          row.names = FALSE)
message("\nSalvato: skewness_stats_group_summary.csv")
message("\n--- Riepilogo per gruppo ---")
print(group_summary)

# --------------------------------------------------------------------------
# TABELLA 3: Test statistici G3 vs SHH (Wilcoxon rank-sum + effect size)
#
# Effect size rank-biserial correlation:
#   r_rb = 1 - (2 * W) / (n1 * n2)
#   Valori: < 0.3 piccolo, 0.3–0.5 medio, > 0.5 grande (in valore assoluto)
# --------------------------------------------------------------------------

wilcox_rows <- lapply(metrics, function(m) {
  g3  <- per_patient[[m]][per_patient$group == "G3"]
  shh <- per_patient[[m]][per_patient$group == "SHH"]
  n1  <- length(g3)
  n2  <- length(shh)

  wt  <- wilcox.test(g3, shh, exact = TRUE, alternative = "two.sided")
  W   <- as.numeric(wt$statistic)

  # Rank-biserial correlation (effect size)
  r_rb <- 1 - (2 * W) / (n1 * n2)

  # Direzione: positivo = G3 > SHH
  direction <- if (mean(g3) > mean(shh)) "G3 > SHH" else "SHH > G3"

  data.frame(
    metric       = m,
    W_statistic  = W,
    p_value      = round(wt$p.value, 4),
    r_rb         = round(r_rb, 3),
    direction    = direction,
    n_G3         = n1,
    n_SHH        = n2,
    stringsAsFactors = FALSE
  )
})
wilcox_df <- do.call(rbind, wilcox_rows)

write.csv(wilcox_df, file.path(OUT_DATA, "skewness_stats_wilcoxon.csv"),
          row.names = FALSE)
message("\nSalvato: skewness_stats_wilcoxon.csv")
message("\n--- Wilcoxon rank-sum tests G3 vs SHH ---")
print(wilcox_df)

message("\n--- Interpretazione ---")
for (i in seq_len(nrow(wilcox_df))) {
  r  <- wilcox_df[i, ]
  es <- abs(r$r_rb)
  es_label <- if (es > 0.5) "grande" else if (es > 0.3) "medio" else "piccolo"
  message(sprintf("  %s: W=%.0f, p=%.4f, r_rb=%.3f (%s), direzione: %s",
                  r$metric, r$W_statistic, r$p_value, r$r_rb, es_label, r$direction))
}

message("\n=== SKEWNESS STATISTICS completato ===")
