# =============================================================================
# mtor_nectin_lc3b_analysis.R
# Correlazione a tre variabili: pMTOR, Nectin2, LC3B
#
# Domanda biologica:
#   Quanto della correlazione pMTOR ~ Nectin2 e' indipendente da LC3B?
#   (cioe': la co-espressione pMTOR-Nectin2 e' mediata dallo stato autofagico?)
#
# Metodo: partial Spearman correlation
#   Per ciascun paziente, a livello single-cell:
#     1. Rank-transform le tre variabili
#     2. Residualizza rank(X) ~ rank(Z) con regressione lineare
#     3. Residualizza rank(Y) ~ rank(Z) con regressione lineare
#     4. Correlazione di Pearson sui residui = partial Spearman(X, Y | Z)
#     5. Significativita': z = atanh(r) * sqrt(n - 3 - k), k = 1 var controllo
#        p-value dalla distribuzione normale standard (valido per n >> 100)
#
# Confronto con correlazione semplice (da 06_spearman_correlations.csv).
#
# Input:  output/data/04_seurat_clustered.rds
#         output/data/06_spearman_correlations.csv
# Output: output/data/mtor_nectin_lc3b_partial_cor.csv
#         output/plots/PartialCor_01_simple_vs_partial.pdf
#         output/plots/PartialCor_02_per_patient.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

.sd <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

message("=== ANALISI CORRELAZIONE A TRE: pMTOR ~ Nectin2 ~ LC3B ===")

# --------------------------------------------------------------------------
# Funzione: partial Spearman correlation
#
# partial_spearman(x, y, z):
#   correlazione tra x e y al netto di z, su scala dei ranghi
#   Equivalente a ppcor::pcor.test(x, y, z, method="spearman")
#   ma senza dipendenza dal pacchetto ppcor.
#
# Returns: list(r, p_value, z_score, n)
# --------------------------------------------------------------------------

partial_spearman <- function(x, y, z) {
  n <- length(x)
  stopifnot(n == length(y), n == length(z), n > 3)

  rx <- rank(x)
  ry <- rank(y)
  rz <- rank(z)

  res_x <- residuals(lm(rx ~ rz))
  res_y <- residuals(lm(ry ~ rz))

  r   <- cor(res_x, res_y)  # partial Spearman coefficient
  # Fisher z-transformation con df = n - 3 - 1 (1 variabile di controllo)
  z_score <- atanh(r) * sqrt(n - 3 - 1)
  p_val   <- 2 * pnorm(-abs(z_score))

  list(r = r, p_value = p_val, z_score = z_score, n = n)
}

# --------------------------------------------------------------------------
# Carica Seurat object
# --------------------------------------------------------------------------

rds_path <- file.path(OUT_DATA, "04_seurat_clustered.rds")
if (!file.exists(rds_path)) stop("File non trovato: ", rds_path)

message("Caricamento Seurat object...")
seurat_obj <- readRDS(rds_path)
mat_all    <- GetAssayData(seurat_obj, layer = "data", assay = "MICS")

PROTEINS <- c("pMTOR", "Nectin2", "LC3B", "P62")
stopifnot(all(PROTEINS %in% rownames(mat_all)))
message(sprintf("  %d cellule, %d proteine disponibili",
                ncol(mat_all), nrow(mat_all)))

# --------------------------------------------------------------------------
# Triplette di correlazione da calcolare
#
# Per ciascuna coppia (X, Y), Z e' la terza variabile di controllo.
# Includiamo anche P62 come quarta variabile per completezza.
# --------------------------------------------------------------------------

TRIPLETS <- list(
  # Domanda principale: pMTOR ~ Nectin2 indipendente da LC3B?
  list(x = "pMTOR",  y = "Nectin2", z = "LC3B",
       label_simple  = "pMTOR ~ Nectin2",
       label_partial = "pMTOR ~ Nectin2 | LC3B"),

  # Domanda complementare: LC3B ~ Nectin2 indipendente da pMTOR?
  list(x = "LC3B",   y = "Nectin2", z = "pMTOR",
       label_simple  = "LC3B ~ Nectin2",
       label_partial = "LC3B ~ Nectin2 | pMTOR"),

  # pMTOR ~ LC3B al netto di Nectin2
  list(x = "pMTOR",  y = "LC3B",    z = "Nectin2",
       label_simple  = "pMTOR ~ LC3B",
       label_partial = "pMTOR ~ LC3B | Nectin2"),

  # pMTOR ~ Nectin2 al netto di P62
  list(x = "pMTOR",  y = "Nectin2", z = "P62",
       label_simple  = "pMTOR ~ Nectin2",
       label_partial = "pMTOR ~ Nectin2 | P62"),

  # P62 ~ Nectin2 al netto di pMTOR
  list(x = "P62",    y = "Nectin2", z = "pMTOR",
       label_simple  = "P62 ~ Nectin2",
       label_partial = "P62 ~ Nectin2 | pMTOR")
)

# --------------------------------------------------------------------------
# Calcolo per paziente
# --------------------------------------------------------------------------

message("\nCalcolo correlazioni per paziente...")

result_rows <- list()

for (pid in sort(unique(seurat_obj$patient_id))) {
  grp   <- unique(seurat_obj$group[seurat_obj$patient_id == pid])
  cells <- colnames(seurat_obj)[seurat_obj$patient_id == pid]
  n     <- length(cells)

  message(sprintf("  %s (%s): %d cellule", pid, grp, n))

  sub_mat <- as.matrix(mat_all[PROTEINS, cells])

  for (tr in TRIPLETS) {
    x_vals <- as.numeric(sub_mat[tr$x, ])
    y_vals <- as.numeric(sub_mat[tr$y, ])
    z_vals <- as.numeric(sub_mat[tr$z, ])

    # Correlazione semplice Spearman
    ct_simple <- cor.test(x_vals, y_vals, method = "spearman", exact = FALSE)

    # Correlazione parziale Spearman
    pc <- partial_spearman(x_vals, y_vals, z_vals)

    result_rows[[length(result_rows) + 1]] <- data.frame(
      patient_id     = pid,
      group          = grp,
      n_cells        = n,
      x_var          = tr$x,
      y_var          = tr$y,
      z_var          = tr$z,
      label_simple   = tr$label_simple,
      label_partial  = tr$label_partial,
      rho_simple     = round(ct_simple$estimate, 4),
      p_simple       = signif(ct_simple$p.value, 3),
      rho_partial    = round(pc$r, 4),
      p_partial      = signif(pc$p_value, 3),
      delta_rho      = round(pc$r - ct_simple$estimate, 4),
      stringsAsFactors = FALSE
    )
  }
}

result_df <- do.call(rbind, result_rows)
rownames(result_df) <- NULL

write.csv(result_df, file.path(OUT_DATA, "mtor_nectin_lc3b_partial_cor.csv"),
          row.names = FALSE)
message("\nSalvato: mtor_nectin_lc3b_partial_cor.csv")

# Stampa riassunto
message("\n--- Correlazioni semplici vs parziali (media per gruppo) ---")
summary_df <- result_df %>%
  group_by(label_partial, group) %>%
  summarise(
    rho_simple_mean  = round(mean(rho_simple), 3),
    rho_partial_mean = round(mean(rho_partial), 3),
    delta_mean       = round(mean(delta_rho), 3),
    .groups = "drop"
  )
print(as.data.frame(summary_df))

# --------------------------------------------------------------------------
# PLOT 1 — Confronto rho semplice vs parziale per paziente
# --------------------------------------------------------------------------

message("\nGenerazione plot...")

plot_df <- result_df %>%
  select(patient_id, group, label_partial, rho_simple, rho_partial) %>%
  pivot_longer(cols = c(rho_simple, rho_partial),
               names_to = "type", values_to = "rho") %>%
  mutate(
    type  = recode(type,
                   rho_simple  = "Simple Spearman",
                   rho_partial = "Partial Spearman"),
    type  = factor(type, levels = c("Simple Spearman", "Partial Spearman")),
    group = factor(group, levels = c("G3", "SHH")),
    patient_id = factor(patient_id,
                        levels = c(sort(unique(patient_id[group == "G3"])),
                                   sort(unique(patient_id[group == "SHH"]))))
  )

p1 <- ggplot(plot_df, aes(x = rho, y = patient_id,
                            shape = type, colour = group)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  geom_line(aes(group = paste(patient_id, label_partial)),
            colour = "grey70", linewidth = 0.5) +
  geom_point(size = 3.5, alpha = 0.9) +
  facet_wrap(~ label_partial, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = PALETTE_GROUP) +
  scale_shape_manual(values = c("Simple Spearman" = 16, "Partial Spearman" = 15)) +
  labs(
    title    = "Simple vs partial Spearman correlations per patient",
    subtitle = "Lines connect the same patient; partial = controlling for third variable",
    x        = "Spearman \u03c1",
    y        = NULL,
    colour   = "Group",
    shape    = "Correlation type"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text       = element_text(size = 9))

ggsave(file.path(OUT_PLOTS_CORR, "PartialCor_01_simple_vs_partial.pdf"),
       p1, width = 12, height = 10)
message("  Salvato: PartialCor_01_simple_vs_partial.pdf")

# --------------------------------------------------------------------------
# PLOT 2 — Delta rho (partial - simple) per triplet, separato per gruppo
# --------------------------------------------------------------------------

delta_df <- result_df %>%
  mutate(
    group = factor(group, levels = c("G3", "SHH")),
    label_partial = factor(label_partial)
  )

p2 <- ggplot(delta_df, aes(x = delta_rho, y = label_partial, colour = group)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  geom_jitter(height = 0.1, size = 3, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.4, colour = "black", linewidth = 0.6,
               position = position_dodge(width = 0.5)) +
  facet_wrap(~ group, ncol = 2) +
  scale_colour_manual(values = PALETTE_GROUP) +
  labs(
    title    = "\u0394\u03c1 = partial \u2212 simple Spearman",
    subtitle = "Negative = correlation reduced after partialling; bar = group mean",
    x        = "\u0394\u03c1 (partial \u2212 simple)",
    y        = NULL,
    colour   = "Group"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position  = "none",
        strip.background = element_rect(fill = "grey90"))

ggsave(file.path(OUT_PLOTS_CORR, "PartialCor_02_delta_rho.pdf"),
       p2, width = 10, height = 6)
message("  Salvato: PartialCor_02_delta_rho.pdf")

message("\n=== ANALISI CORRELAZIONE A TRE completata ===")
message(sprintf("Output in: %s", OUT_DATA))
message(sprintf("Plot in: %s", OUT_PLOTS_CORR))
