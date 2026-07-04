# =============================================================================
# figure_panel.R
# Genera la figura multi-pannello per il paper
#
# Pannello A: Heatmap correlazioni Spearman per paziente (scala colori ottimizzata)
# Pannello B: Distribuzione LC3B skewness per paziente (senza annotazioni testo)
#
# Input:  output/data/06_spearman_correlations.csv
#         output/data/skewness_cells.csv
# Output: output/plots/Figure_Panel_combined.pdf
#         output/plots/Figure_Panel_A_correlations.pdf
#         output/plots/Figure_Panel_B_skewness.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
})

.sd <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

message("=== FIGURE PANEL ===")

# --------------------------------------------------------------------------
# PANNELLO A — Correlazioni Spearman per paziente (scala colori migliorata)
# --------------------------------------------------------------------------

message("Pannello A: correlazioni per paziente...")

corr_df <- read.csv(file.path(OUT_DATA, "06_spearman_correlations.csv"),
                    stringsAsFactors = FALSE)

# Ordine pazienti: G3 in alto, SHH in basso
g3_pats  <- sort(unique(corr_df$patient_id[corr_df$group == "G3"]))
shh_pats <- sort(unique(corr_df$patient_id[corr_df$group == "SHH"]))
pat_order <- c(shh_pats, g3_pats)  # reverse: G3 in cima nel plot (y invertita)

# Etichetta con gruppo
corr_df <- corr_df %>%
  mutate(
    patient_label = paste0(patient_id, " (", group, ")"),
    patient_label = factor(patient_label,
      levels = paste0(pat_order, " (",
                      c(rep("SHH", length(shh_pats)), rep("G3", length(g3_pats))), ")")),
    # Riordina le pair: prima LC3B, poi P62, poi pMTOR; prima Nectin2 poi PVR
    pair = factor(pair, levels = c(
      "LC3B ~ Nectin2", "LC3B ~ PVR",
      "P62 ~ Nectin2",  "P62 ~ PVR",
      "pMTOR ~ Nectin2","pMTOR ~ PVR"
    ))
  )

# Scala colori: range teorico Spearman -1 a +1, bianco ESATTAMENTE a rho=0.
# I breakpoints sono compressi nel range negativo (pochi valori) e
# allargati nel range positivo (0-0.65) dove si concentrano i dati,
# cosi anche correlazioni moderate (es. 0.1-0.3) appaiono chiaramente rosse.
rho_min <- -1.0
rho_max <-  1.0

# Breakpoints in spazio rho: -1, -0.3, 0, 0.20, 0.40, 0.65, 1
# Posizioni normalizzate in [0,1]: (rho - rho_min) / (rho_max - rho_min)
break_rho <- c(-1.0, -0.3, 0.0, 0.20, 0.40, 0.65, 1.0)
vals_norm <- (break_rho - rho_min) / (rho_max - rho_min)
# = c(0.0, 0.35, 0.50, 0.60, 0.70, 0.825, 1.0)

cols_7 <- c("#1A5276",  # blu scuro  (rho = -1.0)
            "#5DADE2",  # blu medio  (rho = -0.3)
            "white",    # bianco     (rho =  0.0)  <-- ZERO ESATTO
            "#FADBD8",  # rosa tenue (rho =  0.2)
            "#E74C3C",  # rosso med  (rho =  0.4)
            "#922B21",  # rosso scuro(rho =  0.65)
            "#641E16")  # bordeaux   (rho =  1.0)

p_A <- ggplot(corr_df,
              aes(x = pair, y = patient_label, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3.2, fontface = "bold") +
  scale_fill_gradientn(
    colours = cols_7,
    values  = vals_norm,
    limits  = c(rho_min, rho_max),
    oob     = squish,   # valori fuori range (nessuno atteso) vengono compressi
    name    = "Spearman rho",
    guide   = guide_colourbar(barwidth = 8, barheight = 0.8,
                               title.position = "top", title.hjust = 0.5)
  ) +
  scale_x_discrete(expand = expansion(add = c(0.6, 0.5))) +
  labs(
    title = "Spearman correlations: autophagy markers vs DNAM-1 ligands",
    x     = NULL,
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1, size = 10),
    axis.text.y      = element_text(size = 10),
    panel.grid       = element_blank(),
    panel.border     = element_blank(),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    plot.title       = element_text(face = "bold", size = 12)
  )

ggsave(file.path(OUT_PLOTS_FIG, "Figure_Panel_A_correlations.pdf"),
       p_A, width = 10, height = 7)
message("  Salvato: Figure_Panel_A_correlations.pdf")

# --------------------------------------------------------------------------
# PANNELLO B — Distribuzione LC3B skewness per paziente (solo titolo, no subtitle)
# Rimosse le annotazioni testuali >0, >0.5, >1
# --------------------------------------------------------------------------

message("Pannello B: distribuzione skewness...")

df_all <- fread(file.path(OUT_DATA, "skewness_cells.csv"), data.table = FALSE)
df_all$group <- factor(df_all$group, levels = c("G3", "SHH"))

# Palette locale per pannello B: arancione (G3) e teal (SHH)
# Diversa da PALETTE_GROUP (rosso/blu) per non confondere con i colori dell'heatmap
PALETTE_SKEW <- c(G3 = "#E67E22", SHH = "#16A085")

# Ordine pazienti: G3 per primi
pat_order_sk <- c(g3_pats, shh_pats)
df_all$patient_id <- factor(df_all$patient_id, levels = pat_order_sk)

set.seed(SEED)
MAX_DENS <- 200000
df_dens  <- df_all[sample(nrow(df_all), min(MAX_DENS, nrow(df_all))), ]

p_B <- ggplot(df_dens, aes(x = LC3B_skewness, fill = group)) +
  geom_density(alpha = 0.55, colour = NA) +
  # Linee verticali per le soglie (senza testo)
  geom_vline(xintercept = 0,   linetype = "dashed",
             colour = "grey40",    linewidth = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             colour = "steelblue", linewidth = 0.5) +
  geom_vline(xintercept = 1,   linetype = "dashed",
             colour = "firebrick", linewidth = 0.5) +
  facet_wrap(~ patient_id, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = PALETTE_SKEW) +
  scale_x_continuous(limits = c(-2, 5), breaks = c(-1, 0, 1, 2, 3)) +
  labs(
    title = "LC3B cytoplasm skewness distribution per patient",
    x     = "LC3B cytoplasm skewness",
    y     = "Density",
    fill  = "Group"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position    = "top",
    strip.background   = element_rect(fill = "grey90"),
    strip.text         = element_text(size = 9),
    plot.title         = element_text(face = "bold", size = 12)
  )

ggsave(file.path(OUT_PLOTS_FIG, "Figure_Panel_B_skewness.pdf"),
       p_B, width = 14, height = 8)
message("  Salvato: Figure_Panel_B_skewness.pdf")

# --------------------------------------------------------------------------
# FIGURA COMBINATA
# --------------------------------------------------------------------------

message("Figura combinata...")

p_combined <- p_A / p_B +
  plot_layout(heights = c(1.1, 1.4)) +
  plot_annotation(tag_levels = "A",
                  tag_prefix = "",
                  theme = theme(plot.tag = element_text(face = "bold", size = 14)))

ggsave(file.path(OUT_PLOTS_FIG, "Figure_Panel_combined.pdf"),
       p_combined, width = 14, height = 16)
message("  Salvato: Figure_Panel_combined.pdf")

message("\n=== FIGURE PANEL completato ===")
message(sprintf("Output in: %s", OUT_PLOTS_FIG))
