# =============================================================================
# generate_figures_english.R
# Regenerate all reviewer-response figures with English labels
#
# Reads pre-computed CSV outputs (no raw data re-processing needed).
# Uses fread() for large cell-level files to avoid macOS errno-22 bug.
#
# Output overwrites the Italian-labelled PDFs in:
#   output/plots/skewness/
#   output/plots/correlations/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)        # percent_format (ggplot2 dependency, always available)
})
# ggridges is optional: Skewness_01 ridge plot is skipped gracefully if not installed
# Install with: install.packages("ggridges")

.sd <- tryCatch(dirname(sys.frame(1)$ofile),
                error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

message("=== GENERATING ENGLISH FIGURES ===")

# ── palette ──────────────────────────────────────────────────────────────────
# PALETTE_GROUP from config.R: c(G3 = ..., SHH = ...)

# ── cofactor ─────────────────────────────────────────────────────────────────
q_df     <- read.csv(file.path(OUT_DATA, "02_cofactor_quantiles.csv"))
COFACTOR <- median(q_df$p50)
message(sprintf("Cofactor: %.0f", COFACTOR))

# =============================================================================
# 1. CELL-LEVEL DATA (LC3B MFI vs Skewness) -fread mandatory
# =============================================================================

message("\nLoading lc3b_mfi_skewness_cells.csv with fread ...")
df_all <- fread(file.path(OUT_DATA, "lc3b_mfi_skewness_cells.csv"),
                data.table = FALSE)
df_all$group <- factor(df_all$group, levels = c("G3", "SHH"))
message(sprintf("  %s cells loaded", format(nrow(df_all), big.mark = ",")))

set.seed(SEED)
df_plot <- df_all %>%
  group_by(patient_id) %>%
  slice_sample(n = 5000) %>%
  ungroup()

mfi_median <- median(df_all$LC3B_MFI_arcsinh)

# ── Figure LC3B_01: scatter all cells, coloured by group ─────────────────────
message("LC3B_01 ...")
p1 <- ggplot(df_plot,
             aes(x = LC3B_MFI_arcsinh, y = LC3B_skew, colour = group)) +
  geom_point(size = 0.15, alpha = 0.25) +
  geom_vline(xintercept = mfi_median, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_density_2d(linewidth = 0.5, alpha = 0.9, bins = 8) +
  scale_colour_manual(values = PALETTE_GROUP) +
  scale_y_continuous(limits = c(-3, 5)) +
  scale_x_continuous(limits = c(0, NA)) +
  annotate("text", x = -Inf, y = 5, hjust = -0.1, vjust = 1,
           label = "SPARSE PUNCTATE\n(few bright spots)", size = 2.8,
           colour = "grey30", fontface = "italic") +
  annotate("text", x = mfi_median * 1.05, y = 5, hjust = 0, vjust = 1,
           label = "DENSE PUNCTATE\n(many spots -> lower skewness)", size = 2.8,
           colour = "grey30", fontface = "italic") +
  annotate("text", x = mfi_median * 1.05, y = -2.8, hjust = 0, vjust = 0,
           label = "UNIFORMLY BRIGHT\n(high MFI, not punctate)", size = 2.8,
           colour = "grey30", fontface = "italic") +
  annotate("text", x = -Inf, y = -2.8, hjust = -0.1, vjust = 0,
           label = "DIFFUSE / BACKGROUND\n(low MFI, low skewness)", size = 2.8,
           colour = "grey30", fontface = "italic") +
  labs(
    title    = "LC3B: Mean Fluorescence Intensity (arcsinh) vs Cytoplasmic Skewness",
    subtitle = sprintf(
      "Vertical line = overall MFI median (%.2f) | Horizontal line = skewness 0\n5,000-cell subsample per patient",
      mfi_median),
    x      = sprintf("LC3B MFI [arcsinh(MFI / %d)]", COFACTOR),
    y      = "LC3B Cytoplasm Intensity Skewness",
    colour = "Group"
  ) +
  theme_bw(base_size = 11) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(OUT_PLOTS_SKEW, "LC3B_01_MFI_vs_Skewness_byGroup.pdf"),
       p1, width = 10, height = 7)
message("  Saved.")

# ── Figure LC3B_02: per-patient facet ────────────────────────────────────────
message("LC3B_02 ...")
p2 <- ggplot(df_plot,
             aes(x = LC3B_MFI_arcsinh, y = LC3B_skew, colour = group)) +
  geom_point(size = 0.2, alpha = 0.2) +
  geom_vline(xintercept = mfi_median, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, colour = "black",
              formula = y ~ x) +
  facet_wrap(~ patient_id, ncol = 3, scales = "free_x") +
  scale_colour_manual(values = PALETTE_GROUP) +
  scale_y_continuous(limits = c(-3, 5)) +
  labs(
    title    = "LC3B MFI vs Skewness -per patient",
    subtitle = "Black line = linear regression per patient | Dashed lines = overall MFI median and skewness 0",
    x        = sprintf("LC3B MFI [arcsinh / %d]", COFACTOR),
    y        = "LC3B Skewness",
    colour   = "Group"
  ) +
  theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "grey90"),
        legend.position  = "bottom")

ggsave(file.path(OUT_PLOTS_SKEW, "LC3B_02_MFI_vs_Skewness_perPatient.pdf"),
       p2, width = 12, height = 14)
message("  Saved.")

# ── Figure LC3B_03: hexbin density G3 vs SHH ────────────────────────────────
message("LC3B_03 ...")
p3_g3 <- ggplot(df_plot[df_plot$group == "G3", ],
                aes(x = LC3B_MFI_arcsinh, y = LC3B_skew)) +
  geom_hex(bins = 60) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "white", linewidth = 0.5) +
  geom_vline(xintercept = mfi_median, linetype = "dashed", colour = "white", linewidth = 0.5) +
  scale_fill_viridis_c(option = "inferno", name = "Cell count") +
  scale_y_continuous(limits = c(-3, 5)) +
  labs(title = "G3 -LC3B MFI vs Skewness (cell density)",
       x = "LC3B MFI arcsinh", y = "LC3B Skewness") +
  theme_bw(base_size = 11)

p3_shh <- ggplot(df_plot[df_plot$group == "SHH", ],
                 aes(x = LC3B_MFI_arcsinh, y = LC3B_skew)) +
  geom_hex(bins = 60) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "white", linewidth = 0.5) +
  geom_vline(xintercept = mfi_median, linetype = "dashed", colour = "white", linewidth = 0.5) +
  scale_fill_viridis_c(option = "viridis", name = "Cell count") +
  scale_y_continuous(limits = c(-3, 5)) +
  labs(title = "SHH -LC3B MFI vs Skewness (cell density)",
       x = "LC3B MFI arcsinh", y = "LC3B Skewness") +
  theme_bw(base_size = 11)

ggsave(file.path(OUT_PLOTS_SKEW, "LC3B_03_density_hexbin_G3vsSHH.pdf"),
       p3_g3 + p3_shh, width = 14, height = 6)
message("  Saved.")

# ── Figure LC3B_04: patient medians scatter ──────────────────────────────────
message("LC3B_04 ...")
pat_summary <- df_all %>%
  group_by(group, patient_id) %>%
  summarise(
    LC3B_MFI_median  = median(LC3B_MFI_arcsinh, na.rm = TRUE),
    LC3B_skew_median = median(LC3B_skew,         na.rm = TRUE),
    n_cells          = n(),
    .groups = "drop"
  )

p4 <- ggplot(pat_summary,
             aes(x = LC3B_MFI_median, y = LC3B_skew_median,
                 colour = group, label = patient_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = median(pat_summary$LC3B_MFI_median),
             linetype = "dashed", colour = "grey50") +
  geom_point(size = 5, alpha = 0.9) +
  geom_text(size = 3, vjust = -0.9, show.legend = FALSE) +
  scale_colour_manual(values = PALETTE_GROUP) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
           label = "Sparse punctate\n(low MFI, high skewness)", size = 3,
           colour = "grey40", fontface = "italic") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3,
           label = "Dense punctate\n(high MFI, moderate skewness)", size = 3,
           colour = "grey40", fontface = "italic") +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.05, vjust = -0.3,
           label = "Diffuse\n(low MFI, low skewness)", size = 3,
           colour = "grey40", fontface = "italic") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.3,
           label = "Uniformly bright\n(high MFI, negative skewness)", size = 3,
           colour = "grey40", fontface = "italic") +
  labs(
    title    = "Patient medians: LC3B MFI vs LC3B Skewness",
    subtitle = "Each point = one patient | Quadrants = biological interpretation",
    x        = sprintf("Median LC3B MFI [arcsinh / %d]", COFACTOR),
    y        = "Median LC3B Cytoplasmic Skewness",
    colour   = "Group"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_PLOTS_SKEW, "LC3B_04_patient_medians_MFI_vs_Skew.pdf"),
       p4, width = 9, height = 7)
message("  Saved.")

# =============================================================================
# 2. PER-PATIENT SUMMARY DATA (skewness CSVs)
# =============================================================================

stats_pp <- read.csv(file.path(OUT_DATA, "skewness_stats_per_patient.csv"),
                     stringsAsFactors = FALSE)
stats_pp$group <- factor(stats_pp$group, levels = c("G3", "SHH"))

prop_df <- read.csv(file.path(OUT_DATA, "skewness_punctate_proportion.csv"),
                    stringsAsFactors = FALSE)
prop_df$group <- factor(prop_df$group, levels = c("G3", "SHH"))

# ── Figure LC3B_05: punctate proportions at multiple thresholds ───────────────
message("LC3B_05 ...")
prop_long <- stats_pp %>%
  select(patient_id, group,
         `skewness > 0`   = prop_LC3B_pos,
         `skewness > 0.5` = prop_LC3B_gt0.5,
         `skewness > 1.0` = prop_LC3B_gt1) %>%
  pivot_longer(cols = c(`skewness > 0`, `skewness > 0.5`, `skewness > 1.0`),
               names_to = "threshold", values_to = "proportion") %>%
  mutate(threshold = factor(threshold,
                            levels = c("skewness > 0", "skewness > 0.5", "skewness > 1.0")))

wilcox_by_thresh <- prop_long %>%
  group_by(threshold) %>%
  summarise(
    p = wilcox.test(proportion[group == "G3"],
                    proportion[group == "SHH"], exact = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = ifelse(p < 0.05,
                        sprintf("p = %.3f *", p),
                        sprintf("p = %.3f", p)),
         y_pos = max(prop_long$proportion) * 1.05)

p5 <- ggplot(prop_long, aes(x = group, y = proportion, colour = group)) +
  geom_jitter(width = 0.12, size = 3, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.35, colour = "black", linewidth = 0.6) +
  geom_text(data = wilcox_by_thresh,
            aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3.2, colour = "grey20") +
  facet_wrap(~ threshold, ncol = 3) +
  scale_colour_manual(values = PALETTE_GROUP) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Proportion of cells with punctate LC3B (skewness above threshold)",
    subtitle = "Each point = one patient median | Bar = group mean | Wilcoxon G3 vs SHH",
    x        = NULL,
    y        = "Proportion of cells",
    colour   = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey90"))

ggsave(file.path(OUT_PLOTS_SKEW, "LC3B_05_punctate_proportion_thresholds.pdf"),
       p5, width = 10, height = 5)
message("  Saved.")

# ── Figure LC3B_06: cumulative skewness distribution per patient ──────────────
message("LC3B_06 ...")
# Use subsampled df_plot for visual clarity
p6 <- ggplot(df_plot, aes(x = LC3B_skew, colour = group,
                           group = patient_id, linetype = group)) +
  stat_ecdf(linewidth = 0.8, alpha = 0.85) +
  geom_vline(xintercept = c(0, 0.5, 1.0), linetype = "dotted",
             colour = "grey50", linewidth = 0.4) +
  annotate("text", x = 0,   y = 0.05, hjust = -0.1, size = 2.8,
           label = "0", colour = "grey40") +
  annotate("text", x = 0.5, y = 0.05, hjust = -0.1, size = 2.8,
           label = "0.5", colour = "grey40") +
  annotate("text", x = 1.0, y = 0.05, hjust = -0.1, size = 2.8,
           label = "1.0", colour = "grey40") +
  scale_colour_manual(values = PALETTE_GROUP) +
  scale_x_continuous(limits = c(-2, 4)) +
  labs(
    title    = "Cumulative distribution of LC3B cytoplasmic skewness",
    subtitle = "Each line = one patient (5,000-cell subsample) | Dotted verticals = thresholds 0, 0.5, 1.0",
    x        = "LC3B Cytoplasm Intensity Skewness",
    y        = "Cumulative proportion of cells",
    colour   = "Group",
    linetype = "Group"
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(OUT_PLOTS_SKEW, "LC3B_06_cumulative_distribution.pdf"),
       p6, width = 9, height = 6)
message("  Saved.")

# =============================================================================
# 3. SKEWNESS SUMMARY PLOTS
# =============================================================================

skew_pm <- read.csv(file.path(OUT_DATA, "skewness_patient_medians.csv"),
                    stringsAsFactors = FALSE)
skew_pm$group <- factor(skew_pm$group, levels = c("G3", "SHH"))

wt_skew <- wilcox.test(
  stats_pp$LC3B_skewness_median[stats_pp$group == "G3"],
  stats_pp$LC3B_skewness_median[stats_pp$group == "SHH"],
  exact = TRUE)

# ── Skewness_01: ridge plot of skewness per patient (requires ggridges) ──────
message("Skewness_01 ...")
tryCatch({
  p_sk1 <- ggplot(df_plot, aes(x = LC3B_skew, y = patient_id, fill = group)) +
    ggridges::geom_density_ridges(alpha = 0.75, scale = 1.2,
                                   rel_min_height = 0.01) +
    geom_vline(xintercept = c(0, 0.5, 1.0), linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = PALETTE_GROUP) +
    scale_x_continuous(limits = c(-3, 5)) +
    labs(
      title    = "LC3B cytoplasmic skewness distribution - per patient",
      subtitle = "Dashed lines: thresholds 0, 0.5, 1.0 | 5,000-cell subsample per patient",
      x        = "LC3B Cytoplasm Intensity Skewness",
      y        = NULL,
      fill     = "Group"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(OUT_PLOTS_SKEW, "Skewness_01_ridge_per_patient.pdf"),
         p_sk1, width = 9, height = 8)
  message("  Saved.")
}, error = function(e) {
  message("  [WARN] ggridges not available, skipping Skewness_01. Install with: install.packages('ggridges')")
})

# ── Skewness_02: violin G3 vs SHH ────────────────────────────────────────────
message("Skewness_02 ...")
p_sk2 <- ggplot(df_plot, aes(x = group, y = LC3B_skew,
                              fill = group, colour = group)) +
  geom_violin(alpha = 0.5, trim = TRUE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(data = stats_pp,
              aes(x = group, y = LC3B_skewness_median),
              width = 0.08, size = 3.5, shape = 21,
              fill = "white", colour = "black", stroke = 1.2,
              inherit.aes = FALSE) +
  scale_fill_manual(values = PALETTE_GROUP) +
  scale_colour_manual(values = PALETTE_GROUP) +
  scale_y_continuous(limits = c(-3, 5)) +
  annotate("text", x = 1.5, y = 4.7,
           label = sprintf("Wilcoxon p = %.3f", wt_skew$p.value),
           size = 4, colour = "grey20") +
  labs(
    title    = "LC3B cytoplasmic skewness: G3 vs SHH",
    subtitle = "Violin = cell distribution (5,000/patient subsample) | Open circles = patient medians\nWilcoxon test on patient medians (n=4 G3, n=4 SHH)",
    x        = NULL,
    y        = "LC3B Cytoplasm Intensity Skewness",
    fill     = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(OUT_PLOTS_SKEW, "Skewness_02_violin_G3vsSHH.pdf"),
       p_sk2, width = 6, height = 7)
message("  Saved.")

# ── Skewness_03: LC3B vs P62 skewness scatter ─────────────────────────────────
message("Skewness_03 ...")
p_sk3 <- ggplot(stats_pp, aes(x = LC3B_skewness_median, y = P62_skewness_median,
                               colour = group, label = patient_id)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text(vjust = -0.9, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = PALETTE_GROUP) +
  labs(
    title    = "Autophagy marker skewness: LC3B vs P62 (patient medians)",
    subtitle = "Each point = one patient | Skewness of cytoplasmic intensity distribution",
    x        = "Median LC3B cytoplasmic skewness",
    y        = "Median P62 cytoplasmic skewness",
    colour   = "Group"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_PLOTS_SKEW, "Skewness_03_scatter_LC3B_vs_P62.pdf"),
       p_sk3, width = 7, height = 6)
message("  Saved.")

# ── Skewness_04: patient medians barplot ─────────────────────────────────────
message("Skewness_04 ...")
sk4_long <- stats_pp %>%
  select(patient_id, group, LC3B = LC3B_skewness_median,
         P62 = P62_skewness_median) %>%
  pivot_longer(cols = c(LC3B, P62), names_to = "marker",
               values_to = "skewness_median")

p_sk4 <- ggplot(sk4_long,
                aes(x = patient_id, y = skewness_median,
                    fill = group)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  facet_wrap(~ marker, ncol = 1) +
  scale_fill_manual(values = PALETTE_GROUP) +
  labs(
    title    = "Median cytoplasmic skewness per patient",
    subtitle = "LC3B and P62 | Higher skewness = more punctate signal",
    x        = NULL,
    y        = "Median skewness",
    fill     = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        strip.background = element_rect(fill = "grey90"))

ggsave(file.path(OUT_PLOTS_SKEW, "Skewness_04_patient_medians.pdf"),
       p_sk4, width = 9, height = 7)
message("  Saved.")

# ── Skewness_05: punctate proportion LC3B ────────────────────────────────────
message("Skewness_05 ...")
p_sk5 <- ggplot(stats_pp,
                aes(x = group, y = prop_LC3B_gt0.5,
                    colour = group)) +
  geom_jitter(width = 0.12, size = 4, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.3, colour = "black", linewidth = 0.6) +
  scale_colour_manual(values = PALETTE_GROUP) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, NA)) +
  annotate("text", x = 1.5,
           y = max(stats_pp$prop_LC3B_gt0.5) * 1.05,
           label = sprintf("Wilcoxon p = %.3f",
                           wilcox.test(
                             stats_pp$prop_LC3B_gt0.5[stats_pp$group == "G3"],
                             stats_pp$prop_LC3B_gt0.5[stats_pp$group == "SHH"],
                             exact = TRUE)$p.value),
           size = 3.5, colour = "grey20") +
  labs(
    title    = "Proportion of cells with punctate LC3B (skewness > 0.5)",
    subtitle = "Each point = one patient | Bar = group mean\nWilcoxon test on patient proportions",
    x        = NULL,
    y        = "Proportion of cells with skewness > 0.5",
    colour   = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(OUT_PLOTS_SKEW, "Skewness_05_punctate_proportion_LC3B.pdf"),
       p_sk5, width = 5, height = 6)
message("  Saved.")

# ── Skewness_06: punctate proportion at multiple thresholds (same as LC3B_05)─
message("Skewness_06 ...")
# Already done as LC3B_05; copy to Skewness_06 path as well
ggsave(file.path(OUT_PLOTS_SKEW, "Skewness_06_punctate_proportion_thresholds.pdf"),
       p5, width = 10, height = 5)
message("  Saved.")

# =============================================================================
# 4. KI67 SENSITIVITY ANALYSIS FIGURES
# =============================================================================

message("\nKi67 sensitivity figures ...")
ki67_df  <- read.csv(file.path(OUT_DATA, "ki67_sensitivity_correlations.csv"),
                     stringsAsFactors = FALSE)
orig_corr <- read.csv(file.path(OUT_DATA, "06_spearman_correlations.csv"),
                      stringsAsFactors = FALSE)

orig_sub <- orig_corr %>%
  select(patient_id, group, pair, rho) %>%
  rename(rho_all = rho)

ki67_sub <- ki67_df %>%
  select(patient_id, pair, rho, n_cells) %>%
  rename(rho_ki67low = rho, n_ki67low = n_cells)

comparison_df <- orig_sub %>%
  left_join(ki67_sub, by = c("patient_id", "pair")) %>%
  mutate(
    delta_rho = round(rho_ki67low - rho_all, 4),
    group     = factor(group, levels = c("G3", "SHH"))
  )

plot_comp <- comparison_df %>%
  select(patient_id, group, pair, rho_all, rho_ki67low) %>%
  pivot_longer(cols = c(rho_all, rho_ki67low),
               names_to = "subset", values_to = "rho") %>%
  mutate(
    subset = recode(subset,
                    rho_all     = "All cells",
                    rho_ki67low = "Ki67-low (excl. top 25% proliferating)"),
    subset = factor(subset,
                    levels = c("All cells", "Ki67-low (excl. top 25% proliferating)"))
  )

p_ki67_1 <- ggplot(plot_comp,
                   aes(x = rho, y = patient_id, colour = group, shape = subset)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  geom_line(aes(group = paste(patient_id, pair)),
            colour = "grey70", linewidth = 0.5) +
  geom_point(size = 3, alpha = 0.9) +
  facet_wrap(~ pair, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = PALETTE_GROUP) +
  scale_shape_manual(values = c("All cells" = 16,
                                "Ki67-low (excl. top 25% proliferating)" = 15)) +
  labs(
    title    = "Sensitivity analysis: Spearman correlations with and without Ki67-high cells",
    subtitle = "Lines connect the same patient; squares = Ki67-low subset (excl. top 25% proliferating cells)",
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
       p_ki67_1, width = 12, height = 10)
message("  Saved: Ki67_01")

corr_mean_ki67low <- ki67_df %>%
  group_by(group, autophagy, ligand) %>%
  summarise(rho_mean = round(mean(rho), 3), .groups = "drop") %>%
  mutate(group = factor(group, levels = c("G3", "SHH")))

p_ki67_2 <- ggplot(corr_mean_ki67low,
                   aes(x = ligand, y = autophagy, fill = rho_mean)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = rho_mean), size = 4) +
  facet_wrap(~ group) +
  scale_fill_gradient2(low = "#1A5276", mid = "white", high = "#C0392B",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman rho") +
  labs(
    title    = "Spearman rho -- Ki67-low cells only (excluding top 25% proliferating)",
    subtitle = "Mean rho across patients per group",
    x        = "DNAM-1 ligand",
    y        = "Autophagy marker"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(OUT_PLOTS_CORR, "Ki67_02_rho_heatmap_ki67low.pdf"),
       p_ki67_2, width = 8, height = 5)
message("  Saved: Ki67_02")

message("\n=== ALL ENGLISH FIGURES GENERATED ===")
message(sprintf("Skewness plots: %s", OUT_PLOTS_SKEW))
message(sprintf("Correlation plots: %s", OUT_PLOTS_CORR))
