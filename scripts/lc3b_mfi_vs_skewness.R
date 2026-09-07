# =============================================================================
# lc3b_mfi_vs_skewness.R
# Analisi diagnostica: LC3B MFI vs LC3B Skewness per cellula
#
# OBIETTIVO (risposta al revisore — Punto 1C):
#   Verificare empiricamente se G3 occupa il quadrante
#   "alta MFI + bassa/moderata skewness" (dense punctate, autofagia intensa)
#   oppure "bassa MFI + bassa skewness" (diffuso, nessuna autofagia).
#
#   Se G3 → alta MFI + skewness bassa/moderata  → interpretazione paper CORRETTA
#   Se G3 → bassa MFI + skewness bassa           → conclusione biologica DA RIVEDERE
#
# METODO:
#   Per ogni ROI, i file *Skewness*.csv e *_WO.csv (MFI) sono
#   row-aligned (stessa cellula alla stessa riga, verificato empiricamente).
#   Si leggono entrambi e si aggiunge LC3B_MFI al dataset skewness.
#
# Input:  data/**/ROI*/ (file Skewness CSV + file MFI CSV)
# Output: output/data/lc3b_mfi_skewness_cells.csv
#         output/plots/skewness/LC3B_MFI_vs_Skewness_*.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

.sd <- tryCatch(dirname(sys.frame(1)$ofile),
                error = function(e) Sys.getenv("MACSIMA_SCRIPTS_DIR"))
source(file.path(.sd, "config.R"), local = TRUE)

message("=== LC3B MFI vs SKEWNESS — Analisi diagnostica per revisore ===")

# --------------------------------------------------------------------------
# Cofactor: letto da CSV generato dal pipeline principale
# --------------------------------------------------------------------------

cofactor_csv <- file.path(OUT_DATA, "02_cofactor_quantiles.csv")
if (file.exists(cofactor_csv)) {
  q_df     <- read.csv(cofactor_csv)
  COFACTOR <- median(q_df$p50)
  message(sprintf("Cofactor riletto da 02_cofactor_quantiles.csv: %.0f", COFACTOR))
} else {
  COFACTOR <- 282
  message(sprintf("02_cofactor_quantiles.csv non trovato — uso cofactor default = %d", COFACTOR))
}

# --------------------------------------------------------------------------
# Nomi colonne attesi nei due file
# --------------------------------------------------------------------------

COL_LC3B_SKEW <- "LC3B Cytoplasm Intensity Skewness"
COL_LC3B_MFI  <- "LC3B Cyto Exp"          # nel file *_WO.csv
COL_P62_SKEW  <- "P62 Cytoplasm Intensity Skewness"
COL_P62_MFI   <- "P62 Cyto Exp"
COL_KI67_MFI  <- "Ki_67 REA1123 Biomarker Exp"  # presente in entrambi i file

# --------------------------------------------------------------------------
# Funzioni helper
# --------------------------------------------------------------------------

# read_csv_header: legge la prima riga di un CSV e restituisce i nomi delle colonne
# Usa readLines (non fread) per evitare errno 22 su macOS con file grandi
read_csv_header <- function(f) {
  tryCatch({
    line1 <- readLines(f, n = 1L, warn = FALSE, encoding = "UTF-8")
    # read.csv con header=TRUE e nrows=0 legge solo la riga header → colnames corretti
    con <- textConnection(line1)
    on.exit(close(con))
    colnames(read.csv(con, nrows = 0, check.names = FALSE))
  }, error = function(e) character(0))
}

find_skewness_file <- function(roi_dir) {
  files <- list.files(roi_dir, pattern = "Sk?ewness", full.names = TRUE,
                      ignore.case = TRUE)
  for (f in files) {
    hdr <- read_csv_header(f)
    if (COL_LC3B_SKEW %in% hdr) return(f)
  }
  NULL
}

find_mfi_file <- function(roi_dir) {
  all_csv <- list.files(roi_dir, pattern = "\\.csv$", full.names = TRUE,
                        ignore.case = FALSE)
  # Seleziona SOLO file che contengono la colonna LC3B Cyto Exp (non Z-Score)
  for (f in all_csv) {
    hdr <- read_csv_header(f)
    if (COL_LC3B_MFI %in% hdr) return(f)
  }
  NULL
}

# --------------------------------------------------------------------------
# Caricamento: per ogni paziente × ROI, leggi entrambi i file e unisci
# --------------------------------------------------------------------------

patient_info <- unique(do.call(rbind, lapply(SAMPLES, function(s) {
  data.frame(patient_id = s$patient_id, group = s$group, stringsAsFactors = FALSE)
})))

df_list <- list()

for (i in seq_len(nrow(patient_info))) {
  pid <- patient_info$patient_id[i]
  grp <- patient_info$group[i]

  pdir <- file.path(DATA_DIR, pid)
  if (!dir.exists(pdir)) next

  roi_dirs <- list.dirs(pdir, recursive = FALSE, full.names = TRUE)
  message(sprintf("\n%s (%s):", pid, grp))

  for (rdir in roi_dirs) {
    roi_name  <- basename(rdir)
    skew_file <- find_skewness_file(rdir)
    mfi_file  <- find_mfi_file(rdir)

    if (is.null(skew_file) || is.null(mfi_file)) {
      message(sprintf("  [SKIP] %s: skew=%s mfi=%s",
                      roi_name, !is.null(skew_file), !is.null(mfi_file)))
      next
    }

    # Rileva dinamicamente il nome della colonna P62 nel file MFI
    # (varia tra ROI: "P62 Cyto Exp", "P62 Cell Intensity Sum", ecc.)
    mfi_hdr <- read_csv_header(mfi_file)
    p62_mfi_col <- mfi_hdr[grepl("p62", mfi_hdr, ignore.case = TRUE) &
                             !grepl("z.score|skewness", mfi_hdr, ignore.case = TRUE)]
    if (length(p62_mfi_col) == 0) {
      message(sprintf("  [WARN] %s: nessuna colonna P62 trovata nel file MFI — LC3B-only", roi_name))
      p62_mfi_col <- NULL
    } else {
      p62_mfi_col <- p62_mfi_col[1]
    }

    # Leggi colonne necessarie
    cols_to_read <- c(COL_LC3B_MFI, p62_mfi_col)
    dt_skew <- tryCatch(
      fread(skew_file, select = c(COL_KI67_MFI, COL_LC3B_SKEW, COL_P62_SKEW),
            data.table = FALSE),
      error = function(e) NULL
    )
    dt_mfi <- tryCatch(
      fread(mfi_file, select = cols_to_read, data.table = FALSE),
      error = function(e) NULL
    )

    if (is.null(dt_skew) || is.null(dt_mfi)) {
      message(sprintf("  [ERR] %s: lettura fallita", roi_name))
      next
    }

    if (nrow(dt_skew) != nrow(dt_mfi)) {
      message(sprintf("  [ERR] %s: riga mismatch skew=%d mfi=%d",
                      roi_name, nrow(dt_skew), nrow(dt_mfi)))
      next
    }

    # Row-bind: i file sono garantiti row-aligned (verificato empiricamente)
    dt <- data.frame(
      patient_id   = pid,
      group        = grp,
      roi_id       = paste0(pid, "_", roi_name),
      Ki67_MFI     = as.numeric(dt_skew[[COL_KI67_MFI]]),
      LC3B_MFI     = as.numeric(dt_mfi[[COL_LC3B_MFI]]),
      LC3B_skew    = as.numeric(dt_skew[[COL_LC3B_SKEW]]),
      P62_MFI      = if (!is.null(p62_mfi_col)) as.numeric(dt_mfi[[p62_mfi_col]]) else NA_real_,
      P62_skew     = as.numeric(dt_skew[[COL_P62_SKEW]]),
      stringsAsFactors = FALSE
    )

    # Rimuovi righe con NA
    dt <- dt[complete.cases(dt), ]
    if (nrow(dt) == 0) next

    df_list[[paste0(pid, "_", roi_name)]] <- dt
    message(sprintf("  %s -> OK (%d cellule)", roi_name, nrow(dt)))
  }
}

if (length(df_list) == 0) stop("Nessun dato caricato.")

df_all <- do.call(rbind, df_list)
rownames(df_all) <- NULL
df_all$group <- factor(df_all$group, levels = c("G3", "SHH"))

# Arcsinh trasformazione su MFI (stesso cofactor della pipeline principale)
df_all$LC3B_MFI_arcsinh <- asinh(df_all$LC3B_MFI / COFACTOR)
df_all$P62_MFI_arcsinh  <- asinh(df_all$P62_MFI  / COFACTOR)

message(sprintf("\nTotale cellule caricate: %s da %d ROI",
                format(nrow(df_all), big.mark = ","), length(df_list)))
message(sprintf("G3:  %s cellule", format(sum(df_all$group == "G3"),  big.mark = ",")))
message(sprintf("SHH: %s cellule", format(sum(df_all$group == "SHH"), big.mark = ",")))

write.csv(df_all, file.path(OUT_DATA, "lc3b_mfi_skewness_cells.csv"),
          row.names = FALSE, quote = FALSE)
message("Salvato: lc3b_mfi_skewness_cells.csv")

# --------------------------------------------------------------------------
# Subsample per i plot (max 5.000 cellule per paziente, riproducibile)
# --------------------------------------------------------------------------

set.seed(SEED)
df_plot <- df_all %>%
  group_by(patient_id) %>%
  slice_sample(n = 5000) %>%
  ungroup()

message(sprintf("Subsample per plot: %s cellule", format(nrow(df_plot), big.mark = ",")))

# --------------------------------------------------------------------------
# PLOT 1 — LC3B MFI (arcsinh) vs LC3B Skewness — colorato per gruppo
#
# CHIAVE INTERPRETATIVA:
#   Quadrante alto-MFI + bassa/mod. skewness → denso puntato  (G3 claim)
#   Quadrante bassa-MFI + bassa skewness     → diffuso        (alternativa)
#   Quadrante bassa-MFI + alta skewness      → sparso puntato (SHH claim)
# --------------------------------------------------------------------------

# Calcola percentili per disegnare i quadranti
mfi_median  <- median(df_all$LC3B_MFI_arcsinh)
skew_zero   <- 0   # soglia biologica naturale (skewness = 0)

p1 <- ggplot(df_plot,
             aes(x = LC3B_MFI_arcsinh, y = LC3B_skew, colour = group)) +
  geom_point(size = 0.15, alpha = 0.25) +
  geom_vline(xintercept = mfi_median, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = skew_zero, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_density_2d(linewidth = 0.5, alpha = 0.9, bins = 8) +
  scale_colour_manual(values = PALETTE_GROUP) +
  scale_y_continuous(limits = c(-3, 5)) +
  scale_x_continuous(limits = c(0, NA)) +

  # Quadrant annotations
  annotate("text", x = -Inf, y = 5,   hjust = -0.1, vjust = 1,
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
      "Vertical line = overall MFI median (%.2f) | Horizontal line = skewness 0 | 5,000-cell subsample per patient",
      mfi_median),
    x        = sprintf("LC3B MFI [arcsinh(MFI / %d)]", COFACTOR),
    y        = "LC3B Cytoplasm Intensity Skewness",
    colour   = "Group"
  ) +
  theme_bw(base_size = 11) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(OUT_PLOTS_SKEW, "LC3B_01_MFI_vs_Skewness_byGroup.pdf"),
       p1, width = 10, height = 7)
message("Salvato: LC3B_01_MFI_vs_Skewness_byGroup.pdf")

# --------------------------------------------------------------------------
# PLOT 2 — Stesso scatter, separato per paziente (facet)
# --------------------------------------------------------------------------

p2 <- ggplot(df_plot,
             aes(x = LC3B_MFI_arcsinh, y = LC3B_skew, colour = group)) +
  geom_point(size = 0.2, alpha = 0.2) +
  geom_vline(xintercept = mfi_median, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, colour = "black") +
  facet_wrap(~ patient_id, ncol = 3, scales = "free_x") +
  scale_colour_manual(values = PALETTE_GROUP) +
  scale_y_continuous(limits = c(-3, 5)) +
  labs(
    title    = "LC3B MFI vs Skewness - per patient",
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
message("Salvato: LC3B_02_MFI_vs_Skewness_perPatient.pdf")

# --------------------------------------------------------------------------
# PLOT 3 — Distribuzione 2D come density heatmap (hexbin) per G3 e SHH
# --------------------------------------------------------------------------

p3_g3  <- ggplot(df_plot[df_plot$group == "G3", ],
                 aes(x = LC3B_MFI_arcsinh, y = LC3B_skew)) +
  geom_hex(bins = 60) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "white", linewidth = 0.5) +
  geom_vline(xintercept = mfi_median, linetype = "dashed", colour = "white", linewidth = 0.5) +
  scale_fill_viridis_c(option = "inferno", name = "Cell count") +
  scale_y_continuous(limits = c(-3, 5)) +
  labs(title = "G3 - LC3B MFI vs Skewness (cell density)",
       x = "LC3B MFI arcsinh", y = "LC3B Skewness") +
  theme_bw(base_size = 11)

p3_shh <- ggplot(df_plot[df_plot$group == "SHH", ],
                 aes(x = LC3B_MFI_arcsinh, y = LC3B_skew)) +
  geom_hex(bins = 60) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "white", linewidth = 0.5) +
  geom_vline(xintercept = mfi_median, linetype = "dashed", colour = "white", linewidth = 0.5) +
  scale_fill_viridis_c(option = "viridis", name = "Cell count") +
  scale_y_continuous(limits = c(-3, 5)) +
  labs(title = "SHH - LC3B MFI vs Skewness (cell density)",
       x = "LC3B MFI arcsinh", y = "LC3B Skewness") +
  theme_bw(base_size = 11)

ggsave(file.path(OUT_PLOTS_SKEW, "LC3B_03_density_hexbin_G3vsSHH.pdf"),
       p3_g3 + p3_shh, width = 14, height = 6)
message("Salvato: LC3B_03_density_hexbin_G3vsSHH.pdf")

# --------------------------------------------------------------------------
# PLOT 4 — Mediane per paziente: MFI vs Skewness (unità biologica)
# --------------------------------------------------------------------------

pat_summary <- df_all %>%
  group_by(group, patient_id) %>%
  summarise(
    LC3B_MFI_median  = median(LC3B_MFI_arcsinh, na.rm = TRUE),
    LC3B_skew_median = median(LC3B_skew,         na.rm = TRUE),
    n_cells          = n(),
    .groups = "drop"
  )

message("\nMediane LC3B per paziente (MFI arcsinh + Skewness):")
print(as.data.frame(pat_summary))

p4 <- ggplot(pat_summary,
             aes(x = LC3B_MFI_median, y = LC3B_skew_median,
                 colour = group, label = patient_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = median(pat_summary$LC3B_MFI_median),
             linetype = "dashed", colour = "grey50") +
  geom_point(size = 5, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20,
                            show.legend = FALSE) +
  scale_colour_manual(values = PALETTE_GROUP) +

  # Quadrant labels
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
           label = "Sparse punctate\n(low MFI, high skewness)", size = 3,
           colour = "grey40", fontface = "italic") +
  annotate("text", x = Inf,  y = Inf, hjust = 1.05, vjust = 1.3,
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

# ggrepel opzionale — usa geom_text se non disponibile
tryCatch({
  ggsave(file.path(OUT_PLOTS_SKEW, "LC3B_04_patient_medians_MFI_vs_Skew.pdf"),
         p4, width = 9, height = 7)
  message("Salvato: LC3B_04_patient_medians_MFI_vs_Skew.pdf")
}, error = function(e) {
  # fallback senza ggrepel
  p4b <- p4
  p4b$layers <- p4b$layers[!sapply(p4b$layers, function(l)
    inherits(l$geom, "GeomTextRepel"))]
  p4b <- p4b + geom_text(size = 3, vjust = -0.8)
  ggsave(file.path(OUT_PLOTS_SKEW, "LC3B_04_patient_medians_MFI_vs_Skew.pdf"),
         p4b, width = 9, height = 7)
  message("Salvato: LC3B_04_patient_medians_MFI_vs_Skew.pdf (senza ggrepel)")
})

# --------------------------------------------------------------------------
# STATISTICHE: dove si concentrano le cellule G3 vs SHH nello spazio MFI×skew?
# --------------------------------------------------------------------------

message("\n=== STATISTICHE DIAGNOSTICHE ===")

for (grp in c("G3", "SHH")) {
  sub <- df_all[df_all$group == grp, ]
  n   <- nrow(sub)
  message(sprintf("\n%s (n=%s cellule):", grp, format(n, big.mark = ",")))
  message(sprintf("  LC3B MFI arcsinh:  median=%.3f  IQR=[%.3f, %.3f]",
                  median(sub$LC3B_MFI_arcsinh),
                  quantile(sub$LC3B_MFI_arcsinh, 0.25),
                  quantile(sub$LC3B_MFI_arcsinh, 0.75)))
  message(sprintf("  LC3B Skewness:     median=%.3f  IQR=[%.3f, %.3f]",
                  median(sub$LC3B_skew),
                  quantile(sub$LC3B_skew, 0.25),
                  quantile(sub$LC3B_skew, 0.75)))

  # % cellule nei 4 quadranti (usando mediana complessiva come soglia MFI)
  hi_mfi   <- sub$LC3B_MFI_arcsinh >= mfi_median
  pos_skew <- sub$LC3B_skew > 0
  message(sprintf("  Quadranti (soglia MFI=mediana complessiva=%.2f):", mfi_median))
  message(sprintf("    Alta MFI + skew>0  (denso puntato):  %.1f%%",
                  100 * mean(hi_mfi &  pos_skew)))
  message(sprintf("    Alta MFI + skew≤0  (uniforme/neg.):  %.1f%%",
                  100 * mean(hi_mfi & !pos_skew)))
  message(sprintf("    Bassa MFI + skew>0 (sparso puntato): %.1f%%",
                  100 * mean(!hi_mfi &  pos_skew)))
  message(sprintf("    Bassa MFI + skew≤0 (diffuso):        %.1f%%",
                  100 * mean(!hi_mfi & !pos_skew)))
}

# Test Wilcoxon G3 vs SHH su MFI e skewness (unità = paziente)
pat_g3  <- pat_summary$LC3B_MFI_median[pat_summary$group == "G3"]
pat_shh <- pat_summary$LC3B_MFI_median[pat_summary$group == "SHH"]
wt_mfi  <- wilcox.test(pat_g3, pat_shh, exact = TRUE)
message(sprintf("\nWilcoxon LC3B MFI G3 vs SHH: W=%.0f, p=%.4f",
                wt_mfi$statistic, wt_mfi$p.value))

pat_g3s  <- pat_summary$LC3B_skew_median[pat_summary$group == "G3"]
pat_shhs <- pat_summary$LC3B_skew_median[pat_summary$group == "SHH"]
wt_skew  <- wilcox.test(pat_g3s, pat_shhs, exact = TRUE)
message(sprintf("Wilcoxon LC3B Skew G3 vs SHH: W=%.0f, p=%.4f",
                wt_skew$statistic, wt_skew$p.value))

message("\n=== ANALISI COMPLETATA ===")
message(sprintf("Plot salvati in: %s", OUT_PLOTS_SKEW))
message(sprintf("Dati salvati in: %s", OUT_DATA))
