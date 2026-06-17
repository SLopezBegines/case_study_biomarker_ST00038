# =============================================================================
# 02_eda_qc.R
# Stage 1 EDA / QC for the dual-cohort lung-adenocarcinoma metabolomics study.
#
# Goals:
#   - Quantify missingness / zeros on the 137 shared metabolites.
#   - Characterise intensity distributions and cross-cohort comparability.
#   - Batch diagnostic: which factor (cohort/batch, organ, class) drives the
#     main axes of variation? Quantified via per-PC ANOVA R^2.
#
# Outputs (output/figures, output/tables, output/reports):
#   eda_missingness.csv, eda_pc_variance_explained.csv
#   eda_intensity_distribution.png, eda_pca_global.png, eda_pca_by_biofluid.png
#   eda_summary.md
# Imputation here is for EDA visualisation only (per-feature median);
# modelling imputation is performed in-fold downstream.
# All comments in English.
# =============================================================================

suppressWarnings(source("config/params.R"))
set.seed(PARAMS$seed)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
dir.create(REP, showWarnings = FALSE, recursive = TRUE)

obj      <- readRDS(file.path(DATA, "analysis_matrices.rds"))
matrices <- obj$matrices
shared   <- obj$shared
meta_cols <- c("sample_id", "cohort", "biofluid", "class", "smoker", "gender")

# Combine all biological samples into one tidy frame -------------------------
all_df <- do.call(rbind, matrices)
rownames(all_df) <- all_df$sample_id
Xall <- as.matrix(all_df[, shared, drop = FALSE])      # samples x metabolites
meta <- all_df[, meta_cols]

# ---- Missingness / zeros per (cohort, biofluid) -----------------------------
miss_tbl <- do.call(rbind, lapply(names(matrices), function(k) {
  fr <- matrices[[k]]; X <- as.matrix(fr[, shared, drop = FALSE])
  data.frame(matrix = k,
             n_samples = nrow(X),
             n_features = ncol(X),
             pct_NA   = round(100 * mean(is.na(X)), 3),
             pct_zero = round(100 * mean(X == 0, na.rm = TRUE), 3),
             feat_max_NA_pct = round(100 * max(colMeans(is.na(X))), 1),
             stringsAsFactors = FALSE)
}))
write.csv(miss_tbl, file.path(TAB, "eda_missingness.csv"), row.names = FALSE)

# ---- EDA-only imputation (per-feature median) + log10 transform -------------
impute_median <- function(X) {
  for (j in seq_len(ncol(X))) {
    v <- X[, j]; v[is.na(v)] <- median(v, na.rm = TRUE); X[, j] <- v
  }
  X
}
Xlog <- log10(impute_median(Xall) + 1)

# ---- Intensity distribution by cohort --------------------------------------
png(file.path(FIG, "eda_intensity_distribution.png"), 1100, 600, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
cols_co <- c(ADC1 = "#1b9e77", ADC2 = "#d95f02")
plot(NA, xlim = range(Xlog), ylim = c(0, 1.1),
     xlab = "log10(intensity + 1)", ylab = "density",
     main = "Intensity density by cohort")
for (co in c("ADC1", "ADC2")) {
  d <- density(as.numeric(Xlog[meta$cohort == co, ]))
  lines(d$x, d$y / max(d$y), col = cols_co[co], lwd = 2)
}
legend("topright", names(cols_co), col = cols_co, lwd = 2, bty = "n")
boxplot(rowMeans(Xlog) ~ meta$cohort, col = cols_co,
        xlab = "cohort", ylab = "mean log10 intensity / sample",
        main = "Per-sample mean intensity")
dev.off()

# ---- PCA (global) -----------------------------------------------------------
pca <- prcomp(Xlog, center = TRUE, scale. = TRUE)
ve  <- (pca$sdev^2) / sum(pca$sdev^2)
scores <- as.data.frame(pca$x[, 1:5])

# Per-PC ANOVA R^2 against each factor (batch diagnostic) --------------------
r2 <- function(pc, g) summary(lm(pc ~ g))$r.squared
pc_r2 <- data.frame(
  PC      = paste0("PC", 1:5),
  var_pct = round(100 * ve[1:5], 1),
  R2_cohort = round(sapply(1:5, function(i) r2(scores[, i], meta$cohort)), 3),
  R2_organ  = round(sapply(1:5, function(i) r2(scores[, i], meta$biofluid)), 3),
  R2_class  = round(sapply(1:5, function(i) r2(scores[, i], meta$class)), 3))
write.csv(pc_r2, file.path(TAB, "eda_pc_variance_explained.csv"), row.names = FALSE)

# Global PCA figure: colour by cohort / organ / class ------------------------
png(file.path(FIG, "eda_pca_global.png"), 1500, 520, res = 130)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
plot_pca <- function(grp, palette, title) {
  g <- factor(grp)
  plot(scores$PC1, scores$PC2, col = palette[g], pch = 19, cex = 0.9,
       xlab = sprintf("PC1 (%.1f%%)", 100 * ve[1]),
       ylab = sprintf("PC2 (%.1f%%)", 100 * ve[2]), main = title)
  legend("topright", levels(g), col = palette[levels(g)], pch = 19, bty = "n")
}
plot_pca(meta$cohort,   c(ADC1 = "#1b9e77", ADC2 = "#d95f02"), "by cohort (batch)")
plot_pca(meta$biofluid, c(Plasma = "#7570b3", Serum = "#e7298a"), "by biofluid")
plot_pca(meta$class,    setNames(c("#377eb8", "#e41a1c"),
         c(PARAMS$negative_class, PARAMS$positive_class)), "by class")
dev.off()

# ---- PCA per biofluid (train+val), colour by class, shape by cohort ---------
png(file.path(FIG, "eda_pca_by_biofluid.png"), 1100, 560, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (bf in PARAMS$biofluids) {
  idx <- meta$biofluid == bf
  p <- prcomp(Xlog[idx, ], center = TRUE, scale. = TRUE)
  v <- (p$sdev^2) / sum(p$sdev^2)
  cl <- factor(meta$class[idx], levels = c(PARAMS$negative_class, PARAMS$positive_class))
  co <- factor(meta$cohort[idx])
  plot(p$x[, 1], p$x[, 2], col = c("#377eb8", "#e41a1c")[cl],
       pch = c(1, 17)[co], cex = 1,
       xlab = sprintf("PC1 (%.1f%%)", 100 * v[1]),
       ylab = sprintf("PC2 (%.1f%%)", 100 * v[2]),
       main = paste0(bf, " (ADC1+ADC2)"))
  legend("topright", c(levels(cl), levels(co)),
         col = c("#377eb8", "#e41a1c", "black", "black"),
         pch = c(15, 15, 1, 17), bty = "n", cex = 0.8)
}
dev.off()

# ---- EDA summary report (Markdown, English) --------------------------------
lines <- c(
  "# Stage 1 EDA — ST000368 (ADC1) / ST000369 (ADC2)",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "## Composition (biological, Adenocarcinoma vs Healthy, per biofluid)",
  "",
  paste(capture.output(print(miss_tbl[, c("matrix","n_samples","n_features")],
                              row.names = FALSE)), collapse = "\n"),
  "",
  "## Missingness / zeros",
  "",
  paste(capture.output(print(miss_tbl, row.names = FALSE)), collapse = "\n"),
  "",
  "## PCA — variance explained and per-PC ANOVA R^2 (batch diagnostic)",
  "",
  "R2_cohort = fraction of each PC explained by cohort/batch (2013 vs 2014);",
  "R2_organ = by biofluid; R2_class = by Adenocarcinoma-vs-Healthy.",
  "",
  paste(capture.output(print(pc_r2, row.names = FALSE)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  sprintf("- PC1 explains %.1f%% of variance; dominant factor: %s.",
          100 * ve[1],
          c("cohort/batch","biofluid","class")[which.max(
            c(pc_r2$R2_cohort[1], pc_r2$R2_organ[1], pc_r2$R2_class[1]))]),
  "- If cohort or biofluid dominate the top PCs, the class signal is weaker",
  "  than the batch/biofluid structure -> reinforces (a) per-biofluid",
  "  modelling and (b) fitting normalisation on ADC1 only and treating the",
  "  ADC1->ADC2 drop as partly batch-driven.",
  "",
  "Figures: eda_intensity_distribution.png, eda_pca_global.png,",
  "eda_pca_by_biofluid.png (output/figures/)."
)
writeLines(lines, file.path(REP, "eda_summary.md"))

cat("EDA done.\n"); print(miss_tbl); cat("\nPer-PC R^2:\n"); print(pc_r2)

