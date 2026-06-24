# =============================================================================
# Metabolomics normalization helpers and limma wrapper.
# Replaces the MetaboAnalystR QC/DE pipeline for Track-A analysis.
# Pipeline: median row-norm → log2 → [Pareto for PCA] → limma.
# All comments in English.
# =============================================================================

# -- 1. Feature filtering -----------------------------------------------------

#' Remove features with coefficient of variation below a threshold.
#' Constant or near-constant features add no information and destabilise limma.
#' @param mat    Numeric matrix, samples in rows, features in columns
#' @param min_cv Minimum CV to retain (default 0.10 = 10 %)
#' @returns Filtered matrix; prints a message listing dropped features
filter_low_var <- function(mat, min_cv = 0.10) {
  cv <- apply(mat, 2, function(x) {
    m <- mean(x, na.rm = TRUE)
    if (is.na(m) || m == 0) 0 else sd(x, na.rm = TRUE) / abs(m)
  })
  keep <- cv >= min_cv
  if (any(!keep))
    message(sprintf("[FILTER] Dropped %d low-variance feature(s) (CV < %.2f): %s",
                    sum(!keep), min_cv,
                    paste(names(cv)[!keep], collapse = ", ")))
  mat[, keep, drop = FALSE]
}

# -- 2. Normalization steps ---------------------------------------------------

#' Median row normalization — corrects sample loading/dilution differences.
#' Each sample is divided by its median metabolite intensity.
#' @param mat Numeric matrix, samples in rows, features in columns
#' @returns Normalized matrix (same dimensions)
median_row_norm <- function(mat) {
  meds <- apply(mat, 1, median, na.rm = TRUE)
  meds[is.na(meds) | meds == 0] <- 1
  sweep(mat, 1, meds, "/")
}

#' Log2 transformation with per-feature pseudocount (min positive / 2).
#' Log2 is chosen so limma's logFC is directly in log2 units (consistent
#' with RNA-seq conventions and the log2(FC_THRESH) threshold).
#' @param mat Numeric matrix, samples in rows, features in columns
#' @returns Log2-transformed matrix
log2_transform <- function(mat) {
  pseudo <- apply(mat, 2, function(x) {
    pos <- x[x > 0 & !is.na(x)]
    if (length(pos) > 0) min(pos) / 2 else 1e-6
  })
  log2(sweep(mat, 2, pseudo, "+"))
}

#' Pareto feature scaling: mean-centre then divide by sqrt(SD).
#' Applied ONLY for PCA/multivariate visualisation.
#' NOT used before limma: Pareto changes the scale of logFC, making
#' fold-change thresholds uninterpretable.
#' @param mat Numeric matrix, samples in rows, features in columns
#' @returns Pareto-scaled matrix
pareto_scale <- function(mat) {
  means <- colMeans(mat, na.rm = TRUE)
  sds   <- apply(mat, 2, sd, na.rm = TRUE)
  sds[is.na(sds) | sds <= 0] <- 1        # protect against zero-variance columns
  sweep(sweep(mat, 2, means, "-"), 2, sqrt(sds), "/")
}

# -- 3. Full pipeline wrapper -------------------------------------------------

#' Prepare one metabolomics matrix: filter → normalize for DE and PCA.
#'
#' Returns a list with four named elements:
#'   $raw     — feature matrix after low-var filter only (no transformation)
#'   $mat_de  — log2-transformed, median row-normalised  (for limma)
#'   $mat_pca — log2 + Pareto-scaled                     (for PCA/EDA)
#'   $meta    — metadata data.frame (META_COLS rows)
#'
#' Row names of all matrices = sample_id.
#'
#' @param df        Data frame with META_COLS + feature columns
#' @param feat_cols Character vector of feature (metabolite) column names
#' @param meta_cols Character vector of metadata column names (default META_COLS)
#' @param min_cv    CV threshold for feature filtering (default 0.10)
#' @returns Named list as described above
prepare_matrix <- function(df, feat_cols,
                           meta_cols = META_COLS,
                           min_cv    = 0.10) {
  raw_mat <- as.matrix(df[, feat_cols, drop = FALSE])
  storage.mode(raw_mat) <- "double"
  rownames(raw_mat) <- df$sample_id

  filtered <- filter_low_var(raw_mat, min_cv)
  log2_mat <- log2_transform(median_row_norm(filtered))

  list(
    raw     = filtered,
    mat_de  = log2_mat,
    mat_pca = pareto_scale(log2_mat),
    meta    = df[, meta_cols, drop = FALSE]
  )
}

# -- 4. Covariate diagnostics -------------------------------------------------

#' Cross-tabulate class vs covariates and check for perfect aliasing.
#' Aliased covariates would make the limma model non-estimable.
#' @param meta Data frame with columns: class, smoker, gender
#' @returns Invisibly: list of cross-tables
check_covariates <- function(meta) {
  cat("\n=== Class × Smoker ===\n")
  print(table(class = meta$class, smoker = meta$smoker, useNA = "ifany"))
  cat("\n=== Class × Gender ===\n")
  print(table(class = meta$class, gender = meta$gender, useNA = "ifany"))

  n_na <- sum(is.na(meta$smoker)) + sum(is.na(meta$gender))
  if (n_na > 0)
    warning("[COVARIATES] ", n_na,
            " NA(s) in smoker/gender — those samples are dropped by model.matrix().")

  fit_chk <- tryCatch(
    lm(as.numeric(factor(meta$class)) ~ meta$smoker + meta$gender),
    error = function(e) { message("[COVARIATES] lm check failed: ", e$message); NULL }
  )
  if (!is.null(fit_chk)) {
    al <- alias(fit_chk)$Complete
    if (!is.null(al) && nrow(al) > 0)
      warning("[COVARIATES] Perfect aliasing: ", paste(rownames(al), collapse = ", "),
              ". Drop the aliased covariate before fitting.")
    else
      message("[COVARIATES] No aliasing detected — model is estimable.")
  }

  invisible(list(
    class_smoker = table(meta$class, meta$smoker, useNA = "ifany"),
    class_gender = table(meta$class, meta$gender, useNA = "ifany")
  ))
}

# -- 5. limma differential analysis -------------------------------------------

#' Run limma moderated t-test (with empirical Bayes shrinkage).
#'
#' Model: ~ class + smoker + gender
#' Contrast: pos_class vs neg_class (intercept = neg_class).
#' trend = TRUE recommended for metabolomics: assumes variance-mean
#' relationship similar to what is seen in RNA-seq.
#'
#' logFC is in log2 units (input = log2-transformed data).
#' Significance threshold: adj_p_value < PVAL_THRESH AND
#'                         |logFC| > log2(FC_THRESH).
#'
#' @param mat_de    Log2 matrix, samples in rows, features in columns
#' @param meta      Metadata data.frame; must have columns sample_id, class,
#'                  smoker, gender; rows aligned with mat_de rows
#' @param label     Label string for the result column (e.g. "ADC1_Plasma")
#' @param pos_class Positive class (numerator of FC; default from PARAMS)
#' @param neg_class Reference class (denominator; default from PARAMS)
#' @returns data.frame with per-feature DE statistics
run_limma_de <- function(mat_de, meta,
                         label     = "",
                         pos_class = PARAMS$positive_class,
                         neg_class = PARAMS$negative_class) {

  # Drop samples with NA covariates
  complete_idx <- !is.na(meta$smoker) & !is.na(meta$gender)
  if (any(!complete_idx)) {
    message("[LIMMA] Dropping ", sum(!complete_idx),
            " sample(s) with NA smoker/gender.")
    meta   <- meta[complete_idx, , drop = FALSE]
    mat_de <- mat_de[meta$sample_id, , drop = FALSE]
  }

  # Align rows (sample order must match between matrix and metadata)
  meta   <- meta[match(rownames(mat_de), meta$sample_id), , drop = FALSE]

  meta$class  <- factor(meta$class,  levels = c(neg_class, pos_class))
  meta$smoker <- droplevels(factor(meta$smoker))
  meta$gender <- droplevels(factor(meta$gender))

  design <- model.matrix(~ class + smoker + gender, data = meta)
  design <- design[, colSums(abs(design)) > 0, drop = FALSE]  # drop empty cols

  fit <- limma::lmFit(t(mat_de), design)
  fit <- limma::eBayes(fit, trend = TRUE)

  coef_name <- paste0("class", pos_class)
  if (!coef_name %in% colnames(fit$coefficients))
    stop("[LIMMA] Coefficient not found: '", coef_name,
         "'. Available: ", paste(colnames(fit$coefficients), collapse = ", "))

  tt <- limma::topTable(fit, coef = coef_name, number = Inf,
                        adjust.method = "BH", sort.by = "P")

  data.frame(
    metabolite  = rownames(tt),
    logFC       = tt$logFC,         # log2 fold change (ADC / Healthy)
    AveExpr     = tt$AveExpr,
    t_stat      = tt$t,
    p_value     = tt$P.Value,
    adj_p_value = tt$adj.P.Val,
    B           = tt$B,             # log-odds of DE
    significant = tt$adj.P.Val < PVAL_THRESH & abs(tt$logFC) > log2(FC_THRESH),
    direction   = ifelse(tt$logFC > 0, "Up", "Down"),
    label       = label,
    row.names   = NULL,
    stringsAsFactors = FALSE
  )
}
