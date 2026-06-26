# =============================================================================
# 01_build_matrices.R
# Build the analysis matrices for the dual-track biomarker case study.
#
# Design (see docs/case_study_ML_biomarker_design.md):
#   - Binary problem: Adenocarcinoma vs Healthy.
#   - Cross-cohort: train on ADC1 (ST000368), validate on ADC2 (ST000369).
#   - Per biofluid (Plasma, Serum) modelled separately -> one sample per
#     patient within a biofluid -> no patient-pairing leakage.
#   - Features restricted to metabolites shared by both cohorts (138).
#
# Outputs (all under data/ and output/tables/, never rawdata/):
#   data/analysis_matrices.rds       : list of per-(cohort,biofluid) frames
#   data/<cohort>_<biofluid>.csv      : tidy samples x (class + metabolites)
#   config/sample_sheet.csv           : all biological samples + factors
#   output/tables/shared_metabolites.csv
#   output/tables/cohort_composition.csv
# =============================================================================


set.seed(PARAMS$seed)
dir.create(DATA, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)

# ---- Load both cohorts ------------------------------------------------------
adc1 <- read_msdata(PARAMS$adc1_file)
adc2 <- read_msdata(PARAMS$adc2_file)

# ---- Normalise cancer-status labels (fix known typo) ------------------------
fix_status <- function(x) {
  x <- trimws(x)
  for (k in names(PARAMS$status_fix)) x[x == k] <- PARAMS$status_fix[[k]]
  x
}
adc1$factors$Cancer_status <- fix_status(adc1$factors$Cancer_status)
adc2$factors$Cancer_status <- fix_status(adc2$factors$Cancer_status)

# ---- Shared metabolites (cross-cohort feature space) ------------------------
shared <- intersect(rownames(adc1$mat), rownames(adc2$mat))
cat(sprintf(
  "ADC1 metabolites: %d | ADC2 metabolites: %d | shared: %d\n",
  nrow(adc1$mat), nrow(adc2$mat), length(shared)
))

shared_tbl <- merge(
  adc1$refmet[adc1$refmet$Metabolite_name %in% shared, ],
  adc2$refmet[adc2$refmet$Metabolite_name %in% shared, ],
  by = "Metabolite_name", all = TRUE, suffixes = c("_ADC1", "_ADC2")
)
write.csv(shared_tbl, file.path(TAB, "shared_metabolites.csv"), row.names = FALSE)

# ---- Build one tidy frame per (cohort, biofluid) ----------------------------
classes <- c(PARAMS$positive_class, PARAMS$negative_class)

build_frame <- function(obj, cohort, biofluid) {
  f <- obj$factors
  keep <- which(f$Organ == biofluid & f$Cancer_status %in% classes)
  if (length(keep) == 0) {
    return(NULL)
  }
  ids <- rownames(f)[keep]
  X <- t(obj$mat[shared, ids, drop = FALSE]) # samples x metabolites
  meta <- data.frame(
    sample_id = ids,
    cohort = cohort,
    biofluid = biofluid,
    class = factor(f$Cancer_status[keep], levels = classes),
    smoker = f$Smoker[keep],
    gender = f$Gender[keep],
    stringsAsFactors = FALSE
  )
  cbind(meta, as.data.frame(X, check.names = FALSE))
}

matrices <- list()
for (cohort in c("ADC1", "ADC2")) {
  obj <- if (cohort == "ADC1") adc1 else adc2
  for (bf in PARAMS$biofluids) {
    key <- paste(cohort, bf, sep = "_")
    fr <- build_frame(obj, cohort, bf)
    if (!is.null(fr)) {
      matrices[[key]] <- fr
      write.csv(fr, file.path(DATA, paste0(key, ".csv")), row.names = FALSE)
    }
  }
}


# ---- Sample sheet + composition summary ------------------------------------
sample_sheet <- do.call(rbind, lapply(matrices, function(fr) {
  fr[, c("sample_id", "cohort", "biofluid", "class", "smoker", "gender")]
}))
rownames(sample_sheet) <- NULL
write.csv(sample_sheet, file.path(DATA, "sample_sheet.csv"), row.names = FALSE)

comp <- as.data.frame(table(
  sample_sheet$cohort, sample_sheet$biofluid,
  sample_sheet$class
))
names(comp) <- c("cohort", "biofluid", "class", "n")
comp <- comp[comp$n > 0, ]
write.csv(comp, file.path(TAB, "cohort_composition.csv"), row.names = FALSE)


analysis_matrices <- list(
  matrices = matrices, shared = shared,
  refmet = shared_tbl, params = PARAMS,
  cohort_composition = comp, samples = sample_sheet
)
saveRDS(
  analysis_matrices, file.path(RDAT, "analysis_matrices.rds")
)


cat("\n=== Composition (biological, Adenocarcinoma vs Healthy) ===\n")
print(comp)
cat("\nMatrices built:", paste(names(matrices), collapse = ", "), "\n")
cat("Feature count (shared metabolites):", length(shared), "\n")


rm(matrices, shared_tbl, shared, adc1, adc2, sample_sheet, comp, fr, obj, bf, classes, build_frame, cohort)
