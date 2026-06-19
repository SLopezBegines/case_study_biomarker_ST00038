# =============================================================================
# Central parameters for the ST000368 (ADC1) / ST000369 (ADC2) biomarker
# case study. Sourced by all analysis scripts. All comments in English.
# Paths are absolute to the project root on this machine.
# =============================================================================

PROJ <- "/home/santi/github_repos/propietary/case_study_biomarker_ST000368"
# Output locations for this notebook
NB_OUT <- file.path(PROJ, "notebooks", "output", "ST000368")

# NB_OUT should be personalized into the notebook
RAW <- file.path(PROJ, "rawdata")
DATA <- file.path(NB_OUT, "data")
FIG <- file.path(NB_OUT, "figures")
TAB <- file.path(NB_OUT, "tables")
REP <- file.path(NB_OUT, "reports")
RDAT <- file.path(NB_OUT, "RData")
MET <- file.path(NB_OUT, "metabo")

for (d in c(NB_OUT, DATA, FIG, TAB, REP, RDAT, MET)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ==============================================================================
# Global Variables
# ==============================================================================

# --- Output file extensions ---------------------------------------------------
tiff_extension <- ".tiff"
pdf_extension <- ".pdf"

# --- Species / organism -------------------------------------------------------
species <- 9606 # NCBI Taxonomy ID for Homo sapiens #9606 for Human, 10090 for mouse, 7955 for zebrafish
organism <- "org.Hs.eg.db" # Bioconductor annotation package #"org.Dr.eg.db". "org.Hs.eg.db"
kegg_organism <- "hsa" # KEGG organism code for Homo sapiens #dre for Danio rerio, hsa for Homo sapiens, mmu for Mus musculus
keyType <- "SYMBOL"
KEGGkeyType <- "kegg"

# --- Memory / hardware settings (i7-7560U: 2 physical / 4 logical cores, 16 GB RAM) --
# Physical cores = 2; sequential plan avoids fork-overhead and protects RAM
# Swap available: 16 GB — can absorb moderate overflow but slows analysis
PARALLEL_WORKERS <- 1 # Sequential: safer than parallel on 2-core laptop
FUTURE_GLOBALS_MAX_MB <- 6000 # 6 GB global size limit for {future} (conservative)
options(future.globals.maxSize = FUTURE_GLOBALS_MAX_MB * 1024^2)

# --- Checkpoint naming --------------------------------------------------------
# Checkpoints are stored as: output_path/RData/checkpoint_<NAME>.rds
# This allows recovery from any step if the session crashes
CHECKPOINT_PREFIX <- "checkpoint_"


# Analysis configuration ------------------------------------------------------
PARAMS <- list(
  # Processed REST data tables (named metabolites + Factors row)
  adc1_file = file.path(RAW, "MSdata_ST000368_1.txt"), # ADC1 = training cohort
  adc2_file = file.path(RAW, "MSdata_ST000369_1.txt"), # ADC2 = validation cohort

  # Class definition: Adenocarcinoma vs Healthy (binary)
  positive_class = "Adenocarcinoma",
  negative_class = "Healthy",

  # Histologies dropped from the binary problem (kept out of all matrices)
  drop_status = c("Adenosquamous"),

  # Known label typo in ADC2 deposited data -> canonical spelling
  status_fix = c("Adenocarcnoma" = "Adenocarcinoma"),

  # QC pools carry no cancer status (NA) and are excluded automatically
  biofluids = c("Plasma", "Serum"),
  seed = 1234L
)


# dpi for Metaboanalyst
if (!exists("default.dpi")) default.dpi <- 300 # MetaboAnalystR default
