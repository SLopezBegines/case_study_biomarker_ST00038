# ==============================================================================
# Package Loading
# Uses pak for all installations: parallel, handles Bioc + GitHub, fast.
# First-time use: run `sudo apt-get install -y libuv1-dev` before this script.
# ==============================================================================

# Bootstrap pak itself --------------------------------------------------------
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = sprintf(
    "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
    .Platform$pkgType, R.Version()$os, R.Version()$arch
  ))
}
library(pak)

# Print system requirements for any missing packages (informational) ----------
check_sysreqs <- function(pkgs) {
  tryCatch(
    {
      reqs <- pak::pkg_sysreqs(pkgs)
      if (length(reqs$packages) > 0) {
        message("[SYSREQS] Missing system packages detected:")
        message("  Run: sudo apt-get install -y ", paste(reqs$packages, collapse = " "))
      }
    },
    error = function(e) invisible(NULL)
  )
}

# Helper: install + load, skipping already-installed -------------------------
# pkgs can use pak prefixes (bioc::Pkg, user/repo) for install;
# library() needs bare names, so strip everything up to and including "::".
load_pkgs <- function(pkgs) {
  bare <- sub("^[^:]+::", "", pkgs) # strip bioc:: / cran:: etc.
  bare <- sub("^.+/", "", bare) # strip user/ from GitHub refs
  missing_idx <- !vapply(bare, requireNamespace, logical(1), quietly = TRUE)
  if (any(missing_idx)) {
    message("[INSTALL] Installing: ", paste(pkgs[missing_idx], collapse = ", "))
    pak::pkg_install(pkgs[missing_idx], ask = FALSE, upgrade = FALSE)
  }
  invisible(lapply(bare, function(p) {
    suppressPackageStartupMessages(suppressWarnings(
      library(p, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
    ))
  }))
}

# =============================================================================
# CRAN packages  (includes Seurat + SeuratObject — they live on CRAN, not Bioc)
# =============================================================================
message("[PACKAGES] Loading CRAN packages...")

cran_packages <- c(
  # Project management
  "renv",
  # Core data manipulation
  "tidyverse", "data.table",
  # Visualization
  "patchwork", "ggrepel", "scales",
  # I/O
  "readxl", "writexl",
  # Utilities
  "remotes", "future", "future.apply",
  # Clustering diagnostics
  "clustree",
  # Counting time
  "tictoc",
  # Install Github
  "devtools",
  # Bioconductor installer
  "BiocManager",
  # Markdown tables
  "kableExtra",
  # Missing Values
  "DataExplorer", "VIM",
  # Differential analysis (non-MetaboAnalystR pipeline)
  "pheatmap",      # heatmaps for significant metabolites
  "factoextra",    # PCA biplots (ggplot2-based)
  "broom"          # tidy model summaries
)

metaboanalyst_packages <- c(
  # MetaboAnalystR dependencies
  "impute", "pcaMethods", "globaltest",
  "GlobalAncova", "Rgraphviz", "preprocessCore",
  "genefilter", "sva", "limma", "KEGGgraph", "siggenes", "BiocParallel", "MSnbase",
  "multtest", "RBGL", "edgeR", "fgsea", "httr", "qs2", "RSclient",
  # MetaboAnalystR itself (GitHub)
  "xia-lab/MetaboAnalystR"
)

load_pkgs(cran_packages)
load_pkgs(metaboanalyst_packages)
# =============================================================================
# Bioconductor packages  (bioc:: prefix — these are genuine Bioc packages)
# =============================================================================
message("[PACKAGES] Loading Bioconductor packages...")

bioc_packages <- paste0("bioc::", c(
  "BiocParallel",
  "SummarizedExperiment",
  "S4Vectors",
  "limma",           # moderated t-test for continuous -omics (metabolomics DE)
  "edgeR",           # also used by limma helpers; kept here for explicit loading
  "EnhancedVolcano",
  "clusterProfiler", # GO enrichment
  "STRINGdb",        # STRING PPI database
  "enrichplot",      # GO enrichment visualization
  "org.Hs.eg.db"
))
load_pkgs(bioc_packages)
# Session info
message("\n[PACKAGES] All packages loaded.")
cat("R version:", R.version$version.string, "\n")
cat("Tidyverse version:", as.character(packageVersion("tidyverse")), "\n")
cat("Bioconductor version:", as.character(packageVersion("BiocManager")), "\n")

# Clean environment
rm(check_sysreqs, load_pkgs, cran_packages, metaboanalyst_packages, bioc_packages)

# renv::snapshot(type = "all") # Uncomment to save package versions to renv.lock
# install.packages("RSclient", repos = "http://www.rforge.net/")
