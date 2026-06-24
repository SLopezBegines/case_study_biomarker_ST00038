# Biomarker Discovery — Lung Adenocarcinoma Metabolomics
### ST000368 (ADC1 training) · ST000369 (ADC2 validation)

> Differential metabolomics of serum and plasma from lung adenocarcinoma vs healthy controls.
> GC-TOF-MS, two independent cohorts, 137 shared metabolites.

---

## Background

This project identifies serum/plasma metabolite biomarkers of lung adenocarcinoma using two GC-TOF-MS cohorts from the NIH Metabolomics Workbench (PR000293):

- **ADC1 (ST000368)** — training cohort: differential analysis, biomarker discovery
- **ADC2 (ST000369)** — independent validation cohort (held out until final evaluation)

The binary comparison is **Adenocarcinoma vs Healthy** (Adenosquamous excluded).  
Analyses run **separately per biofluid** (Plasma, Serum) to avoid within-patient plasma–serum pairing leakage.

---

## Dataset

| Cohort | Study | Samples | Biofluids | Metabolites |
|---|---|---|---|---|
| ADC1 | ST000368 | 86 (43 ADC / 43 Healthy) | Plasma, Serum | 137 shared |
| ADC2 | ST000369 | TBD | Plasma, Serum | 137 shared |

Raw data: NIH Metabolomics Workbench REST API → `rawdata/MSdata_ST000368_1.txt`, `rawdata/MSdata_ST000369_1.txt`  
Processed matrices: `notebooks/output/ST000368/data/ADC1_{Plasma,Serum}.csv`

---

## Analysis Pipeline

### Track A — Univariate Differential (this notebook)

```
Raw MSdata
    │
    ▼
03_build_matrices.R          ← parse, filter, align 137 shared features
    │
    ▼
Missing value imputation     ← kNN (VIM, k=10); features >50% missing removed first
    │
    ▼
Normalization per matrix
  1. Median row normalization (sample loading correction)
  2. Log10 transform          (variance stabilization)
  3. Pareto feature scaling   (÷ √SD; preserves structure, reduces dominance)
    │
    ▼
QC / EDA
  • Density & boxplot before/after normalization
  • CV distribution
  • DataExplorer correlation heatmap
  • PCA score + biplot (PC1×PC2, ellipses by class)
    │
    ▼
Differential analysis — limma
  Model: ~ class + smoker + gender
  Contrast: Adenocarcinoma − Healthy
  FDR: Benjamini-Hochberg (adj.P.Val < 0.05 AND |logFC| > threshold)
    │
    ▼
Results
  • Volcano plot (ggplot2 + ggrepel)
  • Heatmap of significant metabolites (pheatmap)
  • Export: TSV tables, RDS objects, PDF/TIFF figures
```

---

## Design Decisions

### Why limma, not DESeq2/edgeR?
DESeq2 and edgeR assume count data (negative binomial). GC-TOF-MS intensities are continuous
and approximately log-normal after log-transform. limma applied to log-transformed data is
the validated standard for continuous -omics differential analysis and has been widely used
in metabolomics (Smyth 2004; Dunn et al. 2011).

### Why not MetaboAnalystR?
MetaboAnalystR was the original engine but this version (installed from GitHub) stores
intermediate state in fixed-name `.qs` files in the working directory, making it
incompatible with a function-wrapped, multi-matrix pipeline. The equivalent analysis is
implemented here with standard Bioconductor (limma) and CRAN packages.

### Normalization rationale
| Method | Reason |
|---|---|
| Median row norm | Conservative sample loading correction; robust to outlier metabolites |
| Log10 transform | Standard for GC-MS; converts multiplicative noise to additive; near-normality |
| Pareto scaling | Reduces dominance of high-variance features without amplifying noise in low-abundance metabolites (preferred over auto-scaling in metabolomics) |

### Covariate inclusion
- **`smoker`** — essential: smoking is the primary lung-cancer risk factor and independently alters the metabolome. Without this covariate, the cancer signal and the smoking signal are entangled.
- **`gender`** — advisable: sex differences in lipid and amino acid metabolism are well documented. With 86 samples the model has sufficient degrees of freedom.

Model: `~ class + smoker + gender`  
Aliasing is checked before fitting; if perfect confounding is detected a covariate is dropped.

---

## Repository Structure

```
case_study_biomarker_ST000368/
├── README.md                          ← this file
├── config/
│   ├── params.R                       ← paths, thresholds, global constants
│   └── README.md
├── code/
│   ├── 00_packages.R                  ← package loading (pak-based)
│   ├── 01_aux_functions.R             ← save_plot, save_table, checkpoints, logging
│   ├── 02_aux_metabolomics.R          ← MSdata parser (read_msdata)
│   ├── 03_build_matrices.R            ← build ADC1/ADC2 per-biofluid matrices
│   └── README.md
├── notebooks/
│   ├── case_study_biomarker.qmd       ← main analysis (Track A, limma-based)
│   └── output/ST000368/
│       ├── data/                      ← processed matrices (CSV)
│       ├── figures/                   ← PDF + TIFF plots
│       ├── tables/                    ← TSV result tables
│       ├── RData/                     ← checkpoint .rds objects
│       └── reports/
├── rawdata/
│   ├── MSdata_ST000368_1.txt          ← ADC1 raw (git-ignored if large)
│   ├── MSdata_ST000369_1.txt          ← ADC2 raw (git-ignored if large)
│   └── README.md
├── docs/
│   ├── case_study_ML_biomarker_design.md
│   └── README.md
├── envs/
│   └── README.md
└── renv/                              ← renv lockfile + library (reproducibility)
```

---

## How to Run

```r
# 1. Restore R environment
renv::restore()

# 2. Render the notebook (from project root)
quarto::quarto_render("notebooks/case_study_biomarker.qmd")

# Or render from the terminal:
# quarto render notebooks/case_study_biomarker.qmd
```

The notebook is self-contained: it sources `config/params.R`, `code/00_packages.R`, and the
aux function files automatically. All output goes to `notebooks/output/ST000368/`.

---

## Key Parameters (`config/params.R`)

| Parameter | Value | Description |
|---|---|---|
| `FC_THRESH` | 1.2 | Fold-change guide (|logFC| > log2(1.2) ≈ 0.263) |
| `PVAL_THRESH` | 0.05 | FDR threshold (adj.P.Val) |
| `PARALLEL_WORKERS` | 1 | Sequential (2-core laptop) |
| `PARAMS$seed` | 1234 | Random seed for reproducibility |

---

## Dependencies

Key R packages (full list in `code/00_packages.R`):

| Package | Role |
|---|---|
| `limma` | Differential analysis (moderated t-test, eBayes) |
| `VIM` | kNN imputation |
| `DataExplorer` | EDA plots |
| `ggplot2` + `patchwork` + `ggrepel` | Visualization |
| `factoextra` | PCA biplots |
| `pheatmap` | Heatmaps |
| `kableExtra` | Report tables |
| `tidyverse` | Data manipulation |

---

## References

- Smyth GK (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. *Statistical Applications in Genetics and Molecular Biology*, 3(1).
- Dunn WB et al. (2011). Procedures for large-scale metabolic profiling of serum and plasma using gas chromatography and liquid chromatography coupled to mass spectrometry. *Nature Protocols*, 6(7):1060–1083.
- van den Berg RA et al. (2006). Centering, scaling, and transformations: improving the biological information content of metabolomics data. *BMC Genomics*, 7:142.
- Worley B & Powers R (2013). Multivariate Analysis in Metabolomics. *Current Metabolomics*, 1(1):92–107.
