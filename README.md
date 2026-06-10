# case_study_biomarker_ST000369

> Short description of the project.

## Structure

```
case_study_biomarker_ST000369/
├── case_study_biomarker_ST000369.Rproj  # RStudio project file
├── rawdata/        # Immutable raw input data
├── data/           # Processed/intermediate data
├── code/           # Scripts and modules
│   └── utils/      # Shared helper functions
├── notebooks/      # Exploratory notebooks (.ipynb, .Rmd, .qmd)
│   ├── output/     # Rendered reports
│   └── old/        # Archived notebooks
├── output/         # Final results
│   ├── figures/
│   ├── tables/
│   └── reports/
├── docs/           # Documentation and notes
├── logs/           # Pipeline logs (git-ignored)
├── config/         # Parameters and config files
└── envs/           # Conda/pip env files + renv.lock
```

## Setup

**Python (conda):**
\`\`\`bash
conda env create -f envs/environment.yml
conda activate ${PROJECT}
\`\`\`

**R (renv):**
\`\`\`r
# Inside R / RStudio — restores all packages from renv.lock
renv::restore()
\`\`\`

## Usage

_Describe how to run the analysis here._
