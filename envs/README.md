# envs/

Environment and dependency files.

- `environment.yml`   — conda environment (recommended for reproducibility).
- `requirements.txt`  — pip packages.

- `renv.lock`         — R package snapshot (committed to git for reproducibility).
  Restore with: `Rscript -e 'renv::restore()'`
