---
editor_options:
  markdown:
    wrap: 72
---

# Case Study Design — ML Biomarker Discovery (Serum Metabolomics, Lung Adenocarcinoma)

**Santiago López Begines, PhD** — *Junio 2026*

**Context:** Design document for a portfolio **case study** (not a generic
repo) demonstrating Service Line B — applied machine learning for
biomarker discovery. Target audience: **pharma / biotech companies**.
Selected problem type: feature selection + interpretability on
case-control omics. The deliverable is a *narrative analysis with a
business recommendation*, backed by a reproducible repo — not a bare
pipeline.

---

## 0. Why a "case study" and not another repo

A repo (e.g. the existing `snRNAseq_mouse`) is a **tool**: reusable,
generic code that proves technical execution. A case study is a
**business narrative**: a concrete question, justified analytical
decisions, an interpreted result, and an actionable recommendation. The
product is *judgement and communication*, not the script. The repo is
the engine underneath; the case study is the report on top.

**Thesis of this case study (the differentiator):**
> *"I do not sell you an overfit model. I sell you a validated model,
> honest about its uncertainty — and I can audit bio-ML that isn't."*

The single most common failure in bio-ML portfolios is reporting
`AUC ≈ 0.97` on n≈40 samples and 20,000 features, where the number is an
artefact of data leakage or selection bias. For a pharma audience, the
selling point is demonstrating that we **avoid** those traps and can
detect them in others' work.

---

## 1. Decision: Dataset — ST000368 (ADC1) + ST000369 (ADC2)

### Project PR000293 — Fahrmann et al. 2015, *Cancer Epidemiol Biomarkers Prev*

**"Investigation of metabolomic blood biomarkers for detection of
adenocarcinoma lung cancer"** — two independent case-control cohorts:
**ST000368 = ADC1** (discovery/training), **ST000369 = ADC2**
(independent validation, "part II").

#### Verified composition (from the downloaded processed data, not the paper)

| Parameter | ADC1 (ST000368) | ADC2 (ST000369) |
|---|---|---|
| Run date (from sample IDs) | Feb 2014 | Jul–Aug 2013 |
| Total samples | 191 | 181 |
| Adenocarcinoma | 86 | 94 (88 + 6 mislabelled "Adenocarcnoma") |
| Healthy | 86 | 63 |
| Other histology | — | Adenosquamous 6 |
| QC pools | 19 | 18 |
| Organ | Plasma 86 / Serum 86 | Plasma 101 / Serum 80 |
| Identified metabolites | 153 | 182 |

- **DOI:** 10.21228/M85P57 (PR000293) · Platform: GC-TOF-MS (Leco Pegasus).
- **Shared metabolites across cohorts: 138** (153 ∩ 182). Cross-cohort
  modelling (ADC1→ADC2) is restricted to these 138 features.
- **Class definition (decided):** Adenocarcinoma vs Healthy. Drop
  Adenosquamous (6) and QC pools.

### Why this dataset (the decisive feature)

- **It ships its own external validation cohort.** ADC1 and ADC2 are
  independent case-control sets run six months apart. This lets us
  *demonstrate* the thesis — validate on patients the model never saw,
  with no leakage — rather than merely claim it. Most portfolio datasets
  do not offer this.
- **Clinical-scale n** → the regime where multivariate prediction is
  actually licensed (see §1b), while still tractable on a laptop.
- **Clinically and commercially live problem** — early detection of lung
  adenocarcinoma from blood metabolomics, with recent ML literature.
- **Aligns with Line A1** (omics / mass-spec) while sitting squarely in
  Line B (ML).
- **Serum + plasma from the same patients** → enables per-biofluid
  modelling and a biofluid-robustness question, *and* sidesteps the
  patient-pairing leakage problem (see below).

### Verified data-quality findings (resolved, not assumptions)

- **Duplicate-download incident (resolved):** the bulk `.zip` for both
  accessions contained the *same* ADC2 processed table (identical md5);
  ADC1 was not in the zip. The real ADC1 was recovered via the
  Metabolomics Workbench REST API (`MSdata_ST000368_1.txt`,
  `ST000368_AN000602.*`). The mislabelled `ST000368/` zip folder is kept
  untouched and documented; do not analyse it.
- **Paper discrepancy (flagged, not hidden):** the article/portal
  describe ADC1 as 52 adenocarcinoma + 31 controls; the *deposited data*
  shows 86 / 86 + 19 pools. We trust the data and note the discrepancy
  (likely a different analysed subset). To be revisited.
- **Patient pairing:** ADC2 pairing is recoverable (84 patients, plasma+
  serum) from the `.xls` Study-Design metadata; **ADC1 pairing is not
  recoverable** from the REST/mwTab files (subject IDs blank). This is
  the reason per-biofluid modelling is chosen (one sample per patient
  within a biofluid → no pairing leakage).
- **Batch/temporal confound:** ADC1 (2014) vs ADC2 (2013) differ in run
  period → cross-cohort performance drop may partly reflect batch. This
  is reported as a finding, and normalisation is fit on ADC1 only.
- **Imbalance:** ADC1 balanced (86/86), ADC2 imbalanced (94/63) → use
  AUROC + AUPRC + balanced accuracy, never raw accuracy.

---

## 1b. Methodological rationale: two analyses, by question *and* by data regime

The case study runs **two analyses in parallel** because they answer
different questions — and because the right tool is also a function of
how much sample is available.

**Univariate inferential (differential analysis)** answers *which
metabolites change between conditions* — mechanism and hypothesis. It is
the only defensible approach in the typical **basic-research regime
(≈6–8 samples per condition)**, where multivariate prediction is not
licensed: n is too small to estimate a generalising classifier without
overfitting. Output: a ranked, FDR-controlled list with effect
directions.

**Multivariate predictive (ML / penalised regression)** answers *can we
combine metabolites to classify an unseen patient, and how well* —
diagnosis and clinical utility. It requires the **clinical-scale regime
(tens–hundreds of patients)** to be meaningful. This dataset (ADC1
86/86, ADC2 94/63) sits in that regime, which is precisely what makes
the predictive track legitimate here.

Key consequences to surface in the narrative:

- The two feature rankings will **partially disagree** by design: a
  strongly differential metabolite can have ~0 model importance if a
  correlated one already carries the signal; a univariately
  non-significant metabolite can be useful in combination. This is
  information, not contradiction.
- **Significance ≠ predictive value:** a large, significant mean
  difference can still classify poorly at the individual level if
  distributions overlap.
- At p≈n, **penalised regression (LASSO / elastic net)** is the bridge:
  a statistical model that is multivariate, predictive *and*
  interpretable. Tree ensembles (RF/XGBoost) are included only as
  comparators; they are not assumed to win at this n.

---

## 2. Methodology (the workflow that *is* the value)

Stack per user convention: **EDA + univariate inference in R**,
**predictive ML in Python**. Both tracks share Stage 1 and are
contrasted in the synthesis stage.

Modelling is run **per biofluid (plasma-only and serum-only,
separately)**, training on ADC1 and validating on ADC2, over the 138
shared metabolites. Within a biofluid each patient contributes a single
sample, so no patient-pairing grouping is required and pairing leakage
is structurally avoided.

### Stage 1 — EDA & QC (R, shared)
- Load identified-metabolite tables for both cohorts; profile
  dimensions, missingness, value distributions.
- PCA / hierarchical clustering coloured by case/control **and** by
  confounders (sex, smoking, organ, cohort/batch). Explicitly test
  whether cohort (2013 vs 2014) or organ drives separation.
- Document every preprocessing decision (log transform, scaling,
  imputation, filtering) with rationale; fit all transforms on ADC1 only.

### Track A — Univariate inferential / differential (R)
- Per-metabolite test (moderated t / Mann–Whitney as appropriate) for
  Adenocarcinoma vs Healthy, **within ADC1**, per biofluid.
- Multiple-testing control via Benjamini–Hochberg FDR.
- Output: ranked differential list with effect size + direction +
  q-value; volcano plot. Optional: check direction reproduces in ADC2.
- Purpose: mechanism / biological interpretation — *what changes*.

### Track B — Multivariate predictive (Python)
- **Nested cross-validation within ADC1** (outer = unbiased estimate,
  inner = hyperparameter + feature-selection tuning).
- **Feature selection strictly inside each CV fold** — central
  anti-leakage control.
- Models: penalised logistic regression (LASSO / elastic net) as the
  primary interpretable model; random forest / XGBoost as comparators
  only (not assumed to win at this n).
- Metrics: AUROC, AUPRC, balanced accuracy, calibration — with
  confidence intervals, not point estimates.

### Stage 3 — External validation (Track B)
- Lock the final per-biofluid pipeline on ADC1; evaluate **once** on
  ADC2. Report the honest performance drop and attribute it (biology vs
  batch).

### Stage 4 — Interpretability, panel & synthesis
- SHAP + LASSO coefficients on the final model; deliver a compact,
  defensible biomarker panel (not a 138-feature black box).
- **Synthesis (the headline result):** overlay Track A (differential)
  vs Track B (predictive) rankings; show explicitly where they agree and
  disagree, and explain why (redundancy, conditional effects,
  significance ≠ prediction). Cross-check both against literature.

### Stage 5 — Critical read & recommendation
- State limitations plainly (sample size, single-project
  generalisability, batch confound, metabolite coverage).
- Close with a business-framed recommendation: what each track supports,
  which biofluid classifies better, and what a confirmatory study needs.

---

## 3. Anti-pitfall checklist (credibility controls)

- [ ] Per-biofluid modelling (one sample per patient) → no pairing leakage.
- [ ] Cross-cohort features restricted to the 138 shared metabolites.
- [ ] Feature selection performed strictly inside CV folds.
- [ ] Normalisation/imputation fit on ADC1 only (never on ADC2).
- [ ] Nested CV for performance estimation (Track B).
- [ ] FDR control for differential tests (Track A).
- [ ] Metrics appropriate to class imbalance (AUPRC, balanced acc).
- [ ] Confounder check (cohort/batch, sex, smoking, organ) before claiming biology.
- [ ] External cohort (ADC2) touched **once**, at the end.
- [ ] Confidence intervals on all reported performance.
- [ ] Track A vs Track B agreement explicitly reported.
- [ ] Honest limitations section (incl. paper-vs-data discrepancy).

---

## 4. Deliverables

1. **Case study report** (`README` + short report): the narrative —
   problem, decisions, results, interpretation, recommendation. This is
   the primary artefact for the pharma audience.
2. **Reproducible repo**: R (EDA + differential) + Python (predictive),
   fixed seeds, environment files, one-command rerun.
3. **Figures**: PCA/QC coloured by cohort & organ; volcano (Track A);
   nested-CV + external-validation performance per biofluid (Track B);
   SHAP/LASSO panel; **Track A vs Track B agreement plot**.

### Proposed repo structure
```
ml-biomarker-lungADC/
├── README.md                  # the case study narrative
├── data/                      # download instructions only (no raw data committed)
├── R/                         # Stage 1 EDA/QC + Track A differential
│   ├── 01_eda_qc.R
│   └── 02_differential.R
├── python/                    # Track B predictive (Stages 3–4)
│   ├── 03_nested_cv.py
│   ├── 04_external_validation.py
│   └── 05_shap_interpret.py
├── figures/
├── report/                    # short PDF/HTML write-up + A-vs-B synthesis
├── environment.yml
└── renv.lock
```

---

## 5. Data availability

| Resource | Link |
|---|---|
| ADC1 — ST000368 | https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000368 |
| ADC2 — ST000369 | https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000369 |
| DOI (PR000293) | https://doi.org/10.21228/M85P57 |

Working files in repo `rawdata/`: `MSdata_ST000368_1.txt` (ADC1),
`MSdata_ST000369_1.txt` (ADC2), plus `*_AN0006xx.txt/.json` mwTab.

---

## 6. Risks & limitations (honest)

- **Single-project generalisability:** both cohorts come from one project
  (PR000293, one lab). External validation across cohorts is strong but
  not multi-site; do not overclaim clinical readiness.
- **Batch/temporal confound:** ADC1 (2014) vs ADC2 (2013) differ in run
  period; cross-cohort performance drop may reflect batch, not biology.
  Reported as a finding; normalisation fit on ADC1 only.
- **ADC1 patient pairing unavailable:** mitigated by per-biofluid
  modelling, but means a paired plasma–serum analysis within ADC1 is not
  possible without recovering ADC1 subject metadata.
- **Paper-vs-data discrepancy:** deposited ADC1 (86/86) does not match
  the published 52/31. Flagged for investigation; analysis follows the
  data.
- **Feature-selection instability:** at p≈n the predictive panel may
  vary across resamples; reported via selection-frequency, not as a
  single definitive list.

---

## 7. Next steps

- [x] Recover ADC1 (ST000368) + ADC2 (ST000369) processed data; verify
      composition, metabolite overlap (138 shared), batch dates.
- [ ] Build the 138-shared-metabolite, per-biofluid analysis matrices
      (ADC1 train / ADC2 validation), Adenocarcinoma vs Healthy.
- [ ] Stage 1 EDA in R; test cohort/organ as drivers of variance.
- [ ] Track A: differential analysis + FDR + volcano (R).
- [ ] Track B: nested-CV LASSO/elastic-net (+ RF/XGBoost comparators),
      in-fold selection (Python).
- [ ] Single external evaluation on ADC2, per biofluid.
- [ ] Synthesis: Track A vs Track B agreement; panel; draft narrative.

---

*Document generated: junio 2026; revised junio 2026 after data
verification (dual-cohort recovered; dual-track design added).*
*Scope decided with user: ML (Line B) + univariate inference, pharma
audience, biomarker discovery, Adenocarcinoma vs Healthy, per-biofluid
ADC1→ADC2 on 138 shared metabolites.*
