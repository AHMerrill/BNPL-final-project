# BNPL Adoption: A Bayesian Hierarchical Take on TEP-Based Drivers

**OM 386 — Advanced Analytics in Marketing — Spring 2026**
**Team:** Zan Merrill, Andy Kim, Thomas Garner

## What this project is

We extend [Beyond Credit Cards: A TEP Perspective on Buy-Now-Pay-Later Adoption](references/TEP_BNPL_paper.pdf)
using a three-stage pipeline on the same 226-respondent survey:

1. **PCA per construct** — reduce 28 Likert items to 9 construct scores (PC1 per construct).
2. **Model-based clustering** — discover latent consumer segments on the 8 predictor constructs.
   Four methods compared (k-means, hierarchical Ward, `mclust` GMM, hand-rolled EM); winner
   picked by BIC + segment interpretability.
3. **Bayesian Hierarchical Linear Model** — regress Intention to Adopt on the 8 predictors with
   segment as the grouping variable. Hand-rolled Gibbs sampler, cross-validated against `lme4::lmer`.

This differs from the original PLS-SEM paper in three ways: (a) PCA-weighted composites instead
of equal-weight averaging, (b) mixture-based heterogeneity instead of one-size-fits-all
coefficients, and (c) full posterior distributions instead of p-values.

## Repo layout

```
.
├── R/
│   ├── utils.R              shared helpers
│   ├── 00_eda.R             summary stats, correlations, construct distributions
│   ├── 01_pca.R             Step 1 — construct scores via PCA
│   ├── 02_cluster.R         Step 2 — 4-method clustering comparison
│   └── 03_bayes_hlm.R       Step 3 — hand-rolled hierarchical Gibbs + lmer validation
├── data/
│   ├── raw/                 BNPL Intention to use.xlsx
│   └── processed/           construct_scores.csv, segment_labels.csv
├── output/
│   ├── figures/             PNG + PDF for report
│   ├── tables/              CSVs for report tables
│   └── diagnostics/         MCMC traces, EM convergence, model comparison
├── report/
│   └── main.tex             Overleaf-synced LaTeX source (source of truth: Overleaf)
├── references/              source paper, rubric, measurement items
└── README.md
```

## Reproducing

```r
# from repo root in RStudio:
source("R/00_eda.R")        # EDA tables + figures
source("R/01_pca.R")        # -> data/processed/construct_scores.csv
source("R/02_cluster.R")    # -> data/processed/segment_labels.csv
source("R/03_bayes_hlm.R")  # -> output/tables/posterior_summary.csv
```

Required packages: `readxl`, `mclust`, `mvtnorm`, `lme4`, `cluster`.

## Write-up

Source of truth for the paper is Overleaf — **full edit link:**
https://www.overleaf.com/4987795324cwjxkfyqwybs#ff0f4a

A local mirror lives at [report/main.tex](report/main.tex); figures and tables referenced from
the paper come from `output/figures/` and `output/tables/`.
