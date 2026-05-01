# BNPL Adoption: A Bayesian Hierarchical Take on TEP-Based Drivers

**OM 386 — Advanced Analytics in Marketing — Spring 2026**
**Team:** Andy Dohyung Kim, Thomas Garner, Zan Merrill

## What this project is

We extend Dadra, Sonavane, Bachwani & Behera (2024),
[_Beyond Credit Cards: A TEP Perspective on Buy-Now-Pay-Later Adoption_](https://www.tandfonline.com/doi/full/10.1080/21639159.2024.2391282),
_Journal of Global Scholars of Marketing Science_ 34(4), 459–484,
using a three-stage pipeline on the same 226-respondent survey:

1. **PCA per construct** — nine separate PCAs reduce 28 Likert items to
   9 construct scores (PC1 per construct).
2. **Model-based clustering** — discover latent consumer segments on the
   8 predictor constructs via finite Gaussian mixture estimated by EM
   (Lecture 13). Two implementations compared: `mclust::Mclust` (BIC
   selects K and covariance structure) and a hand-coded EM sampler that
   generalizes the Lecture 13 regression-mixture template to
   unsupervised multivariate clustering. Both pick K = 3.
3. **Bayesian Hierarchical Linear Model** — regress Intention to Adopt
   on the 8 predictors with segment-specific random coefficients,
   estimated via a hand-coded Gibbs sampler that extends
   `data examples/BayeianLM.r` to the hierarchical case.
   Cross-validated against `lme4::lmer`.

This differs from the original PLS-SEM paper in three ways: (a) PCA-weighted
composites instead of equal-weight averaging, (b) mixture-based segment
heterogeneity instead of one-size-fits-all coefficients, and (c) full
posterior distributions instead of p-values.

## Repo layout

```
.
├── R/
│   ├── utils.R                  shared helpers (palettes, save.fig, show.table)
│   ├── 00_eda.R                 summary stats, correlations, construct distributions
│   ├── 01_pca.R                 Step 1 — construct scores via PCA
│   ├── 02_cluster.R             Step 2 — mclust + hand-coded EM mixture
│   ├── 03_bayes_hlm.R           Step 3 — hand-coded hierarchical Gibbs + lmer cross-check
│   └── run_all.R                master script: sources the four above end-to-end
├── data/
│   ├── raw/                     BNPL Intention to use.xlsx
│   └── processed/               construct_scores.csv, segment_labels.csv
├── output/
│   ├── figures/                 PNG + PDF figures (reproduced by run_all.R)
│   ├── tables/                  CSVs (summary stats, posteriors, comparisons)
│   ├── diagnostics/             MCMC traces, EM convergence, model comparison
│   └── for_overleaf/            drop-in upload bundle for Overleaf
│                                (main.tex + references.bib + all figure PDFs)
└── README.md
```

## Reproducing

One command from anywhere:

```bash
Rscript /path/to/repo/R/run_all.R
```

The master script self-locates, sets the working directory to the repo
root, and sources the four pipeline scripts in dependency order. Outputs
land in `data/processed/`, `output/tables/`, `output/figures/`, and
`output/diagnostics/`. Total wall time on a recent laptop: ~6 seconds.

You can also source it from RStudio:

```r
source("/path/to/repo/R/run_all.R")
```

In an interactive RStudio session, the master script additionally
renders each figure to the Plots pane (multi-panel figures that don't
fit fall back to the file outputs without breaking the run) and prints
key tables to the console; the per-segment posterior summary is
opened in the data viewer via `View()`.

If you prefer to step through stage by stage:

```r
source("R/01_pca.R")        # -> data/processed/construct_scores.csv
source("R/00_eda.R")        # EDA tables + figures
source("R/02_cluster.R")    # -> data/processed/segment_labels.csv
source("R/03_bayes_hlm.R")  # -> output/tables/posterior_summary.csv
```

Required packages: `readxl`, `mclust`, `mvtnorm`, `MCMCpack`, `lme4`.

Outputs overwrite cleanly on every run; no manual cleanup is needed
between iterations.

## Write-up

Source of truth for the paper is Overleaf — **full edit link:**
https://www.overleaf.com/4987795324cwjxkfyqwybs#ff0f4a

The drop-in upload bundle for Overleaf is at
[output/for_overleaf/](output/for_overleaf/) — `main.tex`,
`references.bib`, and all 12 figure PDFs in one flat folder, with
`\graphicspath{{./}}` so the LaTeX builds in Overleaf's flat layout.
After uploading, choose **Recompile from scratch** in Overleaf's
recompile menu so the bibliography resolves on the first build.
