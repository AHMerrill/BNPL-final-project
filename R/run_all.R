# -----------------------------------------------------------------------------
# run_all.R  --  master script: runs the entire analysis end-to-end
# Zan Merrill, Andy Kim, Thomas Garner  |  OM 386 Spring 2026
#
# Usage from anywhere:
#   Rscript /path/to/repo/R/run_all.R
# Or from R / RStudio:
#   source("/path/to/repo/R/run_all.R")
#
# The script self-locates and sets the working directory to the repository
# root before sourcing the pipeline scripts, so it runs the same way
# regardless of where you launch R from.
#
# No analysis code is duplicated here; each pipeline script remains the
# single source of truth for its stage. If you edit 02_cluster.R, you do
# not edit this file too.
# -----------------------------------------------------------------------------

# ---------- self-locate ------------------------------------------------------
# Find the path to this script whether we were invoked via Rscript or
# source()'d from an interactive R session.
get.this.script.path = function() {
  # (1) Rscript: --file=... appears in command args
  args = commandArgs(trailingOnly = FALSE)
  match.idx = grep("^--file=", args)
  if(length(match.idx) > 0) {
    raw.path = sub("^--file=", "", args[match.idx[1]])
    # macOS Rscript escapes spaces in paths as "~+~"; reverse that here
    # so paths containing spaces (e.g. "marketing analytics") resolve.
    raw.path = gsub("~\\+~", " ", raw.path)
    return(normalizePath(raw.path))
  }
  # (2) source() from R: walk frames looking for $ofile
  for(frame.idx in seq_len(sys.nframe())) {
    ofile = sys.frame(frame.idx)$ofile
    if(!is.null(ofile)) return(normalizePath(ofile))
  }
  # (3) RStudio Source button uses rstudioapi
  if(requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p = tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if(nzchar(p)) return(normalizePath(p))
  }
  return(NA_character_)
}

script.path = get.this.script.path()
if(is.na(script.path)) {
  stop("run_all.R could not locate itself. Run via 'Rscript R/run_all.R' from ",
       "the repo root, or set your working directory to the repo root and ",
       "source('R/run_all.R').")
}

# script lives in <repo>/R/run_all.R  ->  repo root is two dirnames up
repo.root = dirname(dirname(script.path))
setwd(repo.root)
cat("Working directory:", getwd(), "\n")

# ---------- run the pipeline -------------------------------------------------
t.start = Sys.time()

cat("\n=== Step 1: PCA construct scores =============================\n")
source("R/01_pca.R")

cat("\n=== EDA: summary stats, correlations, item distributions ====\n")
source("R/00_eda.R")

cat("\n=== Step 2: clustering (mclust + hand-coded EM) =============\n")
source("R/02_cluster.R")

cat("\n=== Step 3: Bayesian hierarchical Gibbs + lmer cross-check ==\n")
source("R/03_bayes_hlm.R")

t.elapsed = difftime(Sys.time(), t.start, units = "secs")
cat(sprintf("\n=== Done. Total wall time: %.1f seconds ===\n",
            as.numeric(t.elapsed)))
cat("Outputs are in data/processed/, output/tables/, ",
    "output/figures/, and output/diagnostics/.\n", sep = "")
