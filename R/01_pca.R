# -----------------------------------------------------------------------------
# 01_pca.R  --  Step 1: construct scores via PCA
# Zan Merrill, Andy Kim, Thomas Garner  |  OM 386 Spring 2026
#
# Input : data/raw/BNPL Intention to use.xlsx  (226 x 28 Likert items)
# Output: data/processed/construct_scores.csv  (226 x 9 PC1 scores)
#         output/tables/pca_loadings.csv       (PC1 loading per item)
#         output/tables/pca_variance_explained.csv
#         output/diagnostics/step1_pca_diagnostics.txt
# -----------------------------------------------------------------------------

library(readxl)
source("R/utils.R")

set.seed(1)

data.path = "data/raw/BNPL Intention to use.xlsx"
bnpl.df   = as.data.frame(read_excel(data.path))

nn = nrow(bnpl.df)
KK = length(bnpl.constructs)

scores = matrix(NA_real_, nrow=nn, ncol=KK)
colnames(scores) = names(bnpl.constructs)

loadings.long = data.frame(construct=character(), item=character(),
                           loading=numeric(), stringsAsFactors=FALSE)
varexp.tbl = data.frame(construct=character(), n.items=integer(),
                        pc1.varexp=numeric(), stringsAsFactors=FALSE)

diag.lines = c()
for(cname in names(bnpl.constructs)) {
  items = bnpl.constructs[[cname]]
  X     = as.matrix(bnpl.df[, items])

  # PCA on correlation matrix (center + scale)
  pc = prcomp(X, center=TRUE, scale.=TRUE)

  load1  = pc$rotation[, 1]
  score1 = pc$x[, 1]

  # orient loadings so PC1 points "high construct" direction
  if(mean(load1) < 0) {
    load1  = -load1
    score1 = -score1
  }

  scores[, cname] = score1

  varexp = (pc$sdev^2) / sum(pc$sdev^2)
  varexp.tbl = rbind(varexp.tbl, data.frame(
    construct=cname, n.items=length(items), pc1.varexp=varexp[1]))
  loadings.long = rbind(loadings.long, data.frame(
    construct=cname, item=names(load1), loading=as.numeric(load1)))

  diag.lines = c(diag.lines,
    sprintf("=== %s (%d items) ===", cname, length(items)),
    sprintf("  PC1 variance explained: %.3f", varexp[1]),
    "  PC1 loadings:",
    paste0("    ", names(load1), ": ", sprintf("%.3f", load1)),
    "")
}

scores.df = as.data.frame(scores)

ensure.dir("data/processed")
ensure.dir("output/tables")
ensure.dir("output/diagnostics")

write.csv(scores.df, "data/processed/construct_scores.csv", row.names=FALSE)
write.csv(round.df(loadings.long), "output/tables/pca_loadings.csv", row.names=FALSE)
write.csv(round.df(varexp.tbl),    "output/tables/pca_variance_explained.csv", row.names=FALSE)
writeLines(diag.lines, "output/diagnostics/step1_pca_diagnostics.txt")

# visual: PC1 variance explained per construct (one color per construct)
save.fig("01_pca_variance_explained", {
  par(mar=c(9, 4, 3, 1))
  bar.cols = construct.palette[varexp.tbl$construct]
  barplot(varexp.tbl$pc1.varexp, names.arg=varexp.tbl$construct,
          las=2, col=bar.cols, border="white", ylim=c(0,1),
          ylab="PC1 variance explained",
          main="PC1 captures the majority of variance per construct")
  abline(h=0.5, lty=2, col="grey40", lwd=1.5)
}, w=8, h=5)

cat("Wrote data/processed/construct_scores.csv (", nrow(scores.df), "x", ncol(scores.df), ")\n")
cat("Wrote output/tables/pca_loadings.csv\n")
cat("Wrote output/tables/pca_variance_explained.csv\n")
cat("Wrote output/diagnostics/step1_pca_diagnostics.txt\n")
cat("Wrote output/figures/01_pca_variance_explained.{png,pdf}\n")
