# -----------------------------------------------------------------------------
# 00_eda.R  --  Exploratory data analysis for the BNPL project
# Zan Merrill, Andy Kim, Thomas Garner  |  OM 386 Spring 2026
#
# Input : data/raw/BNPL Intention to use.xlsx
#         data/processed/construct_scores.csv  (run 01_pca.R first)
# Output: output/tables/item_summary_stats.csv
#         output/tables/construct_summary_stats.csv
#         output/tables/construct_correlations.csv
#         output/figures/00_item_response_histograms.{png,pdf}
#         output/figures/00_construct_correlation_heatmap.{png,pdf}
#         output/figures/00_construct_pairs.{png,pdf}
# -----------------------------------------------------------------------------

library(readxl)
source("R/utils.R")

set.seed(1)

raw.df    = as.data.frame(read_excel("data/raw/BNPL Intention to use.xlsx"))
scores.df = read.csv("data/processed/construct_scores.csv")

ensure.dir("output/tables")
ensure.dir("output/figures")

# ---------- item-level summary stats (28 Likert items) -----------------------
all.items = unlist(bnpl.constructs, use.names=FALSE)
item.stats = data.frame(
  item   = all.items,
  n      = sapply(all.items, function(it) sum(!is.na(raw.df[[it]]))),
  mean   = sapply(all.items, function(it) mean(raw.df[[it]], na.rm=TRUE)),
  sd     = sapply(all.items, function(it) sd(raw.df[[it]], na.rm=TRUE)),
  min    = sapply(all.items, function(it) min(raw.df[[it]], na.rm=TRUE)),
  median = sapply(all.items, function(it) median(raw.df[[it]], na.rm=TRUE)),
  max    = sapply(all.items, function(it) max(raw.df[[it]], na.rm=TRUE))
)
# attach construct name
item.stats$construct = unlist(lapply(names(bnpl.constructs),
                              function(cn) rep(cn, length(bnpl.constructs[[cn]]))))
item.stats = item.stats[, c("construct","item","n","mean","sd","min","median","max")]
write.csv(round.df(item.stats), "output/tables/item_summary_stats.csv", row.names=FALSE)

# ---------- construct-level summary stats (9 PC1 scores) ---------------------
construct.stats = data.frame(
  construct = colnames(scores.df),
  n         = sapply(scores.df, function(x) sum(!is.na(x))),
  mean      = sapply(scores.df, mean),
  sd        = sapply(scores.df, sd),
  min       = sapply(scores.df, min),
  median    = sapply(scores.df, median),
  max       = sapply(scores.df, max)
)
write.csv(round.df(construct.stats), "output/tables/construct_summary_stats.csv", row.names=FALSE)

# ---------- construct correlation matrix ------------------------------------
cor.mat = cor(scores.df)
write.csv(round.df(as.data.frame(cor.mat), 3),
          "output/tables/construct_correlations.csv", row.names=TRUE)

# ---------- figure: item response distributions (faceted histograms) --------
# color each item's histogram by its parent construct so related items cluster
# visually
item.to.construct = setNames(
  unlist(lapply(names(bnpl.constructs),
                function(cn) rep(cn, length(bnpl.constructs[[cn]])))),
  all.items)

save.fig("00_item_response_histograms", {
  par(mfrow=c(5,6), mar=c(2,2,2,1), oma=c(1,1,2,1))
  for(it in all.items) {
    hist(raw.df[[it]], breaks=seq(0.5, 7.5, by=1),
         main=it, xlab="", ylab="",
         col=construct.palette[item.to.construct[it]],
         border="white", xaxt="n")
    axis(1, at=1:7, cex.axis=0.7)
  }
  mtext("Item response distributions (color = parent construct, 7-point Likert)",
        outer=TRUE, cex=1.1)
}, w=11, h=9)

# ---------- figure: construct correlation heatmap ---------------------------
save.fig("00_construct_correlation_heatmap", {
  par(mar=c(10, 10, 3, 2))
  rng = c(-1, 1)
  cm  = cor.mat[nrow(cor.mat):1, ]  # flip so first construct visually top-left
  image(1:ncol(cm), 1:nrow(cm), t(cm), col=diverging.palette, zlim=rng,
        xaxt="n", yaxt="n", xlab="", ylab="",
        main="Correlations between construct scores (PC1)")
  axis(1, at=1:ncol(cm), labels=colnames(cm), las=2, cex.axis=0.8)
  axis(2, at=1:nrow(cm), labels=rownames(cm), las=1, cex.axis=0.8)
  for(i in 1:nrow(cm)) for(j in 1:ncol(cm)) {
    text(j, i, sprintf("%.2f", cm[i,j]), cex=0.65,
         col=ifelse(abs(cm[i,j]) > 0.6, "white", "black"))
  }
}, w=8, h=7)

# ---------- figure: pairs plot of the 8 predictors + outcome ---------------
save.fig("00_construct_pairs", {
  pairs(scores.df, pch=20, cex=0.4,
        col=adjustcolor("#2E86AB", alpha.f=0.55),
        main="Pairwise construct-score relationships")
}, w=10, h=10)

cat("Wrote output/tables/item_summary_stats.csv\n")
cat("Wrote output/tables/construct_summary_stats.csv\n")
cat("Wrote output/tables/construct_correlations.csv\n")
cat("Wrote output/figures/00_item_response_histograms.{png,pdf}\n")
cat("Wrote output/figures/00_construct_correlation_heatmap.{png,pdf}\n")
cat("Wrote output/figures/00_construct_pairs.{png,pdf}\n")

show.table(construct.stats,                      "Construct-score summary statistics")
show.table(as.data.frame(round(cor.mat, 2)),     "Construct-score correlation matrix")
