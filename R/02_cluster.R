# -----------------------------------------------------------------------------
# 02_cluster.R  --  Step 2: discover consumer segments via finite mixture model
# Zan Merrill, Andy Kim, Thomas Garner  |  OM 386 Spring 2026
#
# The course formally taught one clustering approach -- finite Gaussian mixture
# estimated via EM (Lecture 13, Beer_purchase_data example). We implement it
# two ways and compare:
#   (1) mclust::Mclust     -- "industrial" GMM; BIC picks K + covariance struct
#   (2) hand-coded EM      -- direct generalization of the Lecture 13 code from
#                             a 1-d regression mixture to a p-dim unsupervised
#                             mixture on the 8 predictor construct scores
#
# Winner (lower BIC) is exported as data/processed/segment_labels.csv.
#
# Input : data/processed/construct_scores.csv
# Output: data/processed/segment_labels.csv
#         output/tables/cluster_comparison.csv
#         output/tables/segment_profiles_<method>.csv
#         output/figures/02_*.{png,pdf}
#         output/diagnostics/step2_cluster_diagnostics.txt
# -----------------------------------------------------------------------------

library(mclust)
library(mvtnorm)
source("R/utils.R")

set.seed(1)

ensure.dir("output/figures")
ensure.dir("output/tables")
ensure.dir("output/diagnostics")
ensure.dir("data/processed")

# ---------- data --------------------------------------------------------------
scores.df = read.csv("data/processed/construct_scores.csv")
Y         = as.matrix(scores.df[, bnpl.predictors])  # 226 x 8
nn        = nrow(Y)
pp        = ncol(Y)

# PC1 scores are mean-centered; standardize to unit variance so no single
# predictor dominates the mixture covariance.
Y = scale(Y)

diag.lines = c(sprintf("N = %d respondents, P = %d predictors", nn, pp),
               sprintf("Predictors: %s", paste(bnpl.predictors, collapse=", ")),
               "")


###############################################################################
# (1) mclust::Mclust -- BIC picks K and covariance structure
###############################################################################

mc.fit   = Mclust(Y, G=1:9, verbose=FALSE)
K        = mc.fit$G                # chosen K; used for the hand-EM below too
mc.lab   = mc.fit$classification
mc.bic   = mc.fit$bic              # mclust sign convention: higher = better
mc.model = mc.fit$modelName

diag.lines = c(diag.lines,
  "=== (1) mclust (GMM with BIC search) ===",
  sprintf("  best K: %d   covariance model: %s", K, mc.model),
  sprintf("  BIC (mclust sign, higher=better): %.2f", mc.bic),
  sprintf("  BIC (standard sign, lower=better): %.2f", -mc.bic),
  sprintf("  segment sizes: %s", paste(table(mc.lab), collapse=", ")),
  "")

save.fig("02_mclust_bic", {
  par(mar=c(4,4,3,1))
  plot(mc.fit, what="BIC", legendArgs=list(x="bottomright", cex=0.6))
}, w=9, h=6)


###############################################################################
# (2) Hand-coded EM  --  mirrors Lecture 13 Beer_purchase_data structure
#     Generalization: K components (not 2), multivariate density (not 1-d),
#     clustering on Y (not regression on y|x). Same skeleton: initialize
#     parameters, iterate E-step and M-step, check LL convergence with
#     LL - prevLL < 10e-8, extract hard labels at the end.
###############################################################################

# initialize with k-means seeds (prof used hand-picked scalars for 1-d; this
# is the natural generalization for an 8-d feature matrix)
km.init = kmeans(Y, centers=K, nstart=20)
p.k     = as.numeric(table(km.init$cluster)) / nn
mu.k    = t(sapply(1:K, function(k) colMeans(Y[km.init$cluster == k, , drop=FALSE])))
Sigma.k = lapply(1:K, function(k) cov(Y[km.init$cluster == k, , drop=FALSE]) + diag(1e-4, pp))

gamma.ik = matrix(0, nn, K)
prevLL   = -Inf
converge = F
LL.trace = c()

for(i in 1:1000) {
  if(!converge && i < 1000) {

    #E-step
    log.dens = sapply(1:K, function(k)
      log(p.k[k]) + dmvnorm(Y, mean=mu.k[k, ], sigma=Sigma.k[[k]], log=TRUE))
    row.max  = apply(log.dens, 1, max)
    log.norm = row.max + log(rowSums(exp(log.dens - row.max)))   # log-sum-exp
    gamma.ik = exp(log.dens - log.norm)

    #M-step
    Nk   = colSums(gamma.ik)
    p.k  = Nk / nn
    for(k in 1:K) {
      mu.k[k, ]    = colSums(gamma.ik[, k] * Y) / Nk[k]
      diff.k       = sweep(Y, 2, mu.k[k, ], "-")
      Sigma.k[[k]] = (t(diff.k) %*% (gamma.ik[, k] * diff.k)) / Nk[k] + diag(1e-4, pp)
    }

    LL = sum(log.norm)
    LL.trace = c(LL.trace, LL)
    if(LL - prevLL < 10e-8) {
      converge = T
    }
    prevLL = LL
  } else {
    break
  }
}

#hard segment assignments (prof used the 1-d version:
# seg_ind = (gamma1>0.5)*1 + (gamma1<=0.5)*2)
em.lab = apply(gamma.ik, 1, which.max)

cat("Mixing proportions: ", round(p.k, 3), "\n")
cat("Segment sizes:      ", table(em.lab), "\n")
cat("Converged in:       ", i, "iterations  (LL =", round(LL, 2), ")\n")

# BIC in the standard sign convention (lower is better)
#   params: (K-1) mixing + K*p means + K*p(p+1)/2 full covariance
em.n.params = (K - 1) + K * pp + K * pp * (pp + 1) / 2
em.bic      = -2 * LL + em.n.params * log(nn)

diag.lines = c(diag.lines,
  "=== (2) hand-coded EM (multivariate Gaussian mixture) ===",
  sprintf("  K = %d (same as mclust's BIC-selected K)", K),
  sprintf("  iterations to convergence: %d", i),
  sprintf("  log-likelihood: %.2f", LL),
  sprintf("  parameters: %d  ->  BIC = %.2f (standard sign, lower=better)",
          em.n.params, em.bic),
  sprintf("  segment sizes: %s", paste(table(em.lab), collapse=", ")),
  sprintf("  ARI vs mclust labels: %.3f", adjustedRandIndex(em.lab, mc.lab)),
  "")

save.fig("02_handem_diagnostics", {
  par(mfrow=c(1,2), mar=c(4,4,3,1))
  plot(LL.trace, type="l", col="#264653", lwd=2,
       xlab="EM iteration", ylab="Log-likelihood",
       main=sprintf("hand-EM log-likelihood (K=%d)", K))
  barplot(p.k, names.arg=paste0("Seg", 1:K),
          col=seg.palette[1:K], border="white",
          ylim=c(0,1), ylab="Mixing proportion",
          main="hand-EM mixing proportions")
}, w=9, h=4)


###############################################################################
# Method comparison and winner selection
###############################################################################

comp.tbl = data.frame(
  method        = c("mclust (GMM, BIC search)", "hand-EM (GMM, full covariance)"),
  best.K        = c(K, K),
  BIC           = c(-mc.bic, em.bic),           # lower = better
  ARI.vs.mclust = c(1, adjustedRandIndex(em.lab, mc.lab))
)
write.csv(round.df(comp.tbl, 3),
          "output/tables/cluster_comparison.csv", row.names=FALSE)

winner     = if(em.bic < -mc.bic) "hand-EM" else "mclust"
winner.lab = if(winner == "hand-EM") em.lab else mc.lab

diag.lines = c(diag.lines,
  "=== Winner ===",
  sprintf("  Method: %s   K = %d", winner, K),
  sprintf("  BIC  mclust=%.2f   hand-EM=%.2f   (lower=better)", -mc.bic, em.bic),
  "")


###############################################################################
# Segment profiles: mean standardized PC1 score per segment per method
###############################################################################

profile.matrix = function(Ymat, labels, method.name) {
  Kk = length(unique(labels))
  prof = t(sapply(1:Kk, function(k) colMeans(Ymat[labels == k, , drop=FALSE])))
  rownames(prof) = paste0("Seg", 1:Kk)
  colnames(prof) = colnames(Ymat)
  prof.df = cbind(segment=rownames(prof),
                  size=as.integer(table(labels)[1:Kk]),
                  as.data.frame(prof))
  write.csv(round.df(prof.df, 3),
            sprintf("output/tables/segment_profiles_%s.csv",
                    gsub("[^a-z0-9]","",tolower(method.name))),
            row.names=FALSE)
  prof
}

prof.mc = profile.matrix(Y, mc.lab, "mclust")
prof.em = profile.matrix(Y, em.lab, "handem")
prof.w  = profile.matrix(Y, winner.lab, paste0("winner_", winner))

save.fig("02_segment_profiles_winner", {
  par(mar=c(10, 5, 3, 2))
  rng = c(-max(abs(prof.w)), max(abs(prof.w)))
  image(1:ncol(prof.w), 1:nrow(prof.w), t(prof.w)[, nrow(prof.w):1],
        col=diverging.palette, zlim=rng,
        xaxt="n", yaxt="n", xlab="", ylab="",
        main=sprintf("Segment profiles (%s, K=%d): mean standardized score",
                     winner, K))
  axis(1, at=1:ncol(prof.w), labels=colnames(prof.w), las=2, cex.axis=0.85)
  axis(2, at=1:nrow(prof.w), labels=rev(rownames(prof.w)), las=1, cex.axis=0.95)
  for(r in 1:nrow(prof.w)) for(c in 1:ncol(prof.w)) {
    v = prof.w[nrow(prof.w) - r + 1, c]
    text(c, r, sprintf("%.2f", v), cex=0.7,
         col=ifelse(abs(v) > 0.8, "white", "black"))
  }
}, w=9, h=5)

size.tbl = data.frame(
  method    = c("mclust", "hand-EM"),
  K         = c(K, K),
  seg.sizes = c(paste(table(mc.lab), collapse=", "),
                paste(table(em.lab), collapse=", "))
)
write.csv(size.tbl, "output/tables/segment_sizes_per_method.csv", row.names=FALSE)


###############################################################################
# Export winner labels for Step 3
###############################################################################

seg.out = data.frame(
  respondent.id = seq_len(nn),
  segment       = winner.lab
)
write.csv(seg.out, "data/processed/segment_labels.csv", row.names=FALSE)

writeLines(diag.lines, "output/diagnostics/step2_cluster_diagnostics.txt")

cat("\n=== Winner:", winner, "  K =", K, "===\n")
print(table(winner.lab))
cat("\nWrote data/processed/segment_labels.csv\n")
cat("Wrote output/tables/cluster_comparison.csv\n")
cat("Wrote output/tables/segment_profiles_*.csv\n")
cat("Wrote output/tables/segment_sizes_per_method.csv\n")
cat("Wrote output/figures/02_*.{png,pdf}\n")
cat("Wrote output/diagnostics/step2_cluster_diagnostics.txt\n")
