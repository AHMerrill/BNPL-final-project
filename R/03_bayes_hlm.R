# -----------------------------------------------------------------------------
# 03_bayes_hlm.R  --  Step 3: Bayesian Hierarchical Linear Model (Gibbs)
# Zan Merrill, Andy Kim, Thomas Garner  |  OM 386 Spring 2026
#
# Model: y_i = x_i' beta_{seg[i]} + eps_i,   eps ~ N(0, 1/tau)
#        beta_k ~ N(mu_beta, Sigma_beta)                  (random coefficients)
#        mu_beta ~ N(0, 1e6 * I)                          (diffuse hyperprior)
#        Sigma_beta ~ Inv-Wishart(nu0 = p+2, S0 = I)       (weakly-informative)
#        tau ~ Gamma(1/2, 1/2)
#
# The sampler directly generalizes the Lecture 3 / BayeianLM.r template from
# a single beta vector to K segment-specific beta_k with a Normal-Inverse-
# Wishart hyperprior. Same diagnostic conventions: ts.plot traces,
# hist + lines(density) posteriors, colMeans summaries.
#
# Baselines for comparison:
#   (a) complete-pooling OLS:  one beta for all N=226                     (glm)
#   (b) no-pooling OLS:        separate beta_k per segment, no pooling    (lm)
#   (c) partial-pooling Bayes: our Gibbs output                           (this)
#
# Additional frequentist cross-check via lme4::lmer(y ~ x + (1 | segment)).
#
# Input : data/processed/construct_scores.csv
#         data/processed/segment_labels.csv
# Output: output/tables/posterior_summary.csv
#         output/tables/coefficient_comparison.csv
#         output/figures/03_*.{png,pdf}
#         output/diagnostics/step3_gibbs_diagnostics.txt
# -----------------------------------------------------------------------------

library(mvtnorm)
library(MCMCpack)   # for riwish
library(lme4)
source("R/utils.R")

set.seed(1)

ensure.dir("output/figures")
ensure.dir("output/tables")
ensure.dir("output/diagnostics")

# ---------- data --------------------------------------------------------------
scores.df = read.csv("data/processed/construct_scores.csv")
seg.df    = read.csv("data/processed/segment_labels.csv")

stopifnot(nrow(scores.df) == nrow(seg.df))

# outcome + predictors (PC1 scores are already mean-centered)
Y   = as.numeric(scores.df$IntentionToAdopt)
Xp  = as.matrix(scores.df[, bnpl.predictors])
X   = cbind(1, Xp)  # intercept + 8 predictors -> design matrix 226 x 9
colnames(X)[1] = "(Intercept)"

seg   = seg.df$segment
K     = length(unique(seg))
n.obs = length(Y)
p.dim = ncol(X)

cat("N =", n.obs, " P =", p.dim, " (1 intercept + 8 predictors)", " K =", K, "\n")
cat("Segment sizes:", table(seg), "\n\n")


###############################################################################
# (a) Complete-pooling OLS  (no segment heterogeneity)
###############################################################################
pool.fit = lm(Y ~ Xp)
pool.beta = coef(pool.fit)
names(pool.beta) = colnames(X)
pool.se   = summary(pool.fit)$coefficients[, "Std. Error"]


###############################################################################
# (b) No-pooling OLS  (segment-stratified; no information shared)
###############################################################################
nop.beta = matrix(NA_real_, K, p.dim)
nop.se   = matrix(NA_real_, K, p.dim)
colnames(nop.beta) = colnames(X)
colnames(nop.se)   = colnames(X)
for(k in 1:K) {
  fit.k = lm(Y[seg == k] ~ Xp[seg == k, ])
  nop.beta[k, ] = coef(fit.k)
  nop.se[k, ]   = summary(fit.k)$coefficients[, "Std. Error"]
}


###############################################################################
# (c) Partial-pooling Bayesian HLM via Gibbs sampler
###############################################################################

# ---------- priors (diffuse / weakly-informative, prof's BayeianLM.r style)
mu.0    = rep(0, p.dim)
V.0     = 1E6 * diag(p.dim)     # huge -> diffuse prior on hyper-mean
iV.0    = 1E-6 * diag(p.dim)

nu.0    = p.dim + 2             # minimum df for proper mean
S.0     = diag(p.dim)           # weakly informative scale

a.tau   = 1/2                   # same as prof's BayeianLM.r
b.tau   = 1/2

# ---------- MCMC settings
NIT   = 10000
nBurn = 2000
nThin = 1                       # chain mixes well; no thinning needed

# ---------- posterior sample storage
beta.pos   = array(0, c(NIT, K, p.dim),
                   dimnames=list(NULL, paste0("Seg",1:K), colnames(X)))
mu.pos     = matrix(0, NIT, p.dim, dimnames=list(NULL, colnames(X)))
Sigma.pos  = array(0, c(NIT, p.dim, p.dim))
tau.pos    = rep(0, NIT)

# ---------- precompute per-segment X'X, X'Y for efficiency
X.k  = lapply(1:K, function(k) X[seg == k, , drop=FALSE])
y.k  = lapply(1:K, function(k) Y[seg == k])
XtX.k = lapply(X.k, function(Xk) t(Xk) %*% Xk)
XtY.k = lapply(1:K, function(k) t(X.k[[k]]) %*% y.k[[k]])
nk    = sapply(1:K, function(k) length(y.k[[k]]))

# ---------- initialize the loop
curBeta  = nop.beta            # warm-start at no-pool OLS estimates
curMu    = colMeans(curBeta)
curSigma = var(curBeta) + 0.01 * diag(p.dim)
curTau   = 1

cat("Starting Gibbs sampler: NIT =", NIT, " burn =", nBurn, "\n")
t.start = Sys.time()

# ---------- main Gibbs loop (mirrors BayeianLM.r structure)
for(m in 1:(NIT + nBurn)) {

  # step 1: sample beta_k | rest  (Normal-Normal conjugacy per segment)
  iSigma = solve(curSigma)
  for(k in 1:K) {
    V.k   = solve(curTau * XtX.k[[k]] + iSigma)
    m.k   = V.k %*% (curTau * XtY.k[[k]] + iSigma %*% curMu)
    curBeta[k, ] = as.vector(rmvnorm(1, mean=m.k, sigma=V.k))
  }

  # step 2: sample mu_beta | rest
  V.mu   = solve(iV.0 + K * iSigma)
  m.mu   = V.mu %*% (iV.0 %*% mu.0 + iSigma %*% colSums(curBeta))
  curMu  = as.vector(rmvnorm(1, mean=m.mu, sigma=V.mu))

  # step 3: sample Sigma_beta | rest  (Inverse-Wishart)
  beta.dev = sweep(curBeta, 2, curMu, "-")
  S.post   = S.0 + t(beta.dev) %*% beta.dev
  curSigma = riwish(nu.0 + K, S.post)

  # step 4: sample tau | rest  (Gamma)  -- residuals across all segments
  res = numeric(n.obs)
  for(k in 1:K) {
    res[seg == k] = y.k[[k]] - X.k[[k]] %*% curBeta[k, ]
  }
  curTau = rgamma(1, a.tau + 0.5 * n.obs,
                     b.tau + 0.5 * sum(res^2))

  # store after burn-in
  if(m > nBurn) {
    ix = m - nBurn
    beta.pos[ix, , ] = curBeta
    mu.pos[ix, ]     = curMu
    Sigma.pos[ix, , ]= curSigma
    tau.pos[ix]      = curTau
  }
}

t.elapsed = difftime(Sys.time(), t.start, units="secs")
cat("Gibbs done in", round(as.numeric(t.elapsed), 1), "seconds\n")


###############################################################################
# Posterior summaries
###############################################################################

# per-segment coefficient: posterior mean, SD, 5% / 50% / 95%
post.summary = expand.grid(segment=paste0("Seg", 1:K),
                           coef=colnames(X), stringsAsFactors=FALSE)
post.summary$post.mean = NA
post.summary$post.sd   = NA
post.summary$q05       = NA
post.summary$q50       = NA
post.summary$q95       = NA
post.summary$P.gt0     = NA

for(k in 1:K) for(j in 1:p.dim) {
  samp = beta.pos[, k, j]
  row  = which(post.summary$segment == paste0("Seg", k) &
               post.summary$coef    == colnames(X)[j])
  post.summary$post.mean[row] = mean(samp)
  post.summary$post.sd[row]   = sd(samp)
  post.summary$q05[row]       = quantile(samp, 0.05)
  post.summary$q50[row]       = quantile(samp, 0.50)
  post.summary$q95[row]       = quantile(samp, 0.95)
  post.summary$P.gt0[row]     = mean(samp > 0)
}
write.csv(round.df(post.summary, 3),
          "output/tables/posterior_summary.csv", row.names=FALSE)

# hyper-parameter summary (mu_beta, sigma = sqrt(1/tau))
mu.summary = data.frame(
  coef      = colnames(X),
  post.mean = colMeans(mu.pos),
  post.sd   = apply(mu.pos, 2, sd),
  q05       = apply(mu.pos, 2, quantile, 0.05),
  q95       = apply(mu.pos, 2, quantile, 0.95)
)
write.csv(round.df(mu.summary, 3),
          "output/tables/hyper_mean_summary.csv", row.names=FALSE)

resid.sd = sqrt(1 / tau.pos)
cat("\nResidual SD (posterior):  mean =", round(mean(resid.sd), 3),
    "   90% CI =", round(quantile(resid.sd, 0.05), 3),
    "-", round(quantile(resid.sd, 0.95), 3), "\n")


###############################################################################
# Comparison table: complete-pool / no-pool / partial-pool posteriors
###############################################################################

# reduce Bayesian beta_k to per-segment posterior means for the comparison
bayes.mean = apply(beta.pos, c(2,3), mean)       # K x p

comp.rows = list()
for(j in 1:p.dim) {
  for(k in 1:K) {
    comp.rows[[length(comp.rows)+1]] = data.frame(
      coef         = colnames(X)[j],
      segment      = paste0("Seg", k),
      complete.pool= pool.beta[j],            # same across segments
      no.pool      = nop.beta[k, j],
      partial.pool = bayes.mean[k, j]
    )
  }
}
comp.tbl = do.call(rbind, comp.rows)
write.csv(round.df(comp.tbl, 3),
          "output/tables/coefficient_comparison.csv", row.names=FALSE)


###############################################################################
# lme4 cross-check: random-intercept model  (full random slopes won't
# identify with only K=3 segments; we compare its fixed effects to our
# hyper-mean mu_beta)
###############################################################################

lmer.df = data.frame(Y=Y, segment=factor(seg), Xp)
lmer.form = as.formula(paste("Y ~", paste(bnpl.predictors, collapse=" + "),
                             "+ (1 | segment)"))
lmer.fit = lmer(lmer.form, data=lmer.df, REML=FALSE)
lmer.fixed = fixef(lmer.fit)

xcheck = data.frame(
  coef       = colnames(X),
  bayes.mu   = colMeans(mu.pos),
  lmer.fixed = lmer.fixed[colnames(X)]
)
xcheck$diff = xcheck$bayes.mu - xcheck$lmer.fixed
write.csv(round.df(xcheck, 3),
          "output/tables/lmer_xcheck.csv", row.names=FALSE)


###############################################################################
# Diagnostics / figures
###############################################################################

# trace plots: intercept per segment + one representative slope (Attitude)
save.fig("03_trace_plots", {
  par(mfrow=c(2, K), mar=c(4,4,3,1))
  att.col = which(colnames(X) == "Attitude")
  for(k in 1:K) {
    ts.plot(beta.pos[, k, 1], col=seg.palette[k],
            ylab="Intercept",
            main=sprintf("Seg %d: intercept trace", k))
  }
  for(k in 1:K) {
    ts.plot(beta.pos[, k, att.col], col=seg.palette[k],
            ylab="beta (Attitude)",
            main=sprintf("Seg %d: Attitude slope trace", k))
  }
}, w=10, h=6)

# trace for tau / residual SD
save.fig("03_trace_tau", {
  par(mfrow=c(1,2), mar=c(4,4,3,1))
  ts.plot(tau.pos, col="#264653", lwd=1.2,
          ylab="tau", main="tau (precision) trace")
  ts.plot(sqrt(1/tau.pos), col="#264653", lwd=1.2,
          ylab="Residual SD", main="Residual SD = sqrt(1/tau)")
}, w=9, h=4)

# posterior densities: per-segment coefficient for each predictor
# one panel per predictor, K overlaid densities
save.fig("03_posterior_densities", {
  par(mfrow=c(3,3), mar=c(4,3,3,1))
  for(j in 1:p.dim) {
    dens.list = lapply(1:K, function(k) density(beta.pos[, k, j]))
    x.rng = range(sapply(dens.list, function(d) d$x))
    y.rng = c(0, max(sapply(dens.list, function(d) d$y)))
    plot(NA, xlim=x.rng, ylim=y.rng,
         xlab="coefficient", ylab="",
         main=colnames(X)[j])
    # shaded density + line outline per segment
    for(k in 1:K) {
      dd = dens.list[[k]]
      polygon(c(dd$x, rev(dd$x)), c(dd$y, rep(0, length(dd$y))),
              col=adjustcolor(seg.palette[k], alpha.f=0.25), border=NA)
      lines(dd, col=seg.palette[k], lwd=2)
    }
    abline(v=0, lty=2, col="grey40")
    if(j == 1) legend("topright", paste0("Seg",1:K),
                      col=seg.palette[1:K], lwd=2, bty="n", cex=0.8)
  }
}, w=10, h=9)

# forest plot: per-segment beta with 90% CI for each predictor
save.fig("03_forest_by_segment", {
  par(mar=c(4, 12, 3, 1))
  pred.idx = which(colnames(X) != "(Intercept)")
  n.pred = length(pred.idx)
  y.pos = 1:(n.pred * K)
  plot(NA, xlim=c(-1.5, 1.5), ylim=c(0.5, max(y.pos) + 0.5),
       xlab="Posterior coefficient (90% credible interval)",
       ylab="", yaxt="n",
       main="Per-segment coefficients (Bayesian posterior)")
  # light horizontal bands per predictor for readability
  for(pp in seq_along(pred.idx)) {
    y.mid = (n.pred - pp) * K + (K+1)/2
    if(pp %% 2 == 0) {
      rect(par("usr")[1], y.mid - K/2, par("usr")[2], y.mid + K/2,
           col=adjustcolor("grey", alpha.f=0.1), border=NA)
    }
  }
  for(pp in seq_along(pred.idx)) {
    j = pred.idx[pp]
    for(k in 1:K) {
      y0 = (n.pred - pp) * K + k
      m  = mean(beta.pos[, k, j])
      lo = quantile(beta.pos[, k, j], 0.05)
      hi = quantile(beta.pos[, k, j], 0.95)
      segments(lo, y0, hi, y0, col=seg.palette[k], lwd=2)
      points(m, y0, pch=19, col=seg.palette[k], cex=1.1)
    }
  }
  abline(v=0, lty=2, col="grey40", lwd=1.2)
  axis(2, at=((n.pred - 1):0) * K + (K+1)/2,
       labels=colnames(X)[pred.idx], las=1, cex.axis=0.85, tick=FALSE)
  legend("bottomright", paste0("Seg",1:K),
         col=seg.palette[1:K], pch=19, lwd=2, bty="n", cex=0.85)
}, w=9, h=7)

###############################################################################
# Diagnostics text file
###############################################################################
diag.lines = c(
  sprintf("Gibbs sampler: NIT=%d, burn=%d, elapsed=%.1fs", NIT, nBurn, as.numeric(t.elapsed)),
  sprintf("N=%d, P=%d, K=%d segments", n.obs, p.dim, K),
  sprintf("Segment sizes: %s", paste(table(seg), collapse=", ")),
  "",
  "=== Posterior means by segment ===",
  "",
  capture.output(print(round(bayes.mean, 3))),
  "",
  "=== Hyper-mean mu_beta (posterior) vs lmer fixed effects ===",
  "",
  capture.output(print(round.df(xcheck, 3))),
  "",
  sprintf("Residual SD posterior mean: %.3f  (90%% CI: %.3f - %.3f)",
          mean(resid.sd), quantile(resid.sd, 0.05), quantile(resid.sd, 0.95)),
  "",
  "=== lme4 model summary ===",
  capture.output(print(summary(lmer.fit)))
)
writeLines(diag.lines, "output/diagnostics/step3_gibbs_diagnostics.txt")

cat("\nWrote output/tables/posterior_summary.csv\n")
cat("Wrote output/tables/coefficient_comparison.csv\n")
cat("Wrote output/tables/hyper_mean_summary.csv\n")
cat("Wrote output/tables/lmer_xcheck.csv\n")
cat("Wrote output/figures/03_*.{png,pdf}\n")
cat("Wrote output/diagnostics/step3_gibbs_diagnostics.txt\n")
