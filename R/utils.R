# -----------------------------------------------------------------------------
# utils.R  --  shared helpers for the BNPL project
# Zan Merrill, Andy Kim, Thomas Garner  |  OM 386 Spring 2026
# -----------------------------------------------------------------------------

# construct -> items map (8 predictors + 1 outcome)
bnpl.constructs = list(
  Flexibility           = c("FL1", "FL2", "FL3"),
  PerformanceExpectancy = c("PE1", "PE2", "PE3"),
  EffortExpectancy      = c("EE1", "EE2", "EE3", "EE4"),
  PerceivedAutonomy     = c("PA1", "PA2", "PA3"),
  PerceivedCompetence   = c("PC1", "PC2", "PC3"),
  PerceivedRelatedness  = c("PRE1", "PRE3"),
  PerceivedUsefulness   = c("PU1", "PU2", "PU3"),
  Attitude              = c("AT1", "AT2", "AT3", "AT4"),
  IntentionToAdopt      = c("IA1", "IA2", "IA3")
)

bnpl.predictors = setdiff(names(bnpl.constructs), "IntentionToAdopt")

# -----------------------------------------------------------------------------
# Color palettes. Consistent segment colors across every figure so the reader
# recognizes groups at a glance; a 9-color qualitative palette for constructs;
# a diverging red-white-blue palette for correlations and segment profiles.
# -----------------------------------------------------------------------------

# segments (in index order 1, 2, 3, ...; extend if K > 3)
seg.palette = c("#2E86AB",  # Seg 1 -- blue
                "#E63946",  # Seg 2 -- red
                "#F4A261",  # Seg 3 -- orange
                "#2A9D8F",  # Seg 4 -- teal
                "#6A4C93")  # Seg 5 -- purple

# the 9 constructs (a distinct qualitative palette; predictors + outcome)
construct.palette = c(
  Flexibility           = "#2E86AB",
  PerformanceExpectancy = "#E63946",
  EffortExpectancy      = "#F4A261",
  PerceivedAutonomy     = "#2A9D8F",
  PerceivedCompetence   = "#6A4C93",
  PerceivedRelatedness  = "#C1666B",
  PerceivedUsefulness   = "#5F7470",
  Attitude              = "#D4A373",
  IntentionToAdopt      = "#264653"
)

# diverging red-white-blue (correlations, segment profile heatmaps)
diverging.palette = colorRampPalette(c("#B63B3B","white","#2E86AB"))(101)

# ensure an output dir exists (idempotent)
ensure.dir = function(path) {
  if(!dir.exists(path)) dir.create(path, recursive=TRUE)
  invisible(path)
}

# write a PNG and PDF of the same figure via base R graphics
# usage: save.fig("name", {plot(...); abline(...)}, w=6, h=4)
save.fig = function(name, expr, w=6, h=4, dpi=150, fig.dir="output/figures") {
  ensure.dir(fig.dir)
  png(file.path(fig.dir, paste0(name, ".png")), width=w*dpi, height=h*dpi, res=dpi)
  eval(substitute(expr), envir=parent.frame())
  dev.off()
  pdf(file.path(fig.dir, paste0(name, ".pdf")), width=w, height=h)
  eval(substitute(expr), envir=parent.frame())
  dev.off()
  invisible(NULL)
}

# round + format numeric df for writing to CSV (avoid 17-digit noise in tables)
round.df = function(df, digits=3) {
  num.cols = sapply(df, is.numeric)
  df[, num.cols] = lapply(df[, num.cols, drop=FALSE], round, digits)
  df
}
