# Overleaf upload bundle

Everything Overleaf needs in one flat folder. Drop all of these files into
the root of your Overleaf project (or zip the folder and use Overleaf's
"Upload Project" option).

## Contents

- `main.tex` — the report source, with `\graphicspath{{./}}` so it finds
  the figures in the same directory as the `.tex` file.
- `references.bib` — bibliography (the `dadra2024tep` key and a few
  methodology references).
- 12 figure PDFs referenced by `main.tex`:

```
01_pca_variance_explained.pdf
00_item_response_histograms.pdf
00_construct_correlation_heatmap.pdf
00_construct_pairs.pdf
02_mclust_bic.pdf
02_handem_diagnostics.pdf
02_segment_profiles_winner.pdf
03_trace_plots.pdf
03_trace_tau.pdf
03_posterior_densities.pdf
03_forest_by_segment.pdf
03_shrinkage_illustration.pdf
```

## If citations render as `?`

Two things must be true in Overleaf for `\citep{dadra2024tep}` to resolve:

1. `references.bib` is uploaded alongside `main.tex`.
2. BibTeX has actually run. Overleaf does this automatically when it sees
   a `.bib` file, but the first compile after upload may still show `?`.
   Fix: in Overleaf, click the small arrow next to "Recompile" and choose
   **Recompile from scratch**. After that, citations should resolve.

The compiler should be set to **pdfLaTeX** (Overleaf default) and the
bibliography engine to **BibTeX** (not biber) since `main.tex` uses
`natbib` + `\bibliographystyle{plainnat}`.

## Regenerating this bundle

From the repo root:

```bash
cp output/figures/*.pdf   output/for_overleaf/
cp report/references.bib  output/for_overleaf/
cp report/main.tex        output/for_overleaf/main.tex
# flatten graphicspath so images resolve in Overleaf's flat layout:
sed -i '' 's|\\graphicspath{{\\.\\./output/figures/}}|\\graphicspath{{./}}|' \
    output/for_overleaf/main.tex
```
