# Meta-analysis Project (R + Python)

This repository contains a **publication-ready meta-analysis pipeline** implemented in both **Python** and **R**.

## Quick start

### Python

```bash
pip install pandas numpy matplotlib scipy statsmodels
python meta_analysis_refined.py --input meta_cleaned_egypt.csv --output-dir results
```

### R

```r
install.packages(c("optparse", "dplyr", "readr", "janitor", "meta", "metafor", "ggplot2"))
```

```bash
Rscript meta_analysis_refined.R --input meta.csv --outdir results
```

## What you get

- Random-effects pooled estimate (OR)
- Heterogeneity: Q, I², tau²
- Publication bias: funnel plot + Egger/linreg
- Influence diagnostics: Baujat plot + leave-one-out
- Publication-ready figures (PNG)

## Recommended folder structure

```text
.
├─ data/
│  └─ meta.csv
├─ results/
├─ meta_analysis_refined.py
├─ meta_analysis_refined.R
├─ META_ANALYSIS_PYTHON.md
└─ META_ANALYSIS_R.md
```
