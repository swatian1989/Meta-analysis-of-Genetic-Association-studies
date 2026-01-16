# meta_analysis_refined.R
#
# Publication-ready random/fixed-effects meta-analysis for odds ratios.
#
# Supports input with either:
#   A) logOR + SE columns
#   B) OR + CI bounds columns
#   C) OR + 95% CI in a single string column like "1.20-1.80" or "1.20–1.80"
#
# Outputs (written to --outdir):
#   - meta_cleaned.csv
#   - forest_plot.png
#   - funnel_plot.png
#   - baujat_plot.png
#   - leave_one_out.png
#   - meta_analysis_results.txt
#
# Usage:
#   Rscript meta_analysis_refined.R --input meta.csv --outdir results

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(janitor)
  library(meta)
  library(metafor)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input CSV file"),
  make_option(c("-o", "--outdir"), type="character", default="results", help="Output directory"),
  make_option(c("--ci_col"), type="character", default="x95_percent_ci", help="Optional CI string column"),
  make_option(c("--title"), type="character", default="Forest Plot", help="Plot title")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$input)) stop("--input is required")
outdir <- opt$outdir
if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)

# -----------------------------
# Load & normalize
# -----------------------------

df <- read_csv(opt$input, show_col_types = FALSE) |> janitor::clean_names()

# Trim spaces in character columns
if (nrow(df) > 0) {
  df <- df |> mutate(across(where(is.character), ~trimws(.x)))
}

# Helper to locate first existing column name
first_col <- function(df, candidates) {
  for (c in candidates) if (c %in% colnames(df)) return(c)
  return(NA_character_)
}

col_logor <- first_col(df, c("logor", "log_or", "log_odds_ratio"))
col_se    <- first_col(df, c("se", "std_error", "stderr", "standard_error"))
col_or    <- first_col(df, c("or", "odds_ratio", "oddsratio"))
col_lo    <- first_col(df, c("ci_lower", "lower_ci", "ci_low", "lcl", "lower"))
col_hi    <- first_col(df, c("ci_upper", "upper_ci", "ci_high", "ucl", "upper"))
col_ci    <- if (opt$ci_col %in% colnames(df)) opt$ci_col else NA_character_

# Coerce numeric where relevant
numify <- function(x) suppressWarnings(as.numeric(x))
if (!is.na(col_logor)) df[[col_logor]] <- numify(df[[col_logor]])
if (!is.na(col_se))    df[[col_se]]    <- numify(df[[col_se]])
if (!is.na(col_or))    df[[col_or]]    <- numify(df[[col_or]])
if (!is.na(col_lo))    df[[col_lo]]    <- numify(df[[col_lo]])
if (!is.na(col_hi))    df[[col_hi]]    <- numify(df[[col_hi]])

# Build logOR and SE
if (!is.na(col_logor) && !is.na(col_se)) {
  df <- df |> mutate(logor = .data[[col_logor]], se = .data[[col_se]])
} else if (!is.na(col_or) && !is.na(col_lo) && !is.na(col_hi)) {
  df <- df |> mutate(
    logor = log(.data[[col_or]]),
    se    = (log(.data[[col_hi]]) - log(.data[[col_lo]])) / (2 * 1.96)
  )
} else if (!is.na(col_or) && !is.na(col_ci)) {
  # Parse CI string like "1.2-1.8" (handles en-dash/em-dash/minus)
  df[[col_ci]] <- gsub("[–—−]", "-", df[[col_ci]])
  df[[col_ci]] <- gsub("\\s", "", df[[col_ci]])
  df <- df |>
    mutate(
      ci_lower = numify(sub("^([0-9.]+).*$", "\\1", .data[[col_ci]])),
      ci_upper = numify(sub(".*-([0-9.]+)$", "\\1", .data[[col_ci]])),
      logor    = log(.data[[col_or]]),
      se       = (log(ci_upper) - log(ci_lower)) / (2 * 1.96)
    )
} else {
  stop("Input must contain either (logOR+SE) or (OR+CI bounds) or (OR + CI string column).")
}

# Drop invalid rows
meta_df <- df |> filter(!is.na(logor), !is.na(se), se > 0)
if (nrow(meta_df) == 0) stop("No valid rows after computing logOR and SE.")

# Labels
label <- if ("label" %in% colnames(meta_df)) {
  meta_df$label
} else if (all(c("gene", "snp") %in% colnames(meta_df))) {
  paste(meta_df$gene, meta_df$snp)
} else {
  as.character(seq_len(nrow(meta_df)))
}

meta_df <- meta_df |> mutate(label = make.unique(as.character(label)))

# Save cleaned
write_csv(meta_df, file.path(outdir, "meta_cleaned.csv"))

# -----------------------------
# Meta-analysis
# -----------------------------

m <- metagen(
  TE = meta_df$logor,
  seTE = meta_df$se,
  studlab = meta_df$label,
  sm = "OR",
  comb.random = TRUE,
  comb.fixed = TRUE,
  hakn = TRUE
)

rma_obj <- rma(yi = meta_df$logor, sei = meta_df$se, method = "REML")

# -----------------------------
# Plots
# -----------------------------

# Forest plot (meta)
png(file.path(outdir, "forest_plot.png"), width=2600, height=3300, res=450)
forest(
  m,
  xlab = "Odds Ratio (log scale)",
  leftcols = c("studlab"),
  leftlabs = c("Study"),
  print.tau2 = TRUE,
  comb.random = TRUE,
  comb.fixed = TRUE
)
dev.off()

# Funnel plot
png(file.path(outdir, "funnel_plot.png"), width=1800, height=1800, res=300)
funnel(m, xlab = "Log(OR)", ylab = "Standard Error")
dev.off()

# Baujat plot (metafor)
png(file.path(outdir, "baujat_plot.png"), width=2200, height=1700, res=300)
try(baujat(rma_obj), silent=TRUE)
dev.off()

# Leave-one-out
loo <- leave1out(rma_obj)
loo_df <- data.frame(
  label = rownames(loo),
  estimate = loo$estimate,
  se = loo$se
)
loo_df <- loo_df |>
  mutate(
    or = exp(estimate),
    lo = exp(estimate - 1.96 * se),
    hi = exp(estimate + 1.96 * se)
  )

p <- ggplot(loo_df, aes(y = reorder(label, or), x = or)) +
  geom_point() +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.2) +
  geom_vline(xintercept = exp(coef(rma_obj)[1]), linetype = "dashed") +
  scale_x_log10() +
  labs(title = "Leave-One-Out Analysis", x = "Odds Ratio (log scale)", y = NULL) +
  theme_bw(base_size = 12)

ggsave(file.path(outdir, "leave_one_out.png"), plot = p, width = 9, height = max(5, 0.25*nrow(loo_df)+2), dpi = 300)

# -----------------------------
# Results summary
# -----------------------------

ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
pooled_or <- exp(m$TE.random)
ci_or <- exp(m$lower.random); ci_or_u <- exp(m$upper.random)

# Egger test (meta package)
# Note: requires >= 10 studies for best practice, but we still compute.
regtest_res <- try(metabias(m, method.bias = "linreg"), silent=TRUE)

sink(file.path(outdir, "meta_analysis_results.txt"))
cat("=== Random/Fixed Effects Meta-Analysis Summary (", ts, ") ===\n", sep="")
print(summary(m))
cat("\n=== Pooled (Random) OR ===\n")
cat(sprintf("Pooled OR: %.4f\n", pooled_or))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", ci_or, ci_or_u))
cat("\n=== Heterogeneity ===\n")
cat(sprintf("Q: %.4f (df=%d), p=%.6f\n", m$Q, m$df.Q, m$pval.Q))
cat(sprintf("I^2: %.2f%%\n", m$I2))
cat(sprintf("Tau^2: %.6f\n", m$tau^2))
cat("\n=== Publication Bias (Egger / Linreg) ===\n")
if (inherits(regtest_res, "try-error")) {
  cat("Egger test could not be computed.\n")
} else {
  print(regtest_res)
}
sink()

cat("Done. Outputs written to:", outdir, "\n")
