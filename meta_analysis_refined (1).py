"""meta_analysis_refined.py

Publication-ready random-effects meta-analysis for odds ratios.

Input formats supported (case-insensitive):
  A) logOR + SE
     - columns: logor (or log_or, logOR) and se (or std_error, SE)
  B) OR + CI bounds
     - columns: or (or OR) and ci_lower/ci_upper (or lower_ci/upper_ci)

Outputs (written to --output-dir):
  - forest_plot.png
  - funnel_plot.png
  - baujat_plot.png
  - leave_one_out.png
  - meta_analysis_results.txt
  - meta_analysis_summary_table.csv

Usage:
  python meta_analysis_refined.py --input meta_cleaned_egypt.csv --output-dir results
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2
import statsmodels.api as sm
from statsmodels.stats.meta_analysis import combine_effects


# -----------------------------
# Helpers
# -----------------------------

def _norm_col(s: str) -> str:
    return (
        str(s)
        .strip()
        .lower()
        .replace(" ", "_")
        .replace("-", "_")
    )


def _first_existing(cols: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    colset = set(cols)
    for c in candidates:
        if c in colset:
            return c
    return None


@dataclass
class MetaInputs:
    logor: np.ndarray
    se: np.ndarray
    labels: np.ndarray
    df_meta: pd.DataFrame


def _build_meta_inputs(df: pd.DataFrame) -> MetaInputs:
    df = df.copy()
    # normalize columns
    df.columns = [_norm_col(c) for c in df.columns]

    # identify possible columns
    col_logor = _first_existing(df.columns, ["logor", "log_or", "log_odds_ratio", "log_or_value"]) 
    col_se = _first_existing(df.columns, ["se", "std_error", "stderr", "standard_error"]) 

    col_or = _first_existing(df.columns, ["or", "odds_ratio", "oddsratio"]) 
    col_lo = _first_existing(df.columns, ["ci_lower", "lower_ci", "ci_low", "lcl", "lower"]) 
    col_hi = _first_existing(df.columns, ["ci_upper", "upper_ci", "ci_high", "ucl", "upper"]) 

    # coerce numeric safely
    for c in [col_logor, col_se, col_or, col_lo, col_hi]:
        if c and c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    if col_logor and col_se:
        df["logor"] = df[col_logor]
        df["se"] = df[col_se]
    elif col_or and col_lo and col_hi:
        # Require positivity
        bad = (df[col_or] <= 0) | (df[col_lo] <= 0) | (df[col_hi] <= 0)
        df.loc[bad, [col_or, col_lo, col_hi]] = np.nan
        df["logor"] = np.log(df[col_or])
        df["se"] = (np.log(df[col_hi]) - np.log(df[col_lo])) / (2 * 1.96)
    else:
        raise ValueError(
            "Input must contain either (logOR + SE) or (OR + CI lower/upper). "
            "Examples: logor,se OR or,ci_lower,ci_upper."
        )

    df_meta = df.dropna(subset=["logor", "se"]).copy()
    if df_meta.empty:
        raise ValueError("No valid rows found after computing logOR and SE.")

    # labels
    if "label" in df_meta.columns:
        labels = df_meta["label"].astype(str).to_numpy()
    elif "gene" in df_meta.columns and "snp" in df_meta.columns:
        labels = (df_meta["gene"].astype(str) + " " + df_meta["snp"].astype(str)).to_numpy()
    else:
        labels = df_meta.index.astype(str).to_numpy()

    return MetaInputs(
        logor=df_meta["logor"].to_numpy(dtype=float),
        se=df_meta["se"].to_numpy(dtype=float),
        labels=labels,
        df_meta=df_meta,
    )


def _egger_test(logor: np.ndarray, se: np.ndarray) -> Tuple[float, float]:
    """Egger's regression test.

    Standard normal deviate: SND = effect/SE
    Predictor: precision = 1/SE

    Returns:
      (intercept_p_value, slope_p_value)
    """
    snd = logor / se
    precision = 1.0 / se
    X = sm.add_constant(precision)
    model = sm.OLS(snd, X).fit()
    p_intercept = float(model.pvalues[0])
    p_slope = float(model.pvalues[1])
    return p_intercept, p_slope


# -----------------------------
# Main analysis
# -----------------------------

def run_meta_analysis(
    input_path: Path,
    output_dir: Path,
    title: str = "Forest Plot of SNPs Associated with Type 2 Diabetes",
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path)
    meta = _build_meta_inputs(df)

    effects = meta.logor
    ses = meta.se
    labels = meta.labels

    # Random effects meta-analysis (DerSimonian-Laird style implemented in statsmodels)
    random = combine_effects(effects, ses, method_re="iterated")

    pooled_logor = float(random.effect[0] if isinstance(random.effect, np.ndarray) else random.effect)
    pooled_se = float(random.sd_eff[0] if isinstance(random.sd_eff, np.ndarray) else random.sd_eff)

    pooled_or = float(np.exp(pooled_logor))
    pooled_ci = np.exp([pooled_logor - 1.96 * pooled_se, pooled_logor + 1.96 * pooled_se])

    # Heterogeneity
    Q = float(random.q)
    df_Q = int(random.k - 1)
    p_Q = float(1 - chi2.cdf(Q, df_Q)) if df_Q > 0 else float("nan")
    I2 = float(max(0.0, float(random.i2)) * 100.0)
    tau2 = float(random.tau2)

    # Study ORs and CIs
    or_values = np.exp(effects)
    ci_low = np.exp(effects - 1.96 * ses)
    ci_high = np.exp(effects + 1.96 * ses)

    # -----------------
    # Forest plot
    # -----------------
    all_or = np.append(or_values, pooled_or)
    all_ci_low = np.append(ci_low, pooled_ci[0])
    all_ci_high = np.append(ci_high, pooled_ci[1])
    all_labels = np.append(labels, "Pooled (Random Effects)")
    or_ci_text = [
        f"{o:.2f} ({lo:.2f}-{hi:.2f})" for o, lo, hi in zip(all_or, all_ci_low, all_ci_high)
    ]

    fig, ax = plt.subplots(figsize=(9.5, len(all_labels) * 0.45 + 2.2))
    yvals = np.arange(len(all_labels))

    # studies
    ax.errorbar(
        or_values,
        yvals[:-1],
        xerr=[or_values - ci_low, ci_high - or_values],
        fmt="s",
        color="black",
        ecolor="black",
        elinewidth=2,
        capsize=4,
        label="Study OR (95% CI)",
    )

    # pooled
    ax.errorbar(
        pooled_or,
        yvals[-1],
        xerr=[[pooled_or - pooled_ci[0]], [pooled_ci[1] - pooled_or]],
        fmt="D",
        color="crimson",
        markersize=11,
        elinewidth=3,
        capsize=6,
        label="Pooled OR",
    )

    ax.axvline(1.0, color="red", linestyle="--", linewidth=2)
    ax.set_yticks(yvals)
    ax.set_yticklabels(all_labels, fontsize=11)
    ax.set_xscale("log")
    ax.set_xlabel("Odds Ratio (log scale)", fontsize=13)
    ax.set_title(title, fontsize=15, pad=12)

    # show OR text on the right
    xmin, xmax = ax.get_xlim()
    xtext = xmax / 1.02
    for i, txt in enumerate(or_ci_text):
        ax.text(
            xtext,
            i,
            txt,
            va="center",
            ha="right",
            fontsize=10,
            fontweight="bold" if i == len(all_labels) - 1 else "normal",
        )

    ax.legend(loc="upper left", fontsize=10)
    plt.tight_layout()
    plt.savefig(output_dir / "forest_plot.png", dpi=350)
    plt.close()

    # -----------------
    # Funnel plot
    # -----------------
    plt.figure(figsize=(7.2, 6.2))
    plt.scatter(effects, 1 / ses, s=75, alpha=0.8)
    plt.axvline(pooled_logor, linestyle="--", linewidth=2)
    plt.xlabel("Log Odds Ratio")
    plt.ylabel("1/SE (Precision)")
    plt.title("Funnel Plot")
    plt.tight_layout()
    plt.savefig(output_dir / "funnel_plot.png", dpi=350)
    plt.close()

    # -----------------
    # Baujat plot (approx): residual vs Q contribution
    # -----------------
    residuals = (effects - pooled_logor) / ses
    q_contrib = ((effects - pooled_logor) ** 2) / (ses**2)

    plt.figure(figsize=(8.2, 6.2))
    plt.scatter(residuals, q_contrib, s=80, alpha=0.85, edgecolor="k")
    for i, lab in enumerate(labels):
        plt.annotate(str(lab), (residuals[i], q_contrib[i]), fontsize=8)
    plt.xlabel("Standardized residual")
    plt.ylabel("Q contribution")
    plt.title("Baujat Plot")
    plt.tight_layout()
    plt.savefig(output_dir / "baujat_plot.png", dpi=350)
    plt.close()

    # -----------------
    # Leave-one-out
    # -----------------
    loo_or = []
    loo_lo = []
    loo_hi = []
    for i in range(len(effects)):
        eff = np.delete(effects, i)
        se_ = np.delete(ses, i)
        r = combine_effects(eff, se_, method_re="iterated")
        pl = float(r.effect[0] if isinstance(r.effect, np.ndarray) else r.effect)
        ps = float(r.sd_eff[0] if isinstance(r.sd_eff, np.ndarray) else r.sd_eff)
        loo_or.append(np.exp(pl))
        loo_lo.append(np.exp(pl - 1.96 * ps))
        loo_hi.append(np.exp(pl + 1.96 * ps))

    loo_or = np.array(loo_or)
    loo_lo = np.array(loo_lo)
    loo_hi = np.array(loo_hi)

    plt.figure(figsize=(10.2, max(5.5, len(labels) * 0.45)))
    yvals = np.arange(len(labels))
    plt.errorbar(
        loo_or,
        yvals,
        xerr=[loo_or - loo_lo, loo_hi - loo_or],
        fmt="o",
        elinewidth=2.0,
        capsize=3,
    )
    plt.axvline(pooled_or, linestyle="--", linewidth=2)
    plt.yticks(yvals, labels, fontsize=9)
    plt.xscale("log")
    plt.xlabel("Odds Ratio (log scale)")
    plt.title("Leave-One-Out Analysis")
    plt.tight_layout()
    plt.savefig(output_dir / "leave_one_out.png", dpi=350)
    plt.close()

    # -----------------
    # Summary table
    # -----------------
    out_tbl = meta.df_meta.copy()
    out_tbl["or"] = np.exp(out_tbl["logor"])
    out_tbl["ci_lower"] = np.exp(out_tbl["logor"] - 1.96 * out_tbl["se"])
    out_tbl["ci_upper"] = np.exp(out_tbl["logor"] + 1.96 * out_tbl["se"])

    preferred_cols = [c for c in ["gene", "snp", "label", "or", "ci_lower", "ci_upper", "logor", "se"] if c in out_tbl.columns]
    out_tbl[preferred_cols].to_csv(output_dir / "meta_analysis_summary_table.csv", index=False)

    # -----------------
    # Results text
    # -----------------
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    p_int, p_slope = _egger_test(effects, ses)

    with open(output_dir / "meta_analysis_results.txt", "w", encoding="utf-8") as f:
        f.write(f"=== Random Effects Meta-Analysis Summary ({ts}) ===\n")
        f.write(f"k (studies): {len(effects)}\n")
        f.write(f"Pooled logOR: {pooled_logor:.4f}\n")
        f.write(f"SE: {pooled_se:.4f}\n")
        f.write(f"Pooled OR: {pooled_or:.4f}\n")
        f.write(f"95% CI (OR): [{pooled_ci[0]:.4f}, {pooled_ci[1]:.4f}]\n")
        f.write("\n=== Heterogeneity ===\n")
        f.write(f"Cochran Q: {Q:.4f} (df={df_Q}), p={p_Q:.6f}\n")
        f.write(f"I^2: {I2:.2f}%\n")
        f.write(f"Tau^2: {tau2:.6f}\n")
        f.write("\n=== Publication Bias (Egger) ===\n")
        f.write(f"Intercept p-value: {p_int:.6f}\n")
        f.write(f"Slope p-value: {p_slope:.6f}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Random-effects meta-analysis for odds ratios")
    parser.add_argument("--input", type=Path, required=True, help="Input CSV file")
    parser.add_argument("--output-dir", type=Path, default=Path("results"), help="Output directory")
    parser.add_argument("--title", type=str, default="Forest Plot of SNPs Associated with Type 2 Diabetes")
    args = parser.parse_args()

    run_meta_analysis(args.input, args.output_dir, title=args.title)
    print(f"Meta-analysis complete. Outputs written to: {args.output_dir}")


if __name__ == "__main__":
    main()
