"""Static plots and an auditable HTML report for an analysis run."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


COLORS = {
    "ink": "#18212f",
    "blue": "#1261a0",
    "cyan": "#38a3a5",
    "orange": "#e07a2f",
    "grey": "#8993a4",
    "light": "#edf2f7",
}


def _save_figure(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_rankings(results: pd.DataFrame, path: str | Path) -> None:
    ordered = results.sort_values(["locus_order", "rank"], ascending=[True, False]).copy()
    labels = ordered["locus_label"] + "  ·  " + ordered["gene_symbol"]
    y = np.arange(len(ordered))
    means = ordered["posterior_mean"].to_numpy()
    lower = ordered["posterior_lower"].to_numpy()
    upper = ordered["posterior_upper"].to_numpy()
    colors = np.where(ordered["rank"].to_numpy() == 1, COLORS["blue"], COLORS["grey"])

    fig, ax = plt.subplots(figsize=(10, max(7, 0.36 * len(ordered))))
    ax.barh(y, means, color=colors, alpha=0.9)
    ax.errorbar(
        means,
        y,
        xerr=np.vstack([means - lower, upper - means]),
        fmt="none",
        ecolor=COLORS["ink"],
        elinewidth=0.8,
        capsize=2,
    )
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Posterior probability (including unresolved mass)")
    ax.set_title("BayesV2G LDL candidate-gene prioritization")
    ax.grid(axis="x", color=COLORS["light"], linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    _save_figure(fig, Path(path))


def plot_coefficients(coefficients: pd.DataFrame, path: str | Path) -> None:
    ordered = coefficients.sort_values("posterior_mean")
    y = np.arange(len(ordered))
    means = ordered["posterior_mean"].to_numpy()
    lower = ordered["posterior_lower"].to_numpy()
    upper = ordered["posterior_upper"].to_numpy()
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.axvline(0, color=COLORS["ink"], linewidth=0.8)
    ax.errorbar(
        means,
        y,
        xerr=np.vstack([means - lower, upper - means]),
        fmt="o",
        color=COLORS["blue"],
        ecolor=COLORS["cyan"],
        capsize=3,
    )
    ax.set_yticks(y, ordered["feature"])
    ax.set_xlabel("Log-odds coefficient (90% Laplace interval)")
    ax.set_title("Learned evidence weights")
    ax.grid(axis="x", color=COLORS["light"])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    _save_figure(fig, Path(path))


def plot_calibration(calibration: pd.DataFrame, path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(5.5, 5.2))
    ax.plot([0, 1], [0, 1], linestyle="--", color=COLORS["grey"], label="ideal")
    if not calibration.empty:
        sizes = 25 + 4 * calibration["n"].to_numpy()
        ax.scatter(
            calibration["mean_confidence"],
            calibration["observed_accuracy"],
            s=sizes,
            color=COLORS["blue"],
            alpha=0.85,
            label="held-out loci",
        )
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="Mean top-gene probability", ylabel="Observed top-1 accuracy")
    ax.set_title("Held-out calibration")
    ax.grid(color=COLORS["light"])
    ax.legend(frameon=False)
    _save_figure(fig, Path(path))


def _table(frame: pd.DataFrame, formats: dict[str, str] | None = None) -> str:
    display = frame.copy()
    for column, formatter in (formats or {}).items():
        if column in display:
            display[column] = display[column].map(lambda value: formatter.format(value) if pd.notna(value) else "NA")
    return display.to_html(index=False, escape=True, border=0, classes="data-table")


def write_html_report(
    output_path: str | Path,
    *,
    config: dict[str, Any],
    results: pd.DataFrame,
    loci: pd.DataFrame,
    coefficients: pd.DataFrame,
    metrics: pd.DataFrame,
    sensitivity: pd.DataFrame,
    manifest: dict[str, Any],
    sort1_case: dict[str, Any],
) -> None:
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    top = results.loc[results["rank"] == 1, [
        "locus_label",
        "lead_variant",
        "gene_symbol",
        "surfaced_candidate_count",
        "candidate_set_is_singleton",
        "conditional_posterior_mean",
        "posterior_mean",
        "posterior_lower",
        "posterior_upper",
        "unresolved_probability",
        "distance_to_tss_bp",
        "evidence_coverage",
    ]].copy()
    top["distance_to_tss_kb"] = top.pop("distance_to_tss_bp") / 1_000
    top = top.rename(
        columns={
            "locus_label": "Locus",
            "lead_variant": "Lead variant",
            "gene_symbol": "Top gene",
            "surfaced_candidate_count": "Candidates",
            "candidate_set_is_singleton": "Singleton",
            "conditional_posterior_mean": "Conditional posterior",
            "posterior_mean": "Absolute posterior",
            "posterior_lower": "P05",
            "posterior_upper": "P95",
            "unresolved_probability": "Unresolved",
            "distance_to_tss_kb": "Distance (kb)",
            "evidence_coverage": "Evidence coverage",
        }
    )
    metric_view = metrics.rename(columns={"model": "Model"})
    coefficient_view = coefficients.rename(
        columns={
            "feature": "Feature",
            "posterior_mean": "Mean",
            "posterior_sd": "SD",
            "posterior_lower": "P05",
            "posterior_upper": "P95",
            "odds_ratio": "exp(Mean)",
        }
    )
    unstable = sensitivity.loc[sensitivity.get("top_gene_changed", False) == True]  # noqa: E712
    singleton_loci = int(top["Singleton"].sum())
    provenance = {
        "GWAS study": config["study_id"],
        "Trait": config["trait"],
        "Genome build": config["genome_build"],
        "Open Targets release": config["open_targets_release"],
        "Training release": config["training_release"],
        "Training revision": config["training_revision"],
        "Snapshot SHA-256": manifest["input_files"]["open_targets_snapshot"]["sha256"],
        "Analysis source SHA-256": manifest["analysis_source"]["aggregate_sha256"],
        "Open Targets TLS verified": manifest["source_tls_verification"]["open_targets_snapshot"],
        "SORT1 JSON/metadata TLS verified": manifest["source_tls_verification"]["sort1_case_json_and_metadata"],
    }
    provenance_html = "".join(
        f"<dt>{html.escape(key)}</dt><dd>{html.escape(str(value))}</dd>" for key, value in provenance.items()
    )
    warning = (
        "This is a research prioritization result, not clinical evidence. Probabilities are conditional "
        "on an externally prefiltered candidate set, the one-primary-gene model, the listed data releases and the "
        f"{config['unresolved_mass']:.0%} outside-candidate assumption."
    )
    sort1_genes = sort1_case["genes"][
        [
            "geneSymbol",
            "distanceToLeadBp",
            "nes",
            "pValue",
            "rank",
            "posterior_mean",
            "l2g_score_external",
        ]
    ].rename(
        columns={
            "geneSymbol": "Gene",
            "distanceToLeadBp": "GTEx v39 TSS distance (bp)",
            "nes": "GTEx liver NES",
            "pValue": "GTEx p-value",
            "rank": "BayesV2G rank",
            "posterior_mean": "BayesV2G posterior",
            "l2g_score_external": "External L2G",
        }
    )
    sort1_variants = sort1_case["top_variants"].head(5)[
        ["rsid", "base_pair_location", "beta", "standard_error", "p_value", "approximate_pip"]
    ].rename(
        columns={
            "rsid": "Variant",
            "base_pair_location": "GRCh38 position",
            "beta": "Beta",
            "standard_error": "SE",
            "p_value": "GWAS p-value",
            "approximate_pip": "Approx. PIP",
        }
    )
    sort1_coloc = sort1_case["colocalization"][
        [
            "evidence_type",
            "dataset_id",
            "gene_symbol",
            "shared_snv_count",
            "pp_h3",
            "pp_h4",
            "pp_h4_prior_variance_0_01",
            "pp_h4_prior_variance_0_16",
            "top_shared_variant",
            "top_shared_pip_given_h4",
            "shared_credible_set_size_95",
        ]
    ].rename(
        columns={
            "evidence_type": "Evidence",
            "dataset_id": "Dataset",
            "gene_symbol": "Gene",
            "shared_snv_count": "Shared SNVs",
            "pp_h3": "PP(H3)",
            "pp_h4": "PP(H4)",
            "pp_h4_prior_variance_0_01": "PP(H4), W=.01",
            "pp_h4_prior_variance_0_16": "PP(H4), W=.16",
            "top_shared_variant": "Top shared variant",
            "top_shared_pip_given_h4": "Variant PIP | H4",
            "shared_credible_set_size_95": "Shared CS size",
        }
    )
    paired = manifest["held_out_paired_comparison"]
    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>BayesV2G LDL report</title>
  <style>
    :root {{ --ink:#18212f; --blue:#1261a0; --cyan:#38a3a5; --paper:#f7f9fc; --line:#dce3ec; }}
    body {{ margin:0; color:var(--ink); font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; background:var(--paper); }}
    main {{ max-width:1160px; margin:0 auto; padding:42px 28px 80px; }}
    h1 {{ font-size:2.5rem; line-height:1.08; margin:0 0 10px; letter-spacing:-.03em; }}
    h2 {{ margin-top:42px; border-bottom:1px solid var(--line); padding-bottom:8px; }}
    .lede {{ font-size:1.08rem; max-width:850px; }}
    .warning {{ background:#fff4dd; border-left:4px solid #e07a2f; padding:13px 16px; margin:24px 0; }}
    .card {{ background:white; border:1px solid var(--line); border-radius:12px; padding:20px; margin:18px 0; box-shadow:0 4px 18px rgba(24,33,47,.04); overflow:auto; }}
    dl {{ display:grid; grid-template-columns:210px 1fr; gap:5px 18px; margin:0; }} dt {{ font-weight:700; }} dd {{ margin:0; overflow-wrap:anywhere; }}
    .data-table {{ width:100%; border-collapse:collapse; font-size:.9rem; }}
    .data-table th {{ text-align:left; background:#edf2f7; }} .data-table th,.data-table td {{ padding:8px 9px; border-bottom:1px solid var(--line); white-space:nowrap; }}
    img {{ width:100%; height:auto; display:block; }} code {{ background:#edf2f7; padding:2px 5px; border-radius:4px; }}
    .muted {{ color:#5f6b7a; }}
  </style>
</head>
<body><main>
  <h1>BayesV2G</h1>
  <p class="lede">A Bayesian, locus-conditional ranking of candidate genes at ten real LDL-cholesterol GWAS signals.</p>
  <div class="warning">{html.escape(warning)}</div>
  <h2>Top result at each signal</h2>
  <p>The live API exposes candidates only after an external L2G score threshold of 0.05. The table therefore shows both the within-set conditional probability and the probability after unresolved mass. {singleton_loci} loci have only one surfaced candidate; their conditional value is structurally 1.0 and is not evidence-derived certainty.</p>
  <div class="card">{_table(top, {'Conditional posterior': '{:.3f}', 'Absolute posterior': '{:.3f}', 'P05': '{:.3f}', 'P95': '{:.3f}', 'Unresolved': '{:.3f}', 'Distance (kb)': '{:.1f}', 'Evidence coverage': '{:.0%}'})}</div>
  <div class="card"><img src="figures/posterior_rankings.png" alt="Posterior candidate-gene rankings"></div>
  <h2>Raw-data audit: rs12740374</h2>
  <p>A separately downloaded, Tabix-sliced 1 Mb region from GWAS Catalog study <code>GCST90018961</code> contains {sort1_case['gwas_variant_count']:,} real summary-statistic rows. A transparent single-signal Wakefield-ABF check assigns the lead variant an approximate PIP of {sort1_case['lead_variant_pip']:.3f}. This approximation ignores LD and is a teaching/audit calculation; the main ranking uses Open Targets' SuSiE-inf credible set from <code>GCST90019512</code>.</p>
  <div class="card">{_table(sort1_variants, {'Beta': '{:.4f}', 'SE': '{:.4f}', 'Approx. PIP': '{:.4f}'})}</div>
  <p>Regional, allele-matched summary statistics from eQTL Catalogue were then colocalized with the GWAS using a classic five-hypothesis single-signal ABF calculation. The liver eQTL data support a shared signal for all three genes, while the available plasma SORT1 pQTL does not show comparable H4 support. This calculation assumes at most one causal variant per trait and ignores LD. Because PP(H4) is prior-dependent, the table includes common effect-prior variances W=.01, .04 (the primary column), and .16.</p>
  <div class="card">{_table(sort1_coloc, {'PP(H3)': '{:.4f}', 'PP(H4)': '{:.4f}', 'PP(H4), W=.01': '{:.4f}', 'PP(H4), W=.16': '{:.4f}', 'Variant PIP | H4': '{:.4f}'})}</div>
  <p>At the same lead variant, GTEx v10 liver reports strong eQTL effects for all three highlighted genes. PSRC1 is nearest by GTEx v39 TSS, while integrated evidence changes the ranking—exactly the distinction between proximity and mediation that BayesV2G is designed to expose.</p>
  <div class="card">{_table(sort1_genes, {'GTEx liver NES': '{:.3f}', 'GTEx p-value': '{:.2e}', 'BayesV2G posterior': '{:.3f}', 'External L2G': '{:.3f}'})}</div>
  <h2>Model and held-out validation</h2>
  <p>The model is a conditional categorical regression: genes compete within each independently fine-mapped signal and conditional probabilities sum to one before explicit unresolved mass is added. Coefficients use Normal(0, {config['prior_scale']}²) priors and a Laplace posterior approximation.</p>
  <p>Against the distance-only comparator, the paired top-1 difference is {paired['top1_difference']:.3f} (locus-bootstrap 95% interval {paired['ci_lower_95']:.3f} to {paired['ci_upper_95']:.3f}; {int(paired['bootstrap_iterations']):,} resamples).</p>
  <div class="card">{_table(metric_view, {'top1_accuracy': '{:.3f}', 'recall_at_3': '{:.3f}', 'mean_reciprocal_rank': '{:.3f}', 'multiclass_log_loss': '{:.3f}', 'brier_score': '{:.3f}', 'mean_entropy': '{:.3f}'})}</div>
  <div class="card"><img src="figures/calibration.png" alt="Held-out probability calibration"></div>
  <h2>Learned evidence weights</h2>
  <p>Coefficients are predictive weights for a full one-unit contrast in bounded, correlated evidence features. They are not causal effects and their magnitudes are not directly comparable without considering each feature's distribution.</p>
  <div class="card">{_table(coefficient_view, {'Mean': '{:.3f}', 'SD': '{:.3f}', 'P05': '{:.3f}', 'P95': '{:.3f}', 'exp(Mean)': '{:.2f}'})}</div>
  <div class="card"><img src="figures/coefficient_posteriors.png" alt="Posterior coefficient intervals"></div>
  <h2>Sensitivity</h2>
  <p>Prior-scale and leave-one-evidence-class-out scenarios were refitted at the MAP estimate; unresolved-mass assumptions of 10% and 40% were also evaluated. {len(unstable)} locus-scenario combinations changed the top-ranked gene.</p>
  <div class="card">{_table(sensitivity, {'top_probability': '{:.3f}', 'reference_probability': '{:.3f}', 'max_probability_change': '{:.3f}'})}</div>
  <h2>Data provenance</h2>
  <div class="card"><dl>{provenance_html}</dl></div>
  <p>Raw credible-set variants, PIPs, candidate features, model artifacts and the complete run manifest are beside this report. The Open Targets external L2G score is retained only as a comparison column and is not used by the Bayesian model.</p>
  <h2>Interpretation limits</h2>
  <ul>
    <li>The API surfaces Open Targets candidates whose external L2G score is at least 0.05; omitted genes are represented only through unresolved mass.</li>
    <li>Training labels are curated biological knowledge and may partly rely on the same evidence classes, so held-out accuracy can still be optimistic.</li>
    <li>QTL feature values summarize the maximum across available studies/tissues and zeros do not cleanly distinguish a true null from unavailable or underpowered evidence.</li>
    <li>Out-of-sample LD was used for these SuSiE-inf credible sets. Ancestry/LD mismatch and multiple biological mediators remain model limitations.</li>
  </ul>
  <p class="muted">Generated {html.escape(manifest['generated_at_utc'])}. Reproduce with <code>bayesv2g run</code>.</p>
</main></body></html>"""
    destination.write_text(content, encoding="utf-8")


def write_json(payload: dict[str, Any], path: str | Path) -> None:
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
