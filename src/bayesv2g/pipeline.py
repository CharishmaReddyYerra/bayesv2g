"""End-to-end orchestration for the real-data LDL demonstration."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy

from bayesv2g import __version__
from bayesv2g.case_study import analyse_sort1_case
from bayesv2g.config import load_config, project_root
from bayesv2g.data_sources import (
    fetch_open_targets_loci,
    sha256_file,
    validate_open_targets_snapshot,
    write_json,
)
from bayesv2g.evaluation import (
    baseline_predictions,
    calibration_table,
    evaluate_baselines,
    evaluate_ranking,
    paired_top1_bootstrap,
)
from bayesv2g.features import evidence_coverage, flatten_open_targets
from bayesv2g.model import BayesianSoftmaxModel
from bayesv2g.reporting import (
    plot_calibration,
    plot_coefficients,
    plot_rankings,
    write_html_report,
)


def _coefficient_table(model: BayesianSoftmaxModel) -> pd.DataFrame:
    standard_error = np.sqrt(np.diag(model.covariance_))
    return pd.DataFrame(
        {
            "feature": model.feature_names,
            "posterior_mean": model.coef_,
            "posterior_sd": standard_error,
            "posterior_lower": model.coef_ - 1.6448536269514722 * standard_error,
            "posterior_upper": model.coef_ + 1.6448536269514722 * standard_error,
            "odds_ratio": np.exp(model.coef_),
        }
    )


def _add_absolute_probabilities(
    predictions: pd.DataFrame,
    loci: pd.DataFrame,
    unresolved_mass: float,
) -> pd.DataFrame:
    if not 0 <= unresolved_mass < 1:
        raise ValueError("unresolved_mass must lie in [0, 1)")
    coverage = loci.set_index("studyLocusId")["pip_coverage_95"].clip(0, 1)
    candidate_counts = loci.set_index("studyLocusId")["candidate_count"]
    result = predictions.copy()
    result = result.rename(
        columns={
            "posterior_mean": "conditional_posterior_mean",
            "posterior_lower": "conditional_posterior_lower",
            "posterior_upper": "conditional_posterior_upper",
        }
    )
    result["pip_coverage_95"] = result["studyLocusId"].map(coverage)
    resolvable = (1.0 - unresolved_mass) * result["pip_coverage_95"]
    result["posterior_mean"] = result["conditional_posterior_mean"] * resolvable
    result["posterior_lower"] = result["conditional_posterior_lower"] * resolvable
    result["posterior_upper"] = result["conditional_posterior_upper"] * resolvable
    result["unresolved_probability"] = 1.0 - resolvable
    result["surfaced_candidate_count"] = result["studyLocusId"].map(candidate_counts).astype(int)
    result["candidate_set_is_singleton"] = result["surfaced_candidate_count"].eq(1)
    result["probability_scope"] = "conditional on Open Targets score >= 0.05 candidate ascertainment"
    result["rank"] = (
        result.groupby("studyLocusId")["posterior_mean"].rank(method="first", ascending=False).astype(int)
    )
    result["evidence_coverage"] = evidence_coverage(result)
    entropy = result.groupby("studyLocusId")["conditional_posterior_mean"].transform(
        lambda values: float(-(values * np.log(np.clip(values, 1e-15, 1))).sum())
    )
    result["locus_entropy"] = entropy
    return result


def _configured_path(
    root: Path,
    config: dict[str, Any],
    key: str,
    default: str,
) -> Path:
    configured = Path(config.get(key, default)).expanduser()
    return configured if configured.is_absolute() else root / configured


def _portable_manifest_path(path: Path, root: Path) -> str:
    """Use repository-relative paths for inputs stored within the project."""
    resolved_path = path.resolve()
    resolved_root = root.resolve()
    try:
        return resolved_path.relative_to(resolved_root).as_posix()
    except ValueError:
        return resolved_path.as_posix()


def _analysis_source_provenance(root: Path) -> dict[str, Any]:
    candidates = [
        root / "pyproject.toml",
        root / "config" / "demo_ldl.json",
        root / "requirements.lock",
        *sorted((root / "src" / "bayesv2g").glob("*.py")),
        root / "src" / "bayesv2g" / "default_config.json",
    ]
    existing = [path for path in candidates if path.exists()]
    hashes = {str(path.relative_to(root)): sha256_file(path) for path in existing}
    aggregate = hashlib.sha256()
    for relative_path, digest in sorted(hashes.items()):
        aggregate.update(relative_path.encode("utf-8"))
        aggregate.update(b"\0")
        aggregate.update(digest.encode("ascii"))
        aggregate.update(b"\n")
    return {"aggregate_sha256": aggregate.hexdigest(), "files": hashes}


def _sensitivity_analysis(
    train: pd.DataFrame,
    candidates: pd.DataFrame,
    loci: pd.DataFrame,
    config: dict[str, Any],
    reference: pd.DataFrame,
) -> pd.DataFrame:
    reference_top = (
        reference.loc[reference["rank"] == 1]
        .set_index("studyLocusId")[["gene_symbol", "posterior_mean"]]
        .rename(columns={"gene_symbol": "reference_gene", "posterior_mean": "reference_probability"})
    )
    scenarios: list[tuple[str, list[str], float]] = []
    for scale in (0.5, 1.5):
        scenarios.append((f"prior_scale_{scale}", list(config["features"]), scale))
    for feature in config["features"]:
        scenarios.append((f"without_{feature}", [item for item in config["features"] if item != feature], config["prior_scale"]))

    records: list[dict[str, Any]] = []

    def record_scenario(name: str, absolute: pd.DataFrame, scenario_type: str) -> None:
        top = absolute.loc[absolute["rank"] == 1].set_index("studyLocusId")
        for locus_id, row in top.iterrows():
            reference_row = reference_top.loc[locus_id]
            merged = absolute.loc[absolute["studyLocusId"] == locus_id].merge(
                reference.loc[
                    reference["studyLocusId"] == locus_id,
                    ["geneId", "posterior_mean"],
                ],
                on="geneId",
                suffixes=("_scenario", "_reference"),
            )
            records.append(
                {
                    "scenario": name,
                    "scenario_type": scenario_type,
                    "studyLocusId": locus_id,
                    "locus_label": row["locus_label"],
                    "top_gene": row["gene_symbol"],
                    "top_probability": row["posterior_mean"],
                    "reference_gene": reference_row["reference_gene"],
                    "reference_probability": reference_row["reference_probability"],
                    "top_gene_changed": row["gene_symbol"] != reference_row["reference_gene"],
                    "max_probability_change": float(
                        (
                            merged["posterior_mean_scenario"]
                            - merged["posterior_mean_reference"]
                        )
                        .abs()
                        .max()
                    ),
                }
            )

    for name, features, scale in scenarios:
        model = BayesianSoftmaxModel(
            feature_names=features,
            prior_scale=scale,
            n_posterior_draws=min(600, config["posterior_draws"]),
            random_seed=config["random_seed"],
        ).fit(train)
        predicted = candidates.copy()
        map_probability = model.predict_map(candidates)
        predicted["posterior_mean"] = map_probability
        predicted["posterior_lower"] = map_probability
        predicted["posterior_upper"] = map_probability
        absolute = _add_absolute_probabilities(predicted, loci, config["unresolved_mass"])
        record_scenario(name, absolute, "model_refit_map")

    conditional_reference = reference.drop(
        columns=[
            "posterior_mean",
            "posterior_lower",
            "posterior_upper",
            "conditional_posterior_mean",
            "conditional_posterior_lower",
            "conditional_posterior_upper",
        ],
        errors="ignore",
    )
    conditional_reference["posterior_mean"] = reference[
        "conditional_posterior_mean"
    ].to_numpy()
    conditional_reference["posterior_lower"] = reference[
        "conditional_posterior_lower"
    ].to_numpy()
    conditional_reference["posterior_upper"] = reference[
        "conditional_posterior_upper"
    ].to_numpy()
    for rho in (0.10, 0.40):
        predicted = conditional_reference.drop(
            columns=[
                "pip_coverage_95",
                "unresolved_probability",
                "surfaced_candidate_count",
                "candidate_set_is_singleton",
                "probability_scope",
                "rank",
                "evidence_coverage",
                "locus_entropy",
            ],
            errors="ignore",
        )
        absolute = _add_absolute_probabilities(predicted, loci, rho)
        record_scenario(f"unresolved_mass_{rho}", absolute, "probability_assumption")
    return pd.DataFrame(records)


def run_pipeline(
    *,
    config_path: str | Path | None = None,
    output_directory: str | Path | None = None,
    refresh_data: bool = False,
    verify_tls: bool = True,
    posterior_draws: int | None = None,
) -> dict[str, Path]:
    root = project_root()
    config = load_config(config_path)
    if posterior_draws is not None:
        config["posterior_draws"] = int(posterior_draws)
    output = Path(output_directory) if output_directory else root / "results" / "ldl_demo"
    output.mkdir(parents=True, exist_ok=True)
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    snapshot_path = _configured_path(
        root,
        config,
        "open_targets_snapshot",
        "data/snapshots/open_targets_26.09_ldl_loci.json",
    )
    ids = [item["study_locus_id"] for item in config["demo_loci"]]
    if refresh_data:
        snapshot = fetch_open_targets_loci(
            config["open_targets_endpoint"], ids, verify_tls=verify_tls
        )
        validate_open_targets_snapshot(snapshot, ids, config["open_targets_release"])
        write_json(snapshot, snapshot_path)
    elif not snapshot_path.exists():
        raise FileNotFoundError(
            f"Open Targets snapshot is missing: {snapshot_path}. "
            "Run `bayesv2g fetch` or rerun with `bayesv2g run --refresh`."
        )
    else:
        snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
    validate_open_targets_snapshot(snapshot, ids, config["open_targets_release"])

    labels = {item["study_locus_id"]: item["label"] for item in config["demo_loci"]}
    candidates, loci, variants = flatten_open_targets(snapshot, labels)
    locus_order = {item["study_locus_id"]: index for index, item in enumerate(config["demo_loci"])}
    candidates["locus_order"] = candidates["studyLocusId"].map(locus_order)

    training_dir = _configured_path(
        root,
        config,
        "training_directory",
        "data/raw/open_targets_l2g_26.06",
    )
    train_path = training_dir / "train.parquet"
    test_path = training_dir / "test.parquet"
    if not train_path.exists() or not test_path.exists():
        raise FileNotFoundError(
            "Versioned Open Targets training artifacts are missing. See data/README.md for the pinned hf download command."
        )
    train = pd.read_parquet(train_path)
    test = pd.read_parquet(test_path)
    overlap_mask = train["studyLocusId"].isin(ids)
    excluded_rows = int(overlap_mask.sum())
    excluded_loci = int(train.loc[overlap_mask, "studyLocusId"].nunique())
    train = train.loc[~overlap_mask].reset_index(drop=True)

    model = BayesianSoftmaxModel(
        feature_names=config["features"],
        prior_scale=config["prior_scale"],
        n_posterior_draws=config["posterior_draws"],
        random_seed=config["random_seed"],
    ).fit(train)
    model_path = output / "model.json"
    model.save(model_path)

    predictions = model.predict(candidates)
    results = _add_absolute_probabilities(predictions, loci, config["unresolved_mass"])
    results = results.sort_values(["locus_order", "rank"]).reset_index(drop=True)

    reference_map = candidates.copy()
    reference_map_probability = model.predict_map(candidates)
    reference_map["posterior_mean"] = reference_map_probability
    reference_map["posterior_lower"] = reference_map_probability
    reference_map["posterior_upper"] = reference_map_probability
    reference_map = _add_absolute_probabilities(
        reference_map,
        loci,
        config["unresolved_mass"],
    )

    test_scored = test.copy()
    test_scored["posterior_mean"] = model.predict_map(test)
    model_metrics = {"model": "BayesV2G", **evaluate_ranking(test_scored)}
    baseline_metrics = evaluate_baselines(test)
    metrics = pd.concat([pd.DataFrame([model_metrics]), baseline_metrics], ignore_index=True)
    distance_probability = baseline_predictions(test)["distance_only"]
    paired_frame = test_scored.copy()
    paired_frame["distance_probability"] = distance_probability
    paired_comparison = {
        "first_model": "BayesV2G",
        "second_model": "distance_only",
        **paired_top1_bootstrap(
            paired_frame,
            "posterior_mean",
            "distance_probability",
            iterations=10_000,
            seed=config["random_seed"],
        ),
    }
    calibration = calibration_table(test_scored)
    coefficients = _coefficient_table(model)
    sensitivity = _sensitivity_analysis(train, candidates, loci, config, reference_map)

    candidates.to_csv(output / "candidate_evidence.csv", index=False)
    loci.to_csv(output / "loci.csv", index=False)
    variants.to_csv(output / "credible_set_variants.csv", index=False)
    results.to_csv(output / "prioritized_genes.csv", index=False)
    coefficients.to_csv(output / "posterior_coefficients.csv", index=False)
    metrics.to_csv(output / "validation_metrics.csv", index=False)
    pd.DataFrame([paired_comparison]).to_csv(
        output / "paired_top1_comparison.csv",
        index=False,
    )
    calibration.to_csv(output / "calibration.csv", index=False)
    sensitivity.to_csv(output / "sensitivity.csv", index=False)
    sort1_source = _configured_path(
        root,
        config,
        "sort1_case_directory",
        "data/snapshots/sort1_case",
    )
    sort1_case = analyse_sort1_case(sort1_source, output, results)

    generated_at = datetime.now(timezone.utc).isoformat()
    release = snapshot.get("data", {}).get("meta", {}).get("dataVersion", {})
    manifest = {
        "generated_at_utc": generated_at,
        "bayesv2g_version": __version__,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "analysis_source": _analysis_source_provenance(root),
        "config": config,
        "open_targets_api_release": release,
        "candidate_definition": "Open Targets Platform L2G predictions surfaced at external score >= 0.05",
        "external_l2g_score_used_as_feature": False,
        "probability_scope": "conditional on Open Targets score >= 0.05 candidate ascertainment",
        "demo_loci_excluded_from_training": {"rows": excluded_rows, "loci": excluded_loci},
        "training_fit": model.fit_metadata_,
        "held_out_paired_comparison": paired_comparison,
        "input_files": {
            "open_targets_snapshot": {
                "path": _portable_manifest_path(snapshot_path, root),
                "sha256": sha256_file(snapshot_path),
            },
            "training": {
                "path": _portable_manifest_path(train_path, root),
                "sha256": sha256_file(train_path),
            },
            "test": {
                "path": _portable_manifest_path(test_path, root),
                "sha256": sha256_file(test_path),
            },
            "sort1_case_manifest": {
                "path": _portable_manifest_path(sort1_source / "manifest.json", root),
                "sha256": sha256_file(sort1_source / "manifest.json"),
            },
        },
        "source_tls_verification": {
            "open_targets_snapshot": bool(
                snapshot.get("bayesv2g_snapshot", {}).get("tls_verified", False)
            ),
            "sort1_case_json_and_metadata": bool(sort1_case["source_tls_verified"]),
        },
        "sort1_case_study": {
            key: value
            for key, value in sort1_case.items()
            if key not in {
                "top_variants",
                "genes",
                "colocalization",
                "source_manifest",
            }
        },
    }
    manifest_path = output / "run_manifest.json"
    write_json(manifest, manifest_path)

    plot_rankings(results, figures / "posterior_rankings.png")
    plot_coefficients(coefficients, figures / "coefficient_posteriors.png")
    plot_calibration(calibration, figures / "calibration.png")
    report_path = output / "report.html"
    write_html_report(
        report_path,
        config=config,
        results=results,
        loci=loci,
        coefficients=coefficients,
        metrics=metrics,
        sensitivity=sensitivity,
        manifest=manifest,
        sort1_case=sort1_case,
    )
    checksum_path = output / "SHA256SUMS"
    artifact_files = sorted(
        path for path in output.rglob("*") if path.is_file() and path != checksum_path
    )
    checksum_path.write_text(
        "".join(
            f"{sha256_file(path)}  {path.relative_to(output)}\n"
            for path in artifact_files
        ),
        encoding="utf-8",
    )
    return {
        "output_directory": output,
        "report": report_path,
        "prioritized_genes": output / "prioritized_genes.csv",
        "manifest": manifest_path,
        "model": model_path,
        "checksums": checksum_path,
    }
