"""Grouped ranking metrics and transparent comparison baselines."""

from __future__ import annotations

from collections.abc import Callable
from math import comb

import numpy as np
import pandas as pd
from scipy.special import softmax


def _tie_aware_rank_metrics(
    probabilities: np.ndarray,
    positive: np.ndarray,
    cutoff: int = 3,
) -> tuple[float, float, float]:
    """Expected top-1, recall@k and reciprocal rank under random tie order."""

    order = np.argsort(-probabilities, kind="stable")
    sorted_probabilities = probabilities[order]
    sorted_positive = positive[order]
    blocks: list[tuple[int, int]] = []
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and np.isclose(
            sorted_probabilities[stop], sorted_probabilities[start], rtol=1e-12, atol=1e-15
        ):
            stop += 1
        blocks.append((start, stop))
        start = stop

    first_start, first_stop = blocks[0]
    first_positive = int(sorted_positive[first_start:first_stop].sum())
    top1 = first_positive / (first_stop - first_start)

    recall = 0.0
    reciprocal_rank = 0.0
    for block_start, block_stop in blocks:
        block = sorted_positive[block_start:block_stop]
        block_size = len(block)
        positive_count = int(block.sum())
        if positive_count == 0:
            continue
        slots = max(0, min(cutoff - block_start, block_size))
        if block_start < cutoff and slots:
            recall = 1.0 - (
                comb(block_size - positive_count, slots) / comb(block_size, slots)
                if block_size - positive_count >= slots
                else 0.0
            )
        denominator = comb(block_size, positive_count)
        reciprocal_rank = sum(
            (comb(block_size - local_rank, positive_count - 1) / denominator)
            / (block_start + local_rank)
            for local_rank in range(1, block_size - positive_count + 2)
        )
        break
    return float(top1), float(recall), float(reciprocal_rank)


def evaluate_ranking(
    frame: pd.DataFrame,
    probability_column: str = "posterior_mean",
    *,
    group_column: str = "studyLocusId",
    label_column: str = "goldStandardSet",
) -> dict[str, float]:
    """Evaluate one probability distribution per locus.

    If a locus has multiple accepted genes, their target mass is divided
    equally for the Brier score and any accepted gene counts as a hit.
    """

    required = {probability_column, group_column, label_column}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing evaluation columns: {sorted(missing)}")

    top1: list[float] = []
    recall3: list[float] = []
    reciprocal_ranks: list[float] = []
    log_losses: list[float] = []
    brier_scores: list[float] = []
    entropies: list[float] = []

    for _, group in frame.groupby(group_column, sort=False):
        probabilities = group[probability_column].to_numpy(dtype=float)
        if np.any(probabilities < 0) or not np.all(np.isfinite(probabilities)):
            raise ValueError("Probabilities must be finite and non-negative")
        total = probabilities.sum()
        if total <= 0:
            raise ValueError("Each locus needs positive probability mass")
        probabilities = probabilities / total
        positive = group[label_column].to_numpy(dtype=float) > 0
        if not positive.any():
            raise ValueError("Every evaluation locus must have a positive label")
        expected_top1, expected_recall3, expected_rr = _tie_aware_rank_metrics(
            probabilities, positive, cutoff=3
        )
        top1.append(expected_top1)
        recall3.append(expected_recall3)
        reciprocal_ranks.append(expected_rr)
        accepted_mass = probabilities[positive].sum()
        log_losses.append(float(-np.log(np.clip(accepted_mass, 1e-15, 1.0))))
        target = positive.astype(float) / positive.sum()
        brier_scores.append(float(np.square(probabilities - target).sum()))
        entropies.append(float(-(probabilities * np.log(np.clip(probabilities, 1e-15, 1))).sum()))

    return {
        "n_loci": float(len(top1)),
        "top1_accuracy": float(np.mean(top1)),
        "recall_at_3": float(np.mean(recall3)),
        "mean_reciprocal_rank": float(np.mean(reciprocal_ranks)),
        "multiclass_log_loss": float(np.mean(log_losses)),
        "brier_score": float(np.mean(brier_scores)),
        "mean_entropy": float(np.mean(entropies)),
    }


def probabilities_from_score(
    frame: pd.DataFrame,
    score_column: str,
    *,
    group_column: str = "studyLocusId",
    temperature: float = 12.0,
) -> np.ndarray:
    """Convert a higher-is-better score to within-locus probabilities."""

    if temperature <= 0:
        raise ValueError("temperature must be positive")
    result = np.zeros(len(frame), dtype=float)
    for _, indices in frame.groupby(group_column, sort=False).groups.items():
        positions = frame.index.get_indexer(indices)
        values = frame.loc[indices, score_column].fillna(0).to_numpy(dtype=float)
        if np.allclose(values, values[0]):
            result[positions] = 1.0 / len(values)
        else:
            result[positions] = softmax(temperature * values)
    return result


def baseline_predictions(frame: pd.DataFrame) -> dict[str, np.ndarray]:
    """Return probability vectors for prespecified, non-trained baselines."""

    uniform = np.zeros(len(frame), dtype=float)
    for _, indices in frame.groupby("studyLocusId", sort=False).groups.items():
        positions = frame.index.get_indexer(indices)
        uniform[positions] = 1.0 / len(indices)

    baselines: dict[str, np.ndarray] = {"uniform": uniform}
    columns = {
        "distance_only": "distanceTssMean",
        "eqtl_coloc_only": "eQtlColocH4Maximum",
        "vep_only": "vepMean",
    }
    for name, column in columns.items():
        if column in frame:
            baselines[name] = probabilities_from_score(frame, column)
    return baselines


def evaluate_baselines(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    for name, probabilities in baseline_predictions(frame).items():
        scored = frame.copy()
        scored["baseline_probability"] = probabilities
        rows.append({"model": name, **evaluate_ranking(scored, "baseline_probability")})
    return pd.DataFrame(rows)


def grouped_bootstrap_interval(
    frame: pd.DataFrame,
    metric: Callable[[pd.DataFrame], float],
    *,
    group_column: str = "studyLocusId",
    iterations: int = 500,
    seed: int = 20260925,
) -> tuple[float, float]:
    """Percentile interval from a locus-level bootstrap."""

    groups = [group for _, group in frame.groupby(group_column, sort=False)]
    rng = np.random.default_rng(seed)
    values = np.empty(iterations, dtype=float)
    for iteration in range(iterations):
        chosen = rng.integers(0, len(groups), size=len(groups))
        sample = pd.concat(
            [
                groups[index].assign(**{group_column: f"bootstrap_{offset}"})
                for offset, index in enumerate(chosen)
            ],
            ignore_index=True,
        )
        values[iteration] = metric(sample)
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def paired_top1_bootstrap(
    frame: pd.DataFrame,
    first_probability_column: str,
    second_probability_column: str,
    *,
    group_column: str = "studyLocusId",
    label_column: str = "goldStandardSet",
    iterations: int = 10_000,
    seed: int = 20260925,
) -> dict[str, float]:
    """Paired locus bootstrap for a difference in expected top-1 accuracy."""

    if not isinstance(iterations, int) or isinstance(iterations, bool) or iterations <= 0:
        raise ValueError("iterations must be a positive integer")
    required = {
        group_column,
        label_column,
        first_probability_column,
        second_probability_column,
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing paired-comparison columns: {sorted(missing)}")

    first: list[float] = []
    second: list[float] = []
    for _, group in frame.groupby(group_column, sort=False):
        positive = group[label_column].to_numpy(dtype=float) > 0
        if not positive.any():
            raise ValueError("Every comparison locus must have a positive label")
        group_values: list[float] = []
        for column in (first_probability_column, second_probability_column):
            probabilities = group[column].to_numpy(dtype=float)
            if np.any(probabilities < 0) or not np.all(np.isfinite(probabilities)):
                raise ValueError("Comparison probabilities must be finite and non-negative")
            if probabilities.sum() <= 0:
                raise ValueError("Each comparison locus needs positive probability mass")
            group_values.append(_tie_aware_rank_metrics(probabilities, positive, cutoff=1)[0])
        first.append(group_values[0])
        second.append(group_values[1])

    differences = np.asarray(first) - np.asarray(second)
    rng = np.random.default_rng(seed)
    bootstrap = np.empty(iterations, dtype=float)
    for index in range(iterations):
        bootstrap[index] = differences[
            rng.integers(0, len(differences), size=len(differences))
        ].mean()
    tolerance = 1e-12
    return {
        "n_loci": float(len(differences)),
        "top1_difference": float(differences.mean()),
        "ci_lower_95": float(np.quantile(bootstrap, 0.025)),
        "ci_upper_95": float(np.quantile(bootstrap, 0.975)),
        "improved_loci": float(np.sum(differences > tolerance)),
        "worsened_loci": float(np.sum(differences < -tolerance)),
        "tied_loci": float(np.sum(np.abs(differences) <= tolerance)),
        "bootstrap_iterations": float(iterations),
    }


def calibration_table(
    frame: pd.DataFrame,
    probability_column: str = "posterior_mean",
    *,
    bins: int = 10,
) -> pd.DataFrame:
    """Top-gene confidence calibration by locus."""

    records: list[dict[str, float]] = []
    for _, group in frame.groupby("studyLocusId", sort=False):
        top = group.loc[group[probability_column].idxmax()]
        records.append(
            {
                "confidence": float(top[probability_column]),
                "correct": float(top["goldStandardSet"] > 0),
            }
        )
    top_frame = pd.DataFrame(records)
    edges = np.linspace(0, 1, bins + 1)
    top_frame["bin"] = pd.cut(top_frame["confidence"], edges, include_lowest=True, duplicates="drop")
    result = (
        top_frame.groupby("bin", observed=True)
        .agg(mean_confidence=("confidence", "mean"), observed_accuracy=("correct", "mean"), n=("correct", "size"))
        .reset_index()
    )
    result["bin"] = result["bin"].astype(str)
    return result
