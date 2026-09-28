import numpy as np
import pandas as pd

from bayesv2g.evaluation import (
    baseline_predictions,
    evaluate_ranking,
    grouped_bootstrap_interval,
    paired_top1_bootstrap,
)


def test_grouped_ranking_metrics_and_baselines():
    frame = pd.DataFrame(
        {
            "studyLocusId": ["a", "a", "b", "b", "b"],
            "goldStandardSet": [1, 0, 0, 1, 0],
            "posterior_mean": [0.8, 0.2, 0.1, 0.7, 0.2],
            "distanceTssMean": [0.9, 0.2, 0.1, 0.8, 0.3],
            "eQtlColocH4Maximum": [0.9, 0, 0, 0.8, 0],
            "vepMean": [0.2, 0.1, 0.1, 0.8, 0.1],
        }
    )
    metrics = evaluate_ranking(frame)
    assert metrics["top1_accuracy"] == 1.0
    assert metrics["recall_at_3"] == 1.0
    baselines = baseline_predictions(frame)
    assert {"uniform", "distance_only", "eqtl_coloc_only", "vep_only"} == set(baselines)
    for probabilities in baselines.values():
        sums = pd.Series(probabilities).groupby(frame["studyLocusId"]).sum()
        np.testing.assert_allclose(sums, 1.0)


def test_grouped_bootstrap_keeps_duplicate_draws_as_distinct_loci():
    frame = pd.DataFrame(
        {
            "studyLocusId": ["a", "b"],
            "value": [1.0, 2.0],
        }
    )
    interval = grouped_bootstrap_interval(
        frame,
        lambda sample: float(sample["studyLocusId"].nunique()),
        iterations=40,
        seed=1,
    )
    assert interval == (2.0, 2.0)


def test_paired_top1_bootstrap_uses_locus_level_differences():
    frame = pd.DataFrame(
        {
            "studyLocusId": ["a", "a", "b", "b"],
            "goldStandardSet": [1, 0, 0, 1],
            "model": [0.8, 0.2, 0.2, 0.8],
            "baseline": [0.2, 0.8, 0.8, 0.2],
        }
    )
    result = paired_top1_bootstrap(
        frame,
        "model",
        "baseline",
        iterations=50,
        seed=3,
    )
    assert result["top1_difference"] == 1.0
    assert result["ci_lower_95"] == 1.0
    assert result["improved_loci"] == 2.0
    assert result["worsened_loci"] == 0.0
