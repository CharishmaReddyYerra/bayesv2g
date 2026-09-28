from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from bayesv2g.model import BayesianSoftmaxModel, DEFAULT_FEATURES


FEATURES = ("distanceTssMean", "eQtlColocH4Maximum")


def toy_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "studyLocusId": ["a", "a", "a", "b", "b", "c", "c", "c"],
            "geneId": ["a1", "a2", "a3", "b1", "b2", "c1", "c2", "c3"],
            "goldStandardSet": [1, 0, 0, 1, 0, 1, 1, 0],
            "distanceTssMean": [0.95, 0.25, 0.10, 0.85, 0.15, 0.90, 0.80, 0.10],
            "eQtlColocH4Maximum": [0.90, 0.10, 0.20, 0.80, 0.05, 0.85, 0.75, 0.05],
        }
    )


def test_defaults_match_requested_evidence_features() -> None:
    assert DEFAULT_FEATURES == (
        "distanceTssMean",
        "eQtlColocH4Maximum",
        "pQtlColocH4Maximum",
        "sQtlColocH4Maximum",
        "vepMean",
        "e2gMean",
    )
    assert BayesianSoftmaxModel().prior_scale == pytest.approx(0.75)


def test_multiple_positives_share_one_locus_target_mass() -> None:
    frame = toy_frame()
    model = BayesianSoftmaxModel(feature_names=FEATURES)
    _, codes, n_groups = model._prepare_frame(frame, require_labels=True)
    target = model._make_targets(frame["goldStandardSet"], codes, n_groups)

    np.testing.assert_allclose(
        np.bincount(codes, weights=target, minlength=n_groups), 1.0
    )
    np.testing.assert_allclose(target[frame["studyLocusId"] == "c"], [0.5, 0.5, 0.0])


def test_fit_hessian_and_groupwise_predictions() -> None:
    frame = toy_frame()
    model = BayesianSoftmaxModel(
        feature_names=FEATURES, n_posterior_draws=500, random_seed=17
    ).fit(frame)

    assert model.is_fitted
    assert model.coef_.shape == (2,)
    assert np.all(model.coef_ > 0.0)
    np.testing.assert_allclose(model.hessian_, model.hessian_.T, atol=1e-12)
    np.testing.assert_allclose(model.covariance_, model.covariance_.T, atol=1e-12)
    np.testing.assert_allclose(
        model.hessian_ @ model.covariance_, np.eye(2), atol=1e-9
    )
    assert np.all(np.linalg.eigvalsh(model.hessian_) > 0.0)

    prediction = model.predict(frame)
    sums = prediction.groupby("studyLocusId", sort=False)["posterior_mean"].sum()
    np.testing.assert_allclose(sums.to_numpy(), 1.0, atol=1e-12)
    assert np.all(prediction["posterior_lower"] <= prediction["posterior_mean"])
    assert np.all(prediction["posterior_mean"] <= prediction["posterior_upper"])
    assert prediction.loc[0, "posterior_mean"] > prediction.loc[1, "posterior_mean"]
    assert prediction.loc[3, "posterior_mean"] > prediction.loc[4, "posterior_mean"]


def test_analytic_hessian_matches_finite_difference_gradient() -> None:
    frame = toy_frame()
    model = BayesianSoftmaxModel(feature_names=FEATURES)
    x, codes, n_groups = model._prepare_frame(frame, require_labels=True)
    target = model._make_targets(frame["goldStandardSet"], codes, n_groups)
    beta = np.array([0.31, -0.27])
    analytic = model._analytic_hessian(beta, x, codes, n_groups)

    epsilon = 1e-5
    finite_difference = np.empty_like(analytic)
    for column in range(beta.size):
        step = np.zeros_like(beta)
        step[column] = epsilon
        gradient_plus = model._objective_and_gradient(
            beta + step, x, codes, n_groups, target
        )[1]
        gradient_minus = model._objective_and_gradient(
            beta - step, x, codes, n_groups, target
        )[1]
        finite_difference[:, column] = (
            gradient_plus - gradient_minus
        ) / (2.0 * epsilon)

    np.testing.assert_allclose(analytic, finite_difference, rtol=1e-6, atol=1e-7)


def test_prediction_is_deterministic_and_save_load_round_trips(tmp_path) -> None:
    frame = toy_frame()
    model = BayesianSoftmaxModel(
        feature_names=FEATURES, n_posterior_draws=300, random_seed=123
    ).fit(frame)
    first = model.predict(frame)
    second = model.predict(frame)
    np.testing.assert_array_equal(
        first[["posterior_mean", "posterior_lower", "posterior_upper"]],
        second[["posterior_mean", "posterior_lower", "posterior_upper"]],
    )

    destination = tmp_path / "model.json"
    model.save(destination)
    payload = json.loads(destination.read_text(encoding="utf-8"))
    assert payload["model_type"] == "conditional_bayesian_softmax"

    restored = BayesianSoftmaxModel.load(destination)
    restored_prediction = restored.predict(frame)
    np.testing.assert_array_equal(
        first[["posterior_mean", "posterior_lower", "posterior_upper"]],
        restored_prediction[
            ["posterior_mean", "posterior_lower", "posterior_upper"]
        ],
    )
    np.testing.assert_array_equal(restored.coef_, model.coef_)
    np.testing.assert_array_equal(restored.covariance_, model.covariance_)


@pytest.mark.parametrize(
    ("mutation", "error", "message"),
    [
        (
            lambda data: data.drop(columns=["distanceTssMean"]),
            ValueError,
            "missing required columns",
        ),
        (
            lambda data: data.assign(distanceTssMean=np.nan),
            ValueError,
            "non-finite",
        ),
        (
            lambda data: data.assign(goldStandardSet=0),
            ValueError,
            "at least one positive",
        ),
        (
            lambda data: data.assign(goldStandardSet=2),
            ValueError,
            "only 0 and 1",
        ),
    ],
)
def test_fit_validation(mutation, error, message) -> None:
    model = BayesianSoftmaxModel(feature_names=FEATURES)
    with pytest.raises(error, match=message):
        model.fit(mutation(toy_frame()))


def test_constructor_and_unfitted_validation() -> None:
    with pytest.raises(ValueError, match="prior_scale"):
        BayesianSoftmaxModel(feature_names=FEATURES, prior_scale=0.0)
    with pytest.raises(ValueError, match="unique"):
        BayesianSoftmaxModel(feature_names=(FEATURES[0], FEATURES[0]))

    model = BayesianSoftmaxModel(feature_names=FEATURES)
    with pytest.raises(RuntimeError, match="fitted or loaded"):
        model.predict(toy_frame())

