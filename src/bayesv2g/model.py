"""Conditional Bayesian softmax model for locus-to-gene prioritisation.

The model treats the candidate genes at a study locus as a single choice set.
It fits a conditional multinomial logistic regression with a zero-centred
Gaussian prior and uses a Laplace approximation around the maximum a
posteriori (MAP) estimate for posterior prediction.

Only NumPy, pandas, and SciPy are used so the model remains lightweight and
can be fitted to the public Open Targets L2G feature matrix on a laptop.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import OptimizeResult, minimize


DEFAULT_FEATURES: tuple[str, ...] = (
    "distanceTssMean",
    "eQtlColocH4Maximum",
    "pQtlColocH4Maximum",
    "sQtlColocH4Maximum",
    "vepMean",
    "e2gMean",
)

_MODEL_TYPE = "conditional_bayesian_softmax"
_SERIALIZATION_VERSION = 1


class BayesianSoftmaxModel:
    """Bayesian conditional softmax model with a Laplace posterior.

    Parameters
    ----------
    feature_names:
        Numeric evidence columns used in the linear predictor.
    group_column:
        Column identifying the competing candidate genes at one study locus.
    label_column:
        Binary gold-standard label. If a locus contains multiple positives,
        its unit target mass is divided equally among them.
    prior_scale:
        Standard deviation of the independent Normal(0, prior_scale) priors.
    n_posterior_draws:
        Default number of multivariate-normal coefficient draws used by
        :meth:`predict`.
    random_seed:
        Default seed used for posterior coefficient draws.
    max_iter:
        Maximum L-BFGS-B iterations.
    tolerance:
        Optimizer convergence tolerance.

    Notes
    -----
    There is deliberately no intercept. A common intercept cancels within
    every groupwise softmax and is therefore not identifiable.
    """

    def __init__(
        self,
        feature_names: Sequence[str] = DEFAULT_FEATURES,
        *,
        group_column: str = "studyLocusId",
        label_column: str = "goldStandardSet",
        prior_scale: float = 0.75,
        n_posterior_draws: int = 1_000,
        random_seed: int = 42,
        max_iter: int = 1_000,
        tolerance: float = 1e-9,
    ) -> None:
        self.feature_names = self._validate_feature_names(feature_names)
        self.group_column = self._validate_column_name(
            group_column, "group_column"
        )
        self.label_column = self._validate_column_name(
            label_column, "label_column"
        )
        if self.group_column in self.feature_names:
            raise ValueError("group_column must not also be a feature")
        if self.label_column in self.feature_names:
            raise ValueError("label_column must not also be a feature")
        if self.group_column == self.label_column:
            raise ValueError("group_column and label_column must be different")

        self.prior_scale = self._validate_positive_float(
            prior_scale, "prior_scale"
        )
        self.n_posterior_draws = self._validate_positive_int(
            n_posterior_draws, "n_posterior_draws"
        )
        self.random_seed = self._validate_seed(random_seed)
        self.max_iter = self._validate_positive_int(max_iter, "max_iter")
        self.tolerance = self._validate_positive_float(tolerance, "tolerance")

        self.coef_: np.ndarray | None = None
        self.hessian_: np.ndarray | None = None
        self.covariance_: np.ndarray | None = None
        self.fit_metadata_: dict[str, Any] | None = None
        self.optimization_result_: OptimizeResult | None = None

    @property
    def is_fitted(self) -> bool:
        """Whether MAP coefficients and a Laplace covariance are available."""

        return (
            self.coef_ is not None
            and self.hessian_ is not None
            and self.covariance_ is not None
        )

    def fit(self, frame: pd.DataFrame) -> "BayesianSoftmaxModel":
        """Fit the conditional model and its Laplace approximation.

        Each locus contributes one cross-entropy term, irrespective of its
        number of candidate genes. Loci with multiple positive labels use a
        target distribution that is uniform over their positive genes.
        """

        x, group_codes, n_groups = self._prepare_frame(
            frame, require_labels=True
        )
        target = self._make_targets(
            frame[self.label_column], group_codes, n_groups
        )

        def objective(beta: np.ndarray) -> tuple[float, np.ndarray]:
            return self._objective_and_gradient(
                beta, x, group_codes, n_groups, target
            )

        result = minimize(
            objective,
            x0=np.zeros(len(self.feature_names), dtype=np.float64),
            method="L-BFGS-B",
            jac=True,
            tol=self.tolerance,
            options={"maxiter": self.max_iter, "maxls": 50},
        )
        if not result.success:
            raise RuntimeError(
                "L-BFGS-B failed to fit the conditional softmax model: "
                f"{result.message}"
            )
        if not np.all(np.isfinite(result.x)):
            raise RuntimeError("optimizer returned non-finite MAP coefficients")

        coef = np.asarray(result.x, dtype=np.float64)
        hessian = self._analytic_hessian(coef, x, group_codes, n_groups)
        covariance = self._invert_positive_definite(hessian)

        self.coef_ = coef
        self.hessian_ = hessian
        self.covariance_ = covariance
        self.optimization_result_ = result
        positive_count = int(np.count_nonzero(np.asarray(target) > 0.0))
        self.fit_metadata_ = {
            "n_rows": int(len(frame)),
            "n_groups": int(n_groups),
            "n_positive_rows": positive_count,
            "objective_at_map": float(result.fun),
            "optimizer_iterations": int(result.nit),
            "optimizer_evaluations": int(result.nfev),
            "optimizer_message": str(result.message),
        }
        return self

    def predict(
        self,
        frame: pd.DataFrame,
        *,
        n_draws: int | None = None,
        random_seed: int | None = None,
    ) -> pd.DataFrame:
        """Return posterior target probabilities and 5%--95% intervals.

        A fresh copy of ``frame`` is returned with three appended columns:
        ``posterior_mean``, ``posterior_lower``, and ``posterior_upper``.
        Probabilities are normalized separately within every study locus.
        Calling this method repeatedly with the same seed is deterministic.
        """

        self._require_fitted()
        x, group_codes, n_groups = self._prepare_frame(
            frame, require_labels=False
        )
        draws = self._coefficient_draws(
            self.n_posterior_draws if n_draws is None else n_draws,
            self.random_seed if random_seed is None else random_seed,
        )

        posterior_mean = np.empty(len(frame), dtype=np.float64)
        posterior_lower = np.empty(len(frame), dtype=np.float64)
        posterior_upper = np.empty(len(frame), dtype=np.float64)

        # Work one choice set at a time. This bounds peak memory by
        # n_draws * largest_group_size instead of n_draws * n_rows.
        for group_index in range(n_groups):
            row_index = np.flatnonzero(group_codes == group_index)
            logits = draws @ x[row_index].T
            logits -= np.max(logits, axis=1, keepdims=True)
            probabilities = np.exp(logits)
            probabilities /= np.sum(probabilities, axis=1, keepdims=True)

            posterior_mean[row_index] = np.mean(probabilities, axis=0)
            bounds = np.quantile(probabilities, (0.05, 0.95), axis=0)
            posterior_lower[row_index] = bounds[0]
            posterior_upper[row_index] = bounds[1]

        result = frame.copy()
        result["posterior_mean"] = posterior_mean
        result["posterior_lower"] = posterior_lower
        result["posterior_upper"] = posterior_upper
        return result

    def predict_map(self, frame: pd.DataFrame) -> np.ndarray:
        """Return group-normalized probabilities at the MAP coefficients."""

        self._require_fitted()
        x, group_codes, n_groups = self._prepare_frame(
            frame, require_labels=False
        )
        return self._groupwise_softmax(
            x @ self.coef_, group_codes, n_groups
        )[0]

    def save(self, path: str | Path) -> None:
        """Serialize the fitted model to a portable JSON file."""

        self._require_fitted()
        destination = Path(path)
        if destination.exists() and destination.is_dir():
            raise IsADirectoryError(f"model path is a directory: {destination}")
        payload = {
            "model_type": _MODEL_TYPE,
            "serialization_version": _SERIALIZATION_VERSION,
            "feature_names": list(self.feature_names),
            "group_column": self.group_column,
            "label_column": self.label_column,
            "prior_scale": self.prior_scale,
            "n_posterior_draws": self.n_posterior_draws,
            "random_seed": self.random_seed,
            "max_iter": self.max_iter,
            "tolerance": self.tolerance,
            "coef": self.coef_.tolist(),
            "hessian": self.hessian_.tolist(),
            "covariance": self.covariance_.tolist(),
            "fit_metadata": self.fit_metadata_,
        }
        with destination.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")

    @classmethod
    def load(cls, path: str | Path) -> "BayesianSoftmaxModel":
        """Load a fitted model written by :meth:`save`."""

        source = Path(path)
        with source.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        if not isinstance(payload, dict):
            raise ValueError("model JSON must contain an object")
        if payload.get("model_type") != _MODEL_TYPE:
            raise ValueError("JSON is not a BayesV2G conditional softmax model")
        if payload.get("serialization_version") != _SERIALIZATION_VERSION:
            raise ValueError(
                "unsupported model serialization version: "
                f"{payload.get('serialization_version')!r}"
            )

        required = {
            "feature_names",
            "group_column",
            "label_column",
            "prior_scale",
            "n_posterior_draws",
            "random_seed",
            "max_iter",
            "tolerance",
            "coef",
            "hessian",
            "covariance",
        }
        missing = sorted(required.difference(payload))
        if missing:
            raise ValueError(f"model JSON is missing fields: {', '.join(missing)}")

        model = cls(
            feature_names=payload["feature_names"],
            group_column=payload["group_column"],
            label_column=payload["label_column"],
            prior_scale=payload["prior_scale"],
            n_posterior_draws=payload["n_posterior_draws"],
            random_seed=payload["random_seed"],
            max_iter=payload["max_iter"],
            tolerance=payload["tolerance"],
        )
        n_features = len(model.feature_names)
        model.coef_ = cls._validated_loaded_array(
            payload["coef"], (n_features,), "coef"
        )
        model.hessian_ = cls._validated_loaded_array(
            payload["hessian"], (n_features, n_features), "hessian"
        )
        model.covariance_ = cls._validated_loaded_array(
            payload["covariance"], (n_features, n_features), "covariance"
        )
        if not np.allclose(model.hessian_, model.hessian_.T, atol=1e-10):
            raise ValueError("loaded Hessian is not symmetric")
        if not np.allclose(model.covariance_, model.covariance_.T, atol=1e-10):
            raise ValueError("loaded covariance is not symmetric")
        if np.min(np.linalg.eigvalsh(model.hessian_)) <= 0.0:
            raise ValueError("loaded Hessian is not positive definite")
        if np.min(np.linalg.eigvalsh(model.covariance_)) <= 0.0:
            raise ValueError("loaded covariance is not positive definite")
        metadata = payload.get("fit_metadata")
        if metadata is not None and not isinstance(metadata, dict):
            raise ValueError("fit_metadata must be an object or null")
        model.fit_metadata_ = metadata
        return model

    def _prepare_frame(
        self, frame: pd.DataFrame, *, require_labels: bool
    ) -> tuple[np.ndarray, np.ndarray, int]:
        if not isinstance(frame, pd.DataFrame):
            raise TypeError("frame must be a pandas DataFrame")
        if frame.empty:
            raise ValueError("frame must contain at least one candidate gene")

        required = [self.group_column, *self.feature_names]
        if require_labels:
            required.append(self.label_column)
        missing = [column for column in required if column not in frame.columns]
        if missing:
            raise ValueError(f"missing required columns: {', '.join(missing)}")

        if frame[self.group_column].isna().any():
            raise ValueError(f"{self.group_column} contains missing values")
        try:
            group_codes, group_values = pd.factorize(
                frame[self.group_column], sort=False
            )
        except TypeError as exc:
            raise TypeError(
                f"{self.group_column} values must be hashable"
            ) from exc
        if np.any(group_codes < 0):
            raise ValueError(f"{self.group_column} contains invalid values")

        non_numeric = [
            column
            for column in self.feature_names
            if not pd.api.types.is_numeric_dtype(frame[column].dtype)
        ]
        if non_numeric:
            raise TypeError(
                "feature columns must be numeric: " + ", ".join(non_numeric)
            )
        x = frame.loc[:, self.feature_names].to_numpy(
            dtype=np.float64, copy=True
        )
        if not np.all(np.isfinite(x)):
            raise ValueError("feature matrix contains missing or non-finite values")
        return x, np.asarray(group_codes, dtype=np.int64), len(group_values)

    @staticmethod
    def _make_targets(
        labels: pd.Series, group_codes: np.ndarray, n_groups: int
    ) -> np.ndarray:
        if not pd.api.types.is_numeric_dtype(labels.dtype):
            raise TypeError("gold-standard labels must be numeric or boolean")
        values = labels.to_numpy(dtype=np.float64, copy=True)
        if not np.all(np.isfinite(values)):
            raise ValueError("gold-standard labels contain missing or non-finite values")
        if not np.all((values == 0.0) | (values == 1.0)):
            raise ValueError("gold-standard labels must contain only 0 and 1")

        positive = values == 1.0
        positive_counts = np.bincount(
            group_codes, weights=positive.astype(np.float64), minlength=n_groups
        )
        groups_without_positive = np.flatnonzero(positive_counts == 0.0)
        if groups_without_positive.size:
            raise ValueError(
                "every study locus must contain at least one positive "
                "gold-standard gene"
            )
        target = np.zeros_like(values, dtype=np.float64)
        target[positive] = 1.0 / positive_counts[group_codes[positive]]
        return target

    def _objective_and_gradient(
        self,
        beta: np.ndarray,
        x: np.ndarray,
        group_codes: np.ndarray,
        n_groups: int,
        target: np.ndarray,
    ) -> tuple[float, np.ndarray]:
        logits = x @ beta
        probabilities, log_probabilities = self._groupwise_softmax(
            logits, group_codes, n_groups
        )
        prior_precision = 1.0 / (self.prior_scale**2)
        objective = -float(target @ log_probabilities)
        objective += 0.5 * prior_precision * float(beta @ beta)
        gradient = x.T @ (probabilities - target) + prior_precision * beta
        return objective, np.asarray(gradient, dtype=np.float64)

    def _analytic_hessian(
        self,
        beta: np.ndarray,
        x: np.ndarray,
        group_codes: np.ndarray,
        n_groups: int,
    ) -> np.ndarray:
        probabilities = self._groupwise_softmax(
            x @ beta, group_codes, n_groups
        )[0]
        weighted_second_moment = x.T @ (x * probabilities[:, None])
        group_means = np.zeros(
            (n_groups, x.shape[1]), dtype=np.float64
        )
        np.add.at(group_means, group_codes, x * probabilities[:, None])
        hessian = weighted_second_moment - group_means.T @ group_means
        hessian += np.eye(x.shape[1]) / (self.prior_scale**2)
        return 0.5 * (hessian + hessian.T)

    @staticmethod
    def _groupwise_softmax(
        logits: np.ndarray, group_codes: np.ndarray, n_groups: int
    ) -> tuple[np.ndarray, np.ndarray]:
        maxima = np.full(n_groups, -np.inf, dtype=np.float64)
        np.maximum.at(maxima, group_codes, logits)
        shifted = logits - maxima[group_codes]
        exponentials = np.exp(shifted)
        normalizers = np.zeros(n_groups, dtype=np.float64)
        np.add.at(normalizers, group_codes, exponentials)
        log_probabilities = shifted - np.log(normalizers[group_codes])
        return np.exp(log_probabilities), log_probabilities

    @staticmethod
    def _invert_positive_definite(hessian: np.ndarray) -> np.ndarray:
        try:
            factor = cho_factor(hessian, lower=True, check_finite=True)
            covariance = cho_solve(
                factor, np.eye(hessian.shape[0]), check_finite=True
            )
        except np.linalg.LinAlgError as exc:
            raise RuntimeError(
                "MAP Hessian is not positive definite; cannot form Laplace "
                "covariance"
            ) from exc
        return 0.5 * (covariance + covariance.T)

    def _coefficient_draws(
        self, n_draws: int, random_seed: int
    ) -> np.ndarray:
        n_draws = self._validate_positive_int(n_draws, "n_draws")
        random_seed = self._validate_seed(random_seed)
        try:
            root = np.linalg.cholesky(self.covariance_)
        except np.linalg.LinAlgError as exc:
            raise RuntimeError(
                "posterior covariance is not positive definite"
            ) from exc
        rng = np.random.default_rng(random_seed)
        standard_normal = rng.standard_normal(
            (n_draws, len(self.feature_names))
        )
        return self.coef_[None, :] + standard_normal @ root.T

    def _require_fitted(self) -> None:
        if not self.is_fitted:
            raise RuntimeError("model must be fitted or loaded before prediction")

    @staticmethod
    def _validate_feature_names(feature_names: Sequence[str]) -> tuple[str, ...]:
        if isinstance(feature_names, str):
            raise TypeError("feature_names must be a sequence of column names")
        try:
            names = tuple(feature_names)
        except TypeError as exc:
            raise TypeError(
                "feature_names must be a sequence of column names"
            ) from exc
        if not names:
            raise ValueError("feature_names must not be empty")
        if any(not isinstance(name, str) or not name for name in names):
            raise ValueError("every feature name must be a non-empty string")
        if len(set(names)) != len(names):
            raise ValueError("feature_names must be unique")
        return names

    @staticmethod
    def _validate_column_name(value: str, parameter: str) -> str:
        if not isinstance(value, str) or not value:
            raise ValueError(f"{parameter} must be a non-empty string")
        return value

    @staticmethod
    def _validate_positive_float(value: float, parameter: str) -> float:
        if isinstance(value, bool):
            raise TypeError(f"{parameter} must be a finite positive number")
        try:
            parsed = float(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(
                f"{parameter} must be a finite positive number"
            ) from exc
        if not np.isfinite(parsed) or parsed <= 0.0:
            raise ValueError(f"{parameter} must be a finite positive number")
        return parsed

    @staticmethod
    def _validate_positive_int(value: int, parameter: str) -> int:
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise TypeError(f"{parameter} must be a positive integer")
        parsed = int(value)
        if parsed <= 0:
            raise ValueError(f"{parameter} must be a positive integer")
        return parsed

    @staticmethod
    def _validate_seed(value: int) -> int:
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise TypeError("random seed must be a non-negative integer")
        parsed = int(value)
        if parsed < 0:
            raise ValueError("random seed must be a non-negative integer")
        return parsed

    @staticmethod
    def _validated_loaded_array(
        value: Any, expected_shape: tuple[int, ...], field: str
    ) -> np.ndarray:
        try:
            array = np.asarray(value, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"loaded {field} is not numeric") from exc
        if array.shape != expected_shape:
            raise ValueError(
                f"loaded {field} has shape {array.shape}; expected "
                f"{expected_shape}"
            )
        if not np.all(np.isfinite(array)):
            raise ValueError(f"loaded {field} contains non-finite values")
        return array


# Descriptive aliases make the implementation easy to discover while keeping
# one canonical class for serialization and documentation.
ConditionalBayesianSoftmax = BayesianSoftmaxModel
ConditionalSoftmaxModel = BayesianSoftmaxModel


__all__ = [
    "BayesianSoftmaxModel",
    "ConditionalBayesianSoftmax",
    "ConditionalSoftmaxModel",
    "DEFAULT_FEATURES",
]
