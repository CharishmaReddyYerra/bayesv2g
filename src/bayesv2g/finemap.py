"""Small, dependency-light locus definition and single-signal fine-mapping helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import logsumexp


def define_distance_loci(
    variants: pd.DataFrame,
    *,
    p_threshold: float = 5e-8,
    window_bp: int = 500_000,
) -> pd.DataFrame:
    """Greedy distance clumping for an exploratory fallback.

    This deliberately does not claim to be LD clumping. The production LDL
    demonstration consumes externally fine-mapped Open Targets signals.
    """

    required = {"chromosome", "position", "p_value"}
    missing = required.difference(variants.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    significant = variants.loc[variants["p_value"] <= p_threshold].copy()
    significant = significant.sort_values("p_value")
    assigned: set[int] = set()
    output: list[pd.Series] = []
    for index, lead in significant.iterrows():
        if index in assigned:
            continue
        same_region = (
            significant["chromosome"].eq(lead["chromosome"])
            & ((significant["position"] - lead["position"]).abs() <= window_bp)
        )
        members = significant.index[same_region]
        assigned.update(members.tolist())
        record = lead.copy()
        record["locus_id"] = f"chr{lead['chromosome']}:{int(lead['position'])}"
        record["locus_start"] = max(1, int(significant.loc[members, "position"].min()) - window_bp)
        record["locus_end"] = int(significant.loc[members, "position"].max()) + window_bp
        record["member_count"] = len(members)
        record["definition_method"] = "greedy distance window (not LD clumping)"
        output.append(record)
    return pd.DataFrame(output).reset_index(drop=True)


def wakefield_log_abf(beta: np.ndarray, standard_error: np.ndarray, prior_variance: float = 0.04) -> np.ndarray:
    """Wakefield approximate log Bayes factor for a quantitative-trait effect."""

    beta = np.asarray(beta, dtype=float)
    standard_error = np.asarray(standard_error, dtype=float)
    if beta.ndim != 1 or standard_error.ndim != 1 or len(beta) == 0:
        raise ValueError("beta and standard_error must be non-empty one-dimensional arrays")
    if beta.shape != standard_error.shape:
        raise ValueError("beta and standard_error must have the same length")
    if not np.all(np.isfinite(beta)) or not np.all(np.isfinite(standard_error)):
        raise ValueError("beta and standard_error must contain only finite values")
    if not np.isfinite(prior_variance) or np.any(standard_error <= 0) or prior_variance <= 0:
        raise ValueError("Standard errors and prior_variance must be finite and positive")
    variance = standard_error**2
    z_squared = (beta / standard_error) ** 2
    shrinkage = prior_variance / (variance + prior_variance)
    return 0.5 * (np.log1p(-shrinkage) + z_squared * shrinkage)


def single_signal_pip(
    beta: np.ndarray,
    standard_error: np.ndarray,
    *,
    prior_variance: float = 0.04,
) -> np.ndarray:
    """Normalize approximate Bayes factors under a one-causal-variant model."""

    log_bf = wakefield_log_abf(beta, standard_error, prior_variance)
    return np.exp(log_bf - logsumexp(log_bf))


def credible_set(pip: np.ndarray, coverage: float = 0.95) -> np.ndarray:
    """Return indices in the smallest descending-PIP set reaching coverage."""

    probabilities = np.asarray(pip, dtype=float)
    if probabilities.ndim != 1 or len(probabilities) == 0:
        raise ValueError("pip must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(probabilities)) or np.any(probabilities < 0):
        raise ValueError("pip must contain finite, non-negative values")
    total = probabilities.sum()
    if not np.isfinite(total) or total <= 0:
        raise ValueError("pip must contain positive total probability mass")
    if not np.isfinite(coverage) or not 0 < coverage <= 1:
        raise ValueError("coverage must lie in (0, 1]")
    probabilities = probabilities / total
    order = np.argsort(-probabilities)
    stop = int(np.searchsorted(np.cumsum(probabilities[order]), coverage, side="left")) + 1
    return order[:stop]
