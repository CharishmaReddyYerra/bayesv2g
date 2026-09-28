"""Classic single-causal-variant colocalization with Wakefield ABFs.

The implementation follows the five-hypothesis calculation popularized by
``coloc.abf`` while keeping the result useful for an auditable regional
analysis.  It assumes that the four input arrays describe the same variants,
in the same order, and that each trait has at most one causal variant in the
region.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Hashable, Sequence

import numpy as np
from scipy.special import logsumexp

from bayesv2g.finemap import credible_set, wakefield_log_abf


HYPOTHESES = ("H0", "H1", "H2", "H3", "H4")


@dataclass(frozen=True, slots=True)
class ColocABFResult:
    """Posterior summary from :func:`coloc_abf`.

    ``shared_variant_pip`` is conditional on H4 (both traits share one causal
    variant).  ``shared_variant_posterior`` is unconditional and therefore
    sums to the H4 posterior probability.  The GWAS and QTL PIPs are each
    conditional on a single causal variant for that trait.
    """

    hypothesis_posteriors: dict[str, float]
    log_hypothesis_weights: dict[str, float]
    gwas_log_abf: np.ndarray
    qtl_log_abf: np.ndarray
    gwas_pip: np.ndarray
    qtl_pip: np.ndarray
    shared_variant_pip: np.ndarray
    shared_variant_posterior: np.ndarray
    variant_ids: tuple[Hashable, ...] | None = None

    @property
    def posterior(self) -> dict[str, float]:
        """Alias for the H0--H4 posterior mapping."""

        return self.hypothesis_posteriors

    @property
    def pp_h4(self) -> float:
        """Posterior probability that the two traits share a causal variant."""

        return self.hypothesis_posteriors["H4"]

    @property
    def top_shared_index(self) -> int:
        """Index of the most likely shared variant, conditional on H4."""

        return int(np.argmax(self.shared_variant_pip))

    @property
    def top_shared_variant(self) -> Hashable | int:
        """Identifier of the top shared variant, or its index if IDs were omitted."""

        if self.variant_ids is None:
            return self.top_shared_index
        return self.variant_ids[self.top_shared_index]

    def shared_credible_set(self, coverage: float = 0.95) -> np.ndarray:
        """Indices of the conditional-H4 shared-variant credible set."""

        return credible_set(self.shared_variant_pip, coverage=coverage)


def _validate_effect_arrays(
    gwas_beta: np.ndarray,
    gwas_se: np.ndarray,
    qtl_beta: np.ndarray,
    qtl_se: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    arrays = tuple(np.asarray(values, dtype=float) for values in (gwas_beta, gwas_se, qtl_beta, qtl_se))
    names = ("gwas_beta", "gwas_se", "qtl_beta", "qtl_se")
    for name, values in zip(names, arrays, strict=True):
        if values.ndim != 1 or values.size == 0:
            raise ValueError(f"{name} must be a non-empty one-dimensional array")
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{name} must contain only finite values")
    if len({values.size for values in arrays}) != 1:
        raise ValueError("GWAS and QTL beta/SE arrays must have the same length and aligned order")
    if np.any(arrays[1] <= 0) or np.any(arrays[3] <= 0):
        raise ValueError("gwas_se and qtl_se must be strictly positive")
    return arrays


def _validate_priors(p1: float, p2: float, p12: float) -> tuple[float, float, float]:
    priors = tuple(float(value) for value in (p1, p2, p12))
    if not all(np.isfinite(value) and 0 < value < 1 for value in priors):
        raise ValueError("p1, p2, and p12 must be finite probabilities strictly between 0 and 1")
    p1_value, p2_value, p12_value = priors
    if p12_value > min(p1_value, p2_value):
        raise ValueError("p12 cannot exceed either marginal association prior p1 or p2")
    if p1_value + p2_value - p12_value >= 1:
        raise ValueError("p1, p2, and p12 imply an invalid per-variant union probability")
    return priors


def _validate_variant_ids(
    variant_ids: Sequence[Hashable] | None,
    expected_length: int,
) -> tuple[Hashable, ...] | None:
    if variant_ids is None:
        return None
    identifiers = tuple(variant_ids)
    if len(identifiers) != expected_length:
        raise ValueError("variant_ids must have the same length as the effect arrays")
    if any(identifier is None for identifier in identifiers):
        raise ValueError("variant_ids cannot contain missing values")
    try:
        if len(set(identifiers)) != expected_length:
            raise ValueError("variant_ids must be unique within the region")
    except TypeError as error:
        raise ValueError("variant_ids must contain hashable scalar values") from error
    return identifiers


def _log_distinct_variant_sum(first: np.ndarray, second: np.ndarray) -> float:
    """Return log(sum_i!=j exp(first_i + second_j)) without subtraction.

    Prefix/suffix log-sums avoid the catastrophic cancellation in
    ``sum(BF1) * sum(BF2) - sum(BF1 * BF2)`` when both signals are sharp and
    peak at the same variant.
    """

    size = first.size
    if size == 1:
        return float("-inf")

    prefix = np.full(size + 1, -np.inf, dtype=float)
    suffix = np.full(size + 1, -np.inf, dtype=float)
    for index in range(size):
        prefix[index + 1] = np.logaddexp(prefix[index], second[index])
    for index in range(size - 1, -1, -1):
        suffix[index] = np.logaddexp(suffix[index + 1], second[index])
    second_excluding_index = np.logaddexp(prefix[:-1], suffix[1:])
    return float(logsumexp(first + second_excluding_index))


def _normalized_exp(log_values: np.ndarray) -> np.ndarray:
    return np.exp(log_values - logsumexp(log_values))


def coloc_abf(
    gwas_beta: np.ndarray,
    gwas_se: np.ndarray,
    qtl_beta: np.ndarray,
    qtl_se: np.ndarray,
    *,
    p1: float = 1e-4,
    p2: float = 1e-4,
    p12: float = 1e-5,
    gwas_prior_variance: float = 0.04,
    qtl_prior_variance: float = 0.04,
    variant_ids: Sequence[Hashable] | None = None,
) -> ColocABFResult:
    """Run classic five-hypothesis colocalization on aligned summary data.

    Parameters
    ----------
    gwas_beta, gwas_se, qtl_beta, qtl_se:
        One-dimensional arrays for the same variants in exactly the same
        order. Alleles must already be harmonized. Because Wakefield ABFs use
        squared z-scores, this function cannot detect allele-orientation
        mistakes.
    p1, p2, p12:
        Per-variant prior probabilities of association with trait 1, trait 2,
        or both traits. Defaults match common ``coloc.abf`` usage.
    gwas_prior_variance, qtl_prior_variance:
        Prior variances for true quantitative-trait effects.
    variant_ids:
        Optional unique identifiers retained in the result for audit output.

    Notes
    -----
    H0 means neither trait is associated; H1/H2 mean only the GWAS/QTL is
    associated; H3 means both are associated through different variants; and
    H4 means they share one causal variant. The method assumes at most one
    causal variant per trait in the analyzed region and does not use LD.
    """

    gwas_beta_values, gwas_se_values, qtl_beta_values, qtl_se_values = _validate_effect_arrays(
        gwas_beta, gwas_se, qtl_beta, qtl_se
    )
    p1_value, p2_value, p12_value = _validate_priors(p1, p2, p12)
    identifiers = _validate_variant_ids(variant_ids, gwas_beta_values.size)

    for name, value in (
        ("gwas_prior_variance", gwas_prior_variance),
        ("qtl_prior_variance", qtl_prior_variance),
    ):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and strictly positive")

    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        gwas_log_abf = wakefield_log_abf(gwas_beta_values, gwas_se_values, gwas_prior_variance)
        qtl_log_abf = wakefield_log_abf(qtl_beta_values, qtl_se_values, qtl_prior_variance)
    if not np.all(np.isfinite(gwas_log_abf)) or not np.all(np.isfinite(qtl_log_abf)):
        raise ValueError("Effect sizes and standard errors produced non-finite approximate Bayes factors")

    shared_log_abf = gwas_log_abf + qtl_log_abf
    log_weights = np.array(
        [
            0.0,
            np.log(p1_value) + logsumexp(gwas_log_abf),
            np.log(p2_value) + logsumexp(qtl_log_abf),
            np.log(p1_value)
            + np.log(p2_value)
            + _log_distinct_variant_sum(gwas_log_abf, qtl_log_abf),
            np.log(p12_value) + logsumexp(shared_log_abf),
        ],
        dtype=float,
    )
    hypothesis_values = _normalized_exp(log_weights)
    hypothesis_posteriors = {
        hypothesis: float(probability)
        for hypothesis, probability in zip(HYPOTHESES, hypothesis_values, strict=True)
    }
    shared_variant_pip = _normalized_exp(shared_log_abf)

    return ColocABFResult(
        hypothesis_posteriors=hypothesis_posteriors,
        log_hypothesis_weights={
            hypothesis: float(weight)
            for hypothesis, weight in zip(HYPOTHESES, log_weights, strict=True)
        },
        gwas_log_abf=gwas_log_abf.copy(),
        qtl_log_abf=qtl_log_abf.copy(),
        gwas_pip=_normalized_exp(gwas_log_abf),
        qtl_pip=_normalized_exp(qtl_log_abf),
        shared_variant_pip=shared_variant_pip,
        shared_variant_posterior=shared_variant_pip * hypothesis_posteriors["H4"],
        variant_ids=identifiers,
    )


# A descriptive alias for callers that want to make the model assumption
# explicit at the call site.
single_signal_coloc_abf = coloc_abf

