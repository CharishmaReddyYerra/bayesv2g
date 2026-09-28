import numpy as np
import pytest

from bayesv2g.coloc import coloc_abf


def test_same_sharp_signal_supports_h4_and_identifies_shared_variant():
    gwas_beta = np.array([0.01, 0.02, 0.32, 0.01])
    qtl_beta = np.array([0.01, 0.03, -0.30, 0.02])
    standard_error = np.full(4, 0.03)

    result = coloc_abf(
        gwas_beta,
        standard_error,
        qtl_beta,
        standard_error,
        variant_ids=["rs1", "rs2", "rs3", "rs4"],
    )

    assert sum(result.posterior.values()) == pytest.approx(1.0)
    assert result.pp_h4 > 0.99
    assert result.top_shared_index == 2
    assert result.top_shared_variant == "rs3"
    assert result.gwas_pip.sum() == pytest.approx(1.0)
    assert result.qtl_pip.sum() == pytest.approx(1.0)
    assert result.shared_variant_pip.sum() == pytest.approx(1.0)
    assert result.shared_variant_posterior.sum() == pytest.approx(result.pp_h4)
    assert np.array_equal(result.shared_credible_set(0.95), np.array([2]))


def test_distinct_sharp_signals_support_h3_not_h4():
    standard_error = np.full(4, 0.025)
    result = coloc_abf(
        np.array([0.35, 0.01, 0.01, 0.01]),
        standard_error,
        np.array([0.01, -0.34, 0.01, 0.01]),
        standard_error,
    )

    assert result.posterior["H3"] > 0.99
    assert result.posterior["H4"] < 1e-6
    assert result.gwas_pip.argmax() == 0
    assert result.qtl_pip.argmax() == 1


def test_extreme_z_scores_remain_numerically_stable():
    result = coloc_abf(
        np.array([10.0, 0.1, -0.1]),
        np.array([0.001, 0.01, 0.01]),
        np.array([-9.0, 0.1, 0.1]),
        np.array([0.001, 0.01, 0.01]),
    )

    assert all(np.isfinite(list(result.posterior.values())))
    assert sum(result.posterior.values()) == pytest.approx(1.0)
    assert result.pp_h4 == pytest.approx(1.0)
    assert result.shared_variant_pip[0] == pytest.approx(1.0)


def test_one_variant_region_makes_h3_impossible():
    result = coloc_abf(
        np.array([0.1]),
        np.array([0.02]),
        np.array([0.1]),
        np.array([0.02]),
    )

    assert result.posterior["H3"] == 0.0
    assert result.log_hypothesis_weights["H3"] == -np.inf
    assert result.shared_variant_pip[0] == 1.0


@pytest.mark.parametrize(
    ("gwas_beta", "gwas_se", "qtl_beta", "qtl_se", "message"),
    [
        ([], [], [], [], "non-empty"),
        ([0.1, 0.2], [0.01], [0.1, 0.2], [0.01, 0.01], "same length"),
        ([0.1], [0.0], [0.1], [0.01], "strictly positive"),
        ([0.1], [0.01], [np.nan], [0.01], "finite"),
        ([[0.1]], [[0.01]], [[0.1]], [[0.01]], "one-dimensional"),
    ],
)
def test_invalid_effect_inputs_are_rejected(gwas_beta, gwas_se, qtl_beta, qtl_se, message):
    with pytest.raises(ValueError, match=message):
        coloc_abf(gwas_beta, gwas_se, qtl_beta, qtl_se)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"p1": 0.0},
        {"p2": np.nan},
        {"p12": 2e-4},
        {"p1": 0.8, "p2": 0.8, "p12": 0.1},
        {"gwas_prior_variance": 0.0},
        {"qtl_prior_variance": np.inf},
    ],
)
def test_invalid_priors_are_rejected(kwargs):
    with pytest.raises(ValueError):
        coloc_abf([0.1, 0.2], [0.01, 0.01], [0.1, 0.2], [0.01, 0.01], **kwargs)


@pytest.mark.parametrize("variant_ids", [["rs1"], ["rs1", "rs1"], ["rs1", None]])
def test_invalid_variant_ids_are_rejected(variant_ids):
    with pytest.raises(ValueError, match="variant_ids"):
        coloc_abf(
            [0.1, 0.2],
            [0.01, 0.01],
            [0.1, 0.2],
            [0.01, 0.01],
            variant_ids=variant_ids,
        )


def test_effect_sign_is_irrelevant_after_allele_harmonization():
    common = dict(gwas_se=np.array([0.02, 0.02]), qtl_se=np.array([0.02, 0.02]))
    positive = coloc_abf(np.array([0.2, 0.01]), qtl_beta=np.array([0.18, 0.01]), **common)
    opposite = coloc_abf(np.array([0.2, 0.01]), qtl_beta=np.array([-0.18, 0.01]), **common)

    assert positive.posterior == pytest.approx(opposite.posterior)
    assert positive.shared_variant_pip == pytest.approx(opposite.shared_variant_pip)
