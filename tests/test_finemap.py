import numpy as np
import pandas as pd
import pytest

from bayesv2g.finemap import credible_set, define_distance_loci, single_signal_pip


def test_single_signal_pips_are_normalized_and_favor_large_z_score():
    beta = np.array([0.02, 0.10, 0.01])
    standard_error = np.array([0.02, 0.02, 0.02])
    pip = single_signal_pip(beta, standard_error)
    assert np.isclose(pip.sum(), 1.0)
    assert int(np.argmax(pip)) == 1
    assert credible_set(pip, 0.95)[0] == 1


def test_distance_locus_definition_groups_nearby_hits():
    frame = pd.DataFrame(
        {
            "chromosome": [1, 1, 1, 2],
            "position": [100, 150, 2_000_000, 100],
            "p_value": [1e-10, 1e-9, 1e-12, 1e-11],
        }
    )
    loci = define_distance_loci(frame, window_bp=1_000)
    assert len(loci) == 3
    assert sorted(loci["member_count"]) == [1, 1, 2]


@pytest.mark.parametrize(
    ("beta", "standard_error"),
    [
        ([0.1], [0.01, 0.02]),
        ([np.nan], [0.01]),
        ([0.1], [0.0]),
        ([], []),
    ],
)
def test_single_signal_pip_rejects_invalid_effect_arrays(beta, standard_error):
    with pytest.raises(ValueError):
        single_signal_pip(np.asarray(beta), np.asarray(standard_error))


@pytest.mark.parametrize("pip", [[0.0, 0.0], [0.2, -0.1], [0.2, np.nan]])
def test_credible_set_rejects_invalid_probabilities(pip):
    with pytest.raises(ValueError):
        credible_set(np.asarray(pip))
