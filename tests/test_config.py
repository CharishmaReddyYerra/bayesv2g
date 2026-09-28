import json
from importlib.resources import files
from pathlib import Path

import pytest

from bayesv2g.config import load_config


ROOT = Path(__file__).resolve().parents[1]


def test_packaged_and_checkout_default_configs_match():
    checkout = json.loads((ROOT / "config/demo_ldl.json").read_text(encoding="utf-8"))
    packaged = json.loads(
        files("bayesv2g").joinpath("default_config.json").read_text(encoding="utf-8")
    )
    assert packaged == checkout
    assert load_config()["study_id"] == "GCST90019512"


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("unresolved_mass", 1.1, "unresolved_mass"),
        ("prior_scale", 0, "prior_scale"),
        ("posterior_draws", 0, "posterior_draws"),
        ("features", ["vepMean", "vepMean"], "unique"),
    ],
)
def test_invalid_config_values_are_rejected(tmp_path, field, value, message):
    config = json.loads((ROOT / "config/demo_ldl.json").read_text(encoding="utf-8"))
    config[field] = value
    path = tmp_path / "invalid.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        load_config(path)
