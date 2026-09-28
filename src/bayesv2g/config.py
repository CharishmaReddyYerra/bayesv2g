"""Configuration loading and validation."""

from __future__ import annotations

import json
from importlib.resources import files
from pathlib import Path
from typing import Any


def project_root() -> Path:
    """Return the source checkout root, or the working directory when installed."""

    checkout = Path(__file__).resolve().parents[2]
    if (checkout / "pyproject.toml").exists():
        return checkout
    return Path.cwd().resolve()


def default_config_path() -> Path:
    checkout_config = project_root() / "config" / "demo_ldl.json"
    if checkout_config.exists():
        return checkout_config
    return Path(str(files("bayesv2g").joinpath("default_config.json")))


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    config_path = Path(path) if path is not None else default_config_path()
    with config_path.open(encoding="utf-8") as handle:
        config = json.load(handle)

    required = {
        "study_id",
        "trait",
        "genome_build",
        "open_targets_endpoint",
        "open_targets_release",
        "training_revision",
        "training_release",
        "prior_scale",
        "posterior_draws",
        "random_seed",
        "unresolved_mass",
        "features",
        "demo_loci",
    }
    missing = required.difference(config)
    if missing:
        raise ValueError(f"Configuration is missing keys: {sorted(missing)}")
    string_keys = (
        "study_id",
        "trait",
        "genome_build",
        "open_targets_endpoint",
        "open_targets_release",
        "training_revision",
        "training_release",
    )
    for key in string_keys:
        if not isinstance(config[key], str) or not config[key].strip():
            raise ValueError(f"Configuration field {key!r} must be a non-empty string")

    features = config["features"]
    if (
        not isinstance(features, list)
        or not features
        or any(not isinstance(item, str) or not item.strip() for item in features)
    ):
        raise ValueError("Configuration features must be a non-empty list of names")
    if len(features) != len(set(features)):
        raise ValueError("Configuration features must be unique")

    if not isinstance(config["prior_scale"], (int, float)) or config["prior_scale"] <= 0:
        raise ValueError("Configuration prior_scale must be positive")
    if not isinstance(config["posterior_draws"], int) or isinstance(config["posterior_draws"], bool) or config["posterior_draws"] <= 0:
        raise ValueError("Configuration posterior_draws must be a positive integer")
    if not isinstance(config["random_seed"], int) or isinstance(config["random_seed"], bool) or not 0 <= config["random_seed"] < 2**32:
        raise ValueError("Configuration random_seed must be an integer in [0, 2^32)")
    if not isinstance(config["unresolved_mass"], (int, float)) or isinstance(config["unresolved_mass"], bool) or not 0 <= config["unresolved_mass"] < 1:
        raise ValueError("Configuration unresolved_mass must lie in [0, 1)")

    if not isinstance(config["demo_loci"], list) or not config["demo_loci"]:
        raise ValueError("Configuration must contain at least one demo locus")
    locus_ids: list[str] = []
    for index, locus in enumerate(config["demo_loci"]):
        if not isinstance(locus, dict):
            raise ValueError(f"demo_loci[{index}] must be an object")
        for key in ("label", "study_locus_id"):
            if not isinstance(locus.get(key), str) or not locus[key].strip():
                raise ValueError(f"demo_loci[{index}].{key} must be a non-empty string")
        locus_ids.append(locus["study_locus_id"])
    if len(locus_ids) != len(set(locus_ids)):
        raise ValueError("Configured study_locus_id values must be unique")

    for key in ("open_targets_snapshot", "training_directory", "sort1_case_directory"):
        if key in config and (not isinstance(config[key], str) or not config[key].strip()):
            raise ValueError(f"Configuration field {key!r} must be a non-empty path string")
    return config
