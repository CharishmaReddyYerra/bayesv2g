import json
from pathlib import Path

import numpy as np
import pandas as pd

from bayesv2g.pipeline import run_pipeline


ROOT = Path(__file__).resolve().parents[1]


def _first_loci(frame: pd.DataFrame, count: int) -> pd.DataFrame:
    locus_ids = frame["studyLocusId"].drop_duplicates().head(count)
    return frame.loc[frame["studyLocusId"].isin(locus_ids)].reset_index(drop=True)


def test_offline_pipeline_smoke_with_real_snapshots(tmp_path):
    source_training = ROOT / "data/raw/open_targets_l2g_26.06"
    training = tmp_path / "training"
    training.mkdir()
    _first_loci(pd.read_parquet(source_training / "train.parquet"), 50).to_parquet(
        training / "train.parquet",
        index=False,
    )
    _first_loci(pd.read_parquet(source_training / "test.parquet"), 30).to_parquet(
        training / "test.parquet",
        index=False,
    )

    config = json.loads((ROOT / "config/demo_ldl.json").read_text(encoding="utf-8"))
    config.update(
        {
            "training_directory": str(training),
            "open_targets_snapshot": str(
                ROOT / "data/snapshots/open_targets_26.09_ldl_loci.json"
            ),
            "sort1_case_directory": str(ROOT / "data/snapshots/sort1_case"),
            "posterior_draws": 20,
        }
    )
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    output = tmp_path / "results"

    paths = run_pipeline(config_path=config_path, output_directory=output)

    assert paths["report"].exists()
    assert (output / "sort1_colocalization.csv").exists()
    assert (output / "paired_top1_comparison.csv").exists()
    assert paths["checksums"].exists()
    rankings = pd.read_csv(paths["prioritized_genes"])
    for _, locus in rankings.groupby("studyLocusId"):
        total = locus["posterior_mean"].sum() + locus["unresolved_probability"].iloc[0]
        assert np.isclose(total, 1.0)
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    assert manifest["config"]["posterior_draws"] == 20
    assert manifest["sort1_case_study"]["regional_colocalization_count"] == 4
    assert manifest["source_tls_verification"] == {
        "open_targets_snapshot": True,
        "sort1_case_json_and_metadata": True,
    }
