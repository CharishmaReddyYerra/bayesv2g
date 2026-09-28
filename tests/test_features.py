import json
from pathlib import Path

from bayesv2g.features import evidence_coverage, flatten_open_targets


ROOT = Path(__file__).resolve().parents[1]


def test_real_snapshot_flattens_to_normalized_tables():
    snapshot = json.loads((ROOT / "data/snapshots/open_targets_26.09_ldl_loci.json").read_text())
    candidates, loci, variants = flatten_open_targets(snapshot)
    assert len(loci) == 10
    assert len(candidates) == int(loci["candidate_count"].sum())
    assert not variants.empty
    assert {"gene_symbol", "eQtlColocH4Maximum", "distance_to_tss_bp"}.issubset(candidates)
    assert evidence_coverage(candidates).between(0, 1).all()
    assert loci["pip_coverage_95"].between(0.94, 1.01).all()
    sort_region = candidates.loc[candidates["lead_variant"] == "rs12740374"].set_index("gene_symbol")
    assert int(sort_region.loc["CELSR2", "distance_to_tss_bp"]) == 25_429
    assert int(sort_region.loc["PSRC1", "distance_to_tss_bp"]) == 9_139
    assert int(sort_region.loc["SORT1", "distance_to_tss_bp"]) == 122_999
