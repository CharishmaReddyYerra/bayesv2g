from pathlib import Path

import pandas as pd

from bayesv2g.case_study import analyse_sort1_case


ROOT = Path(__file__).resolve().parents[1]


def test_real_sort1_case_study(tmp_path):
    rankings = pd.DataFrame(
        {
            "studyLocusId": ["5eb8878b948ed248f0832a6a1b27d324"] * 3,
            "gene_symbol": ["SORT1", "CELSR2", "PSRC1"],
            "rank": [1, 2, 3],
            "posterior_mean": [0.6, 0.3, 0.1],
            "posterior_lower": [0.5, 0.2, 0.05],
            "posterior_upper": [0.7, 0.4, 0.2],
            "conditional_posterior_mean": [0.6, 0.3, 0.1],
            "l2g_score_external": [0.98, 0.95, 0.22],
            "eQtlColocH4Maximum": [1.0, 1.0, 1.0],
            "pQtlColocH4Maximum": [0.99, 1.0, 0.0],
            "sQtlColocH4Maximum": [1.0, 1.0, 0.99],
            "vepMean": [0.2, 0.37, 0.2],
            "e2gMean": [0.77, 0.98, 0.98],
        }
    )
    result = analyse_sort1_case(ROOT / "data/snapshots/sort1_case", tmp_path, rankings)
    assert result["gwas_variant_count"] == 5772
    assert result["lead_variant"] == "rs12740374"
    assert result["lead_variant_pip"] > 0.9
    assert set(result["genes"]["geneSymbol"]) == {"SORT1", "CELSR2", "PSRC1"}
