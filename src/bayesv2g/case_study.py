"""Deep raw-data audit for the rs12740374 LDL locus."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from bayesv2g.coloc import coloc_abf
from bayesv2g.data_sources import sha256_file
from bayesv2g.finemap import credible_set, single_signal_pip, wakefield_log_abf


SORT1_STUDY_LOCUS_ID = "5eb8878b948ed248f0832a6a1b27d324"
GENE_SYMBOLS = {
    "ENSG00000134243": "SORT1",
    "ENSG00000143126": "CELSR2",
    "ENSG00000134222": "PSRC1",
}


def _snv_key(position: pd.Series, first: pd.Series, second: pd.Series) -> pd.Series:
    allele_one = first.astype(str).str.upper()
    allele_two = second.astype(str).str.upper()
    low = allele_one.where(allele_one <= allele_two, allele_two)
    high = allele_two.where(allele_one <= allele_two, allele_one)
    return position.astype(int).astype(str) + ":" + low + "/" + high


def _deduplicate_aligned_variants(
    frame: pd.DataFrame,
    keys: list[str],
    numeric_columns: list[str],
    label: str,
) -> pd.DataFrame:
    duplicates = frame.duplicated(keys, keep=False)
    if duplicates.any():
        inconsistent = (
            frame.loc[duplicates]
            .groupby(keys, dropna=False)[numeric_columns]
            .nunique(dropna=False)
            .gt(1)
            .any(axis=1)
        )
        if inconsistent.any():
            raise ValueError(f"{label} contains conflicting duplicate variant records")
    return frame.drop_duplicates(keys).copy()


def _regional_colocalization(
    gwas: pd.DataFrame,
    qtl_path: Path,
    *,
    evidence_type: str,
    dataset_id: str,
) -> tuple[list[dict[str, Any]], list[pd.DataFrame]]:
    qtl = pd.read_csv(qtl_path, sep="\t")
    required = {"gene_id", "position", "ref", "alt", "type", "beta", "se", "rsid"}
    missing = required.difference(qtl.columns)
    if missing:
        raise ValueError(f"{dataset_id} QTL snapshot is missing columns: {sorted(missing)}")
    qtl = qtl.loc[
        qtl["type"].eq("SNP")
        & qtl["ref"].astype(str).str.fullmatch("[ACGTacgt]")
        & qtl["alt"].astype(str).str.fullmatch("[ACGTacgt]")
    ].copy()
    qtl["variant_key"] = _snv_key(qtl["position"], qtl["ref"], qtl["alt"])
    qtl = _deduplicate_aligned_variants(
        qtl,
        ["gene_id", "variant_key"],
        ["beta", "se"],
        dataset_id,
    )

    summaries: list[dict[str, Any]] = []
    variant_tables: list[pd.DataFrame] = []
    for gene_id, molecular_trait in qtl.groupby("gene_id", sort=True):
        aligned = gwas.merge(
            molecular_trait,
            on="variant_key",
            how="inner",
            suffixes=("_gwas", "_qtl"),
            validate="one_to_one",
        )
        if len(aligned) < 2:
            continue
        result = coloc_abf(
            aligned["beta_gwas"].to_numpy(),
            aligned["standard_error"].to_numpy(),
            aligned["beta_qtl"].to_numpy(),
            aligned["se"].to_numpy(),
            variant_ids=aligned["variant_key"].tolist(),
        )
        prior_sensitivity = {
            variance: coloc_abf(
                aligned["beta_gwas"].to_numpy(),
                aligned["standard_error"].to_numpy(),
                aligned["beta_qtl"].to_numpy(),
                aligned["se"].to_numpy(),
                gwas_prior_variance=variance,
                qtl_prior_variance=variance,
            ).pp_h4
            for variance in (0.01, 0.16)
        }
        top_index = result.top_shared_index
        top_rsid = aligned.iloc[top_index].get("rsid_gwas")
        if pd.isna(top_rsid) or not str(top_rsid):
            top_rsid = aligned.iloc[top_index].get("rsid_qtl")
        credible_indices = result.shared_credible_set(0.95)
        summaries.append(
            {
                "evidence_type": evidence_type,
                "dataset_id": dataset_id,
                "gene_id": gene_id,
                "gene_symbol": GENE_SYMBOLS.get(gene_id, gene_id),
                "shared_snv_count": len(aligned),
                "pp_h0": result.posterior["H0"],
                "pp_h1": result.posterior["H1"],
                "pp_h2": result.posterior["H2"],
                "pp_h3": result.posterior["H3"],
                "pp_h4": result.posterior["H4"],
                "pp_h4_prior_variance_0_01": prior_sensitivity[0.01],
                "pp_h4_prior_variance_0_16": prior_sensitivity[0.16],
                "pp_h4_sensitivity_min": min(
                    result.posterior["H4"], *prior_sensitivity.values()
                ),
                "pp_h4_sensitivity_max": max(
                    result.posterior["H4"], *prior_sensitivity.values()
                ),
                "top_shared_variant": top_rsid,
                "top_shared_variant_key": result.top_shared_variant,
                "top_shared_pip_given_h4": result.shared_variant_pip[top_index],
                "shared_credible_set_size_95": len(credible_indices),
            }
        )
        detail = aligned[
            [
                "variant_key",
                "rsid_gwas",
                "beta_gwas",
                "standard_error",
                "beta_qtl",
                "se",
            ]
        ].copy()
        detail.insert(0, "gene_symbol", GENE_SYMBOLS.get(gene_id, gene_id))
        detail.insert(0, "gene_id", gene_id)
        detail.insert(0, "dataset_id", dataset_id)
        detail.insert(0, "evidence_type", evidence_type)
        detail["gwas_pip_single_signal"] = result.gwas_pip
        detail["qtl_pip_single_signal"] = result.qtl_pip
        detail["shared_variant_pip_given_h4"] = result.shared_variant_pip
        detail["shared_variant_posterior"] = result.shared_variant_posterior
        detail["in_shared_95_credible_set"] = False
        detail.loc[credible_indices, "in_shared_95_credible_set"] = True
        variant_tables.append(detail)
    return summaries, variant_tables


def analyse_sort1_case(
    snapshot_directory: str | Path,
    output_directory: str | Path,
    prioritized_genes: pd.DataFrame,
) -> dict[str, Any]:
    """Fine-map the real 1 Mb GWAS slice and join GTEx liver summaries."""

    source = Path(snapshot_directory)
    output = Path(output_directory)
    gwas_path = source / "gwas_gcst90018961_chr1_108774968_109774968.tsv.gz"
    gtex_path = source / "gtex_v10_liver_rs12740374_eqtl.json"
    eqtl_path = source / "eqtl_catalogue_qtd000266_liver_sort1_genes.tsv.gz"
    pqtl_path = source / "eqtl_catalogue_qtd000584_plasma_sort1_genes.tsv.gz"
    source_manifest = source / "manifest.json"
    if not all(
        path.exists()
        for path in (gwas_path, gtex_path, eqtl_path, pqtl_path, source_manifest)
    ):
        raise FileNotFoundError("SORT1 case-study snapshots are incomplete; run scripts/fetch_sort1_case_study.py")

    manifest = json.loads(source_manifest.read_text(encoding="utf-8"))
    recorded_files = manifest.get("files", {})
    if not isinstance(recorded_files, dict) or not recorded_files:
        raise ValueError("SORT1 source manifest contains no file checksums")
    for filename, metadata in recorded_files.items():
        path = source / filename
        if not path.exists():
            raise FileNotFoundError(f"SORT1 source file listed in manifest is missing: {filename}")
        expected = metadata.get("sha256") if isinstance(metadata, dict) else None
        actual = sha256_file(path)
        if not expected or actual != expected:
            raise ValueError(
                f"SORT1 source checksum mismatch for {filename}: "
                f"expected {expected!r}, observed {actual}"
            )

    variants = pd.read_csv(gwas_path, sep="\t", dtype={"p_value": str})
    variants["log_approximate_bayes_factor"] = wakefield_log_abf(
        variants["beta"].to_numpy(), variants["standard_error"].to_numpy()
    )
    variants["approximate_pip"] = single_signal_pip(
        variants["beta"].to_numpy(), variants["standard_error"].to_numpy()
    )
    chosen = set(credible_set(variants["approximate_pip"].to_numpy(), 0.95).tolist())
    variants["in_95_credible_set"] = [index in chosen for index in range(len(variants))]
    variants = variants.sort_values("approximate_pip", ascending=False).reset_index(drop=True)
    variants["pip_rank"] = variants.index + 1
    variants.to_csv(output / "sort1_finemapping_all.csv.gz", index=False, compression="gzip")
    top_variants = variants.head(25).copy()
    top_variants.to_csv(output / "sort1_finemapping_top25.csv", index=False)

    gwas_snvs = variants.loc[
        variants["effect_allele"].astype(str).str.fullmatch("[ACGTacgt]")
        & variants["other_allele"].astype(str).str.fullmatch("[ACGTacgt]")
    ].copy()
    gwas_snvs["variant_key"] = _snv_key(
        gwas_snvs["base_pair_location"],
        gwas_snvs["effect_allele"],
        gwas_snvs["other_allele"],
    )
    gwas_snvs = _deduplicate_aligned_variants(
        gwas_snvs,
        ["variant_key"],
        ["beta", "standard_error"],
        "GWAS",
    )
    coloc_summaries: list[dict[str, Any]] = []
    coloc_variants: list[pd.DataFrame] = []
    for qtl_file, evidence_type, dataset_id in (
        (eqtl_path, "liver_eQTL", "QTD000266"),
        (pqtl_path, "plasma_pQTL", "QTD000584"),
    ):
        summaries, details = _regional_colocalization(
            gwas_snvs,
            qtl_file,
            evidence_type=evidence_type,
            dataset_id=dataset_id,
        )
        coloc_summaries.extend(summaries)
        coloc_variants.extend(details)
    colocalization = pd.DataFrame(coloc_summaries).sort_values(
        ["evidence_type", "pp_h4"], ascending=[True, False]
    )
    colocalization.to_csv(output / "sort1_colocalization.csv", index=False)
    pd.concat(coloc_variants, ignore_index=True).to_csv(
        output / "sort1_coloc_shared_variants.csv.gz",
        index=False,
        compression="gzip",
    )

    gtex = pd.DataFrame(json.loads(gtex_path.read_text(encoding="utf-8")))
    ranking = prioritized_genes.loc[
        prioritized_genes["studyLocusId"] == SORT1_STUDY_LOCUS_ID,
        [
            "gene_symbol",
            "rank",
            "posterior_mean",
            "posterior_lower",
            "posterior_upper",
            "conditional_posterior_mean",
            "l2g_score_external",
            "eQtlColocH4Maximum",
            "pQtlColocH4Maximum",
            "sQtlColocH4Maximum",
            "vepMean",
            "e2gMean",
        ],
    ]
    genes = gtex.merge(ranking, left_on="geneSymbol", right_on="gene_symbol", how="left")
    genes = genes.sort_values("rank")
    genes.to_csv(output / "sort1_gene_evidence.csv", index=False)

    return {
        "gwas_study_id": "GCST90018961",
        "main_ranking_study_id": "GCST90019512",
        "lead_variant": "rs12740374",
        "gwas_variant_count": int(len(variants)),
        "approximate_95_credible_set_size": int(len(chosen)),
        "approximate_95_credible_set_coverage": float(variants.loc[variants["in_95_credible_set"], "approximate_pip"].sum()),
        "lead_variant_pip": float(variants.iloc[0]["approximate_pip"]),
        "source_tls_verified": bool(manifest.get("tls_verified_for_json_and_metadata")),
        "verified_source_file_count": len(recorded_files),
        "regional_colocalization_count": len(colocalization),
        "maximum_eqtl_pp_h4": float(
            colocalization.loc[
                colocalization["evidence_type"] == "liver_eQTL", "pp_h4"
            ].max()
        ),
        "sort1_pqtl_pp_h4": float(
            colocalization.loc[
                (colocalization["evidence_type"] == "plasma_pQTL")
                & (colocalization["gene_symbol"] == "SORT1"),
                "pp_h4",
            ].iloc[0]
        ),
        "source_manifest": manifest,
        "top_variants": top_variants,
        "genes": genes,
        "colocalization": colocalization,
    }
