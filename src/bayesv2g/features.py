"""Transform Open Targets GraphQL records into analysis tables."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def _rsid(variant: dict[str, Any] | None) -> str:
    if not variant:
        return ""
    rsids = variant.get("rsIds") or []
    return rsids[0] if rsids else variant.get("id", "")


def _tss(location: dict[str, Any] | None) -> int | None:
    if not location:
        return None
    strand = str(location.get("strand", "")).upper()
    if strand in {"FORWARD", "POSITIVE", "1", "+", "+1"}:
        return int(location["start"])
    return int(location["end"])


def flatten_open_targets(
    payload: dict[str, Any],
    locus_labels: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return candidate-gene, locus, and credible-set-variant tables."""

    labels = locus_labels or {}
    credible_sets = payload.get("data", {}).get("credibleSets", {})
    rows = credible_sets.get("rows", [])
    if not rows:
        raise ValueError("Snapshot contains no credible-set rows")

    candidate_rows: list[dict[str, Any]] = []
    locus_rows: list[dict[str, Any]] = []
    variant_rows: list[dict[str, Any]] = []

    for locus in rows:
        locus_id = locus["studyLocusId"]
        lead = locus.get("variant") or {}
        locus_variants = locus.get("locus", {}).get("rows", [])
        coverage_95 = float(
            sum(v.get("posteriorProbability") or 0.0 for v in locus_variants if v.get("is95CredibleSet"))
        )
        coverage_99 = float(
            sum(v.get("posteriorProbability") or 0.0 for v in locus_variants if v.get("is99CredibleSet"))
        )
        study = locus.get("study") or {}
        locus_label = labels.get(locus_id, _rsid(lead))

        locus_rows.append(
            {
                "studyLocusId": locus_id,
                "locus_label": locus_label,
                "study_id": locus.get("studyId"),
                "trait": study.get("traitFromSource"),
                "chromosome": locus.get("chromosome"),
                "position": locus.get("position"),
                "lead_variant": _rsid(lead),
                "lead_variant_id": lead.get("id"),
                "p_value": (locus.get("pValueMantissa") or np.nan)
                * (10.0 ** (locus.get("pValueExponent") or 0)),
                "beta": locus.get("beta"),
                "standard_error": locus.get("standardError"),
                "sample_size": locus.get("sampleSize") or study.get("nSamples"),
                "finemapping_method": locus.get("finemappingMethod"),
                "credible_set_confidence": locus.get("confidence"),
                "credible_set_size_95": sum(bool(v.get("is95CredibleSet")) for v in locus_variants),
                "pip_coverage_95": coverage_95,
                "credible_set_size_99": sum(bool(v.get("is99CredibleSet")) for v in locus_variants),
                "pip_coverage_99": coverage_99,
                "candidate_count": locus.get("l2GPredictions", {}).get("count", 0),
                "publication_title": study.get("publicationTitle"),
                "pubmed_id": study.get("pubmedId"),
            }
        )

        for variant in locus_variants:
            annotation = variant.get("variant") or {}
            consequence = annotation.get("mostSevereConsequence") or {}
            variant_rows.append(
                {
                    "studyLocusId": locus_id,
                    "locus_label": locus_label,
                    "variant_id": annotation.get("id"),
                    "rsid": _rsid(annotation),
                    "chromosome": annotation.get("chromosome"),
                    "position": annotation.get("position"),
                    "reference_allele": annotation.get("referenceAllele"),
                    "alternate_allele": annotation.get("alternateAllele"),
                    "posterior_inclusion_probability": variant.get("posteriorProbability"),
                    "is_95_credible_set": variant.get("is95CredibleSet"),
                    "is_99_credible_set": variant.get("is99CredibleSet"),
                    "log_bayes_factor": variant.get("logBF"),
                    "most_severe_consequence": consequence.get("label"),
                }
            )

        for prediction in locus.get("l2GPredictions", {}).get("rows", []):
            target = prediction.get("target") or {}
            location = target.get("genomicLocation") or {}
            tss = _tss(location)
            lead_position = locus.get("position")
            feature_values = {item["name"]: item["value"] for item in prediction.get("features", [])}
            row = {
                "studyLocusId": locus_id,
                "locus_label": locus_label,
                "lead_variant": _rsid(lead),
                "chromosome": locus.get("chromosome"),
                "position": lead_position,
                "geneId": target.get("id"),
                "gene_symbol": target.get("approvedSymbol"),
                "gene_name": target.get("approvedName"),
                "gene_biotype": target.get("biotype"),
                "gene_tss": tss,
                "distance_to_tss_bp": abs(int(lead_position) - tss)
                if lead_position is not None and tss is not None
                else np.nan,
                "l2g_score_external": prediction.get("score"),
                "candidate_source": "Open Targets L2G surfaced candidates (score >= 0.05)",
            }
            row.update(feature_values)
            candidate_rows.append(row)

    candidates = pd.DataFrame(candidate_rows)
    loci = pd.DataFrame(locus_rows)
    variants = pd.DataFrame(variant_rows)
    order = {label: index for index, label in enumerate(labels)}
    if order:
        candidates["_order"] = candidates["studyLocusId"].map(order)
        loci["_order"] = loci["studyLocusId"].map(order)
        variants["_order"] = variants["studyLocusId"].map(order)
        candidates = candidates.sort_values(["_order", "gene_symbol"]).drop(columns="_order")
        loci = loci.sort_values("_order").drop(columns="_order")
        variants = variants.sort_values(["_order", "posterior_inclusion_probability"], ascending=[True, False]).drop(columns="_order")
    return candidates.reset_index(drop=True), loci.reset_index(drop=True), variants.reset_index(drop=True)


def evidence_coverage(candidates: pd.DataFrame) -> pd.Series:
    """Fraction of the five non-distance evidence classes with a non-zero value."""

    columns = [
        "eQtlColocH4Maximum",
        "pQtlColocH4Maximum",
        "sQtlColocH4Maximum",
        "vepMean",
        "e2gMean",
    ]
    available = [column for column in columns if column in candidates]
    if not available:
        return pd.Series(0.0, index=candidates.index)
    return candidates[available].fillna(0).gt(0).mean(axis=1)
