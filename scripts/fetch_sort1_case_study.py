#!/usr/bin/env python3
"""Fetch a compact real-data SORT1 locus case study from official sources.

The GWAS file is queried remotely through its Tabix index, so only a 1 Mb
region is transferred from the approximately 411 MB genome-wide file.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import ssl
import tempfile
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pysam


GWAS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
    "GCST90018001-GCST90019000/GCST90018961/harmonised/GCST90018961.h.tsv.gz"
)
GWAS_META_URL = GWAS_URL + "-meta.yaml"
GTEX_ROOT = "https://gtexportal.org/api/v2"
ENSEMBL_ROOT = "https://rest.ensembl.org"
EQTL_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
    "QTS000015/QTD000266/QTD000266.all.tsv.gz"
)
PQTL_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
    "QTS000035/QTD000584/QTD000584.all.tsv.gz"
)
LEAD_POSITION = 109_274_968
REGION_START = LEAD_POSITION - 500_000
REGION_END = LEAD_POSITION + 500_000
GWAS_COLUMNS = [
    "chromosome",
    "base_pair_location",
    "effect_allele",
    "other_allele",
    "beta",
    "standard_error",
    "effect_allele_frequency",
    "p_value",
    "variant_id",
    "hm_coordinate_conversion",
    "hm_code",
    "rsid",
]
QTL_COLUMNS = [
    "molecular_trait_id",
    "chromosome",
    "position",
    "ref",
    "alt",
    "variant",
    "ma_samples",
    "maf",
    "pvalue",
    "beta",
    "se",
    "type",
    "ac",
    "an",
    "r2",
    "molecular_trait_object_id",
    "gene_id",
    "median_tpm",
    "rsid",
]
GENES = {
    "SORT1": "ENSG00000134243",
    "CELSR2": "ENSG00000143126",
    "PSRC1": "ENSG00000134222",
}


def _context(verify_tls: bool) -> ssl.SSLContext | None:
    return None if verify_tls else ssl._create_unverified_context()  # noqa: SLF001


def _get(url: str, verify_tls: bool, accept: str = "application/json") -> bytes:
    request = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "BayesV2G/0.1"})
    with urllib.request.urlopen(request, context=_context(verify_tls), timeout=120) as response:
        return response.read()


def _get_json(url: str, verify_tls: bool) -> Any:
    return json.loads(_get(url, verify_tls).decode("utf-8"))


def _json_url(root: str, path: str, parameters: dict[str, str]) -> str:
    return f"{root}{path}?{urllib.parse.urlencode(parameters)}"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _remote_tabix_rows(url: str, *, genes: set[str] | None = None) -> list[str]:
    """Fetch the configured region without leaving downloaded index files."""

    previous = Path.cwd()
    with tempfile.TemporaryDirectory(prefix="bayesv2g-tabix-") as temporary:
        try:
            os.chdir(temporary)
            tabix = pysam.TabixFile(url)
            rows = []
            for row in tabix.fetch("1", REGION_START - 1, REGION_END):
                if genes is not None:
                    fields = row.split("\t")
                    if len(fields) != len(QTL_COLUMNS) or fields[16] not in genes:
                        continue
                rows.append(row)
            tabix.close()
        finally:
            os.chdir(previous)
    return rows


def _write_gzip_table(path: Path, columns: list[str], rows: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write(row + "\n")


def fetch(output: Path, *, verify_tls: bool) -> None:
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    rows = _remote_tabix_rows(GWAS_URL)
    gwas_path = output / "gwas_gcst90018961_chr1_108774968_109774968.tsv.gz"
    _write_gzip_table(gwas_path, GWAS_COLUMNS, rows)

    eqtl_rows = _remote_tabix_rows(EQTL_URL, genes=set(GENES.values()))
    eqtl_path = output / "eqtl_catalogue_qtd000266_liver_sort1_genes.tsv.gz"
    _write_gzip_table(eqtl_path, QTL_COLUMNS, eqtl_rows)

    pqtl_rows = _remote_tabix_rows(PQTL_URL, genes=set(GENES.values()))
    pqtl_path = output / "eqtl_catalogue_qtd000584_plasma_sort1_genes.tsv.gz"
    _write_gzip_table(pqtl_path, QTL_COLUMNS, pqtl_rows)

    metadata_path = output / "gwas_gcst90018961_metadata.yaml"
    metadata_path.write_bytes(_get(GWAS_META_URL, verify_tls, "text/yaml"))

    variant_lookup = _get_json(
        _json_url(
            GTEX_ROOT,
            "/dataset/variant",
            {"snpId": "rs12740374", "datasetId": "gtex_v10"},
        ),
        verify_tls,
    )
    variant_id = variant_lookup["data"][0]["variantId"]
    gtex_records: list[dict[str, Any]] = []
    for symbol, ensembl_id in GENES.items():
        gene_lookup = _get_json(
            _json_url(
                GTEX_ROOT,
                "/reference/gene",
                {
                    "geneId": ensembl_id,
                    "gencodeVersion": "v39",
                    "genomeBuild": "GRCh38/hg38",
                },
            ),
            verify_tls,
        )
        gene = gene_lookup["data"][0]
        association = _get_json(
            _json_url(
                GTEX_ROOT,
                "/association/dyneqtl",
                {
                    "tissueSiteDetailId": "Liver",
                    "gencodeId": gene["gencodeId"],
                    "variantId": variant_id,
                    "datasetId": "gtex_v10",
                },
            ),
            verify_tls,
        )
        # Deliberately discard donor-level expression/genotype arrays.
        association.pop("data", None)
        association.pop("genotypes", None)
        association["ensemblGeneId"] = ensembl_id
        association["tss"] = gene["tss"]
        association["distanceToLeadBp"] = abs(int(gene["tss"]) - LEAD_POSITION)
        association["expectedGeneSymbol"] = symbol
        gtex_records.append(association)
    gtex_path = output / "gtex_v10_liver_rs12740374_eqtl.json"
    gtex_path.write_text(json.dumps(gtex_records, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    ensembl_genes = _get_json(
        f"{ENSEMBL_ROOT}/overlap/region/human/1:{REGION_START}-{REGION_END}"
        "?feature=gene;biotype=protein_coding",
        verify_tls,
    )
    ensembl_vep = _get_json(
        f"{ENSEMBL_ROOT}/vep/human/id/rs12740374?canonical=1;pick=1",
        verify_tls,
    )
    # Ensembl can return colocated records sourced from redistribution-restricted
    # databases such as HGMD and COSMIC. They are not used by this analysis.
    for consequence in ensembl_vep:
        consequence.pop("colocated_variants", None)
    ensembl_ld = _get_json(
        f"{ENSEMBL_ROOT}/ld/human/rs12740374/1000GENOMES:phase_3:EUR?window_size=500;r2=0.8",
        verify_tls,
    )
    for filename, payload in (
        ("ensembl_protein_coding_genes.json", ensembl_genes),
        ("ensembl_vep_rs12740374.json", ensembl_vep),
        ("ensembl_1000g_eur_ld_rs12740374.json", ensembl_ld),
        ("gtex_variant_lookup_rs12740374.json", variant_lookup),
    ):
        (output / filename).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    files = sorted(path for path in output.iterdir() if path.is_file() and path.name != "manifest.json")
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "tls_verified_for_json_and_metadata": verify_tls,
        "trait": "LDL cholesterol",
        "study_id": "GCST90018961",
        "lead_variant": "rs12740374",
        "lead_variant_grch38": "1_109274968_G_T",
        "genome_build": "GRCh38",
        "window": f"1:{REGION_START}-{REGION_END}",
        "gwas_rows": len(rows),
        "eqtl_rows": len(eqtl_rows),
        "pqtl_rows": len(pqtl_rows),
        "sources": {
            "gwas_summary_statistics": GWAS_URL,
            "gwas_metadata": GWAS_META_URL,
            "eqtl_catalogue_liver_qtd000266": EQTL_URL,
            "eqtl_catalogue_plasma_pqtl_qtd000584": PQTL_URL,
            "gtex_api": GTEX_ROOT,
            "ensembl_rest": ENSEMBL_ROOT,
        },
        "files": {path.name: {"bytes": path.stat().st_size, "sha256": _sha256(path)} for path in files},
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Wrote {len(files)} files and {len(rows):,} real GWAS rows to {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "snapshots" / "sort1_case",
    )
    parser.add_argument("--insecure", action="store_true", help="disable TLS verification for intercepted local networks")
    args = parser.parse_args()
    fetch(args.output, verify_tls=not args.insecure)


if __name__ == "__main__":
    main()
