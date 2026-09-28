"""Versioned access to the public data used by BayesV2G."""

from __future__ import annotations

import hashlib
import json
import shutil
import ssl
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


OPEN_TARGETS_QUERY = r"""
query BayesV2GLoci($ids: [String!]) {
  meta { dataVersion { year month iteration } }
  credibleSets(studyLocusIds: $ids, page: {index: 0, size: 20}) {
    count
    rows {
      studyLocusId studyId chromosome position locusStart locusEnd
      pValueMantissa pValueExponent beta standardError sampleSize
      finemappingMethod confidence credibleSetIndex qualityControls
      variant {
        id rsIds chromosome position referenceAllele alternateAllele
        mostSevereConsequence { id label }
      }
      study {
        id traitFromSource studyType projectId nSamples hasSumstats
        summarystatsLocation publicationTitle publicationFirstAuthor
        publicationDate pubmedId analysisFlags
        ldPopulationStructure { ldPopulation relativeSampleSize }
        diseases { id name }
      }
      locus(page: {index: 0, size: 150}) {
        count
        rows {
          posteriorProbability is95CredibleSet is99CredibleSet logBF
          beta standardError pValueMantissa pValueExponent
          variant {
            id rsIds chromosome position referenceAllele alternateAllele
            mostSevereConsequence { id label }
          }
        }
      }
      l2GPredictions(page: {index: 0, size: 50}) {
        count
        rows {
          score shapBaseValue
          target {
            id approvedSymbol approvedName biotype
            genomicLocation { chromosome start end strand }
          }
          features { name value shapValue }
        }
      }
    }
  }
}
"""


def graphql_request(
    endpoint: str,
    query: str,
    variables: dict[str, Any],
    *,
    verify_tls: bool = True,
    timeout: int = 120,
) -> dict[str, Any]:
    """Issue one GraphQL request and return validated JSON.

    TLS verification is on by default. ``verify_tls=False`` exists only for
    environments with a locally intercepted certificate chain and is recorded
    in every generated snapshot.
    """

    body = json.dumps({"query": query, "variables": variables}).encode("utf-8")
    request = urllib.request.Request(
        endpoint,
        data=body,
        headers={"Content-Type": "application/json", "User-Agent": "BayesV2G/0.1"},
        method="POST",
    )
    context = None if verify_tls else ssl._create_unverified_context()  # noqa: SLF001
    try:
        with urllib.request.urlopen(request, context=context, timeout=timeout) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as error:
        details = error.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"Open Targets returned HTTP {error.code}: {details}") from error
    except urllib.error.URLError as error:
        reason = error.reason
        hint = ""
        if isinstance(reason, ssl.SSLCertVerificationError):
            hint = (
                " Configure the correct CA bundle (for example with SSL_CERT_FILE). "
                "Use --insecure only on a trusted intercepted network; that choice is recorded."
            )
        raise RuntimeError(f"Could not reach Open Targets: {reason}.{hint}") from error
    except TimeoutError as error:
        raise RuntimeError(f"Open Targets request timed out after {timeout} seconds") from error
    except json.JSONDecodeError as error:
        raise RuntimeError("Open Targets returned a response that was not valid JSON") from error
    if payload.get("errors"):
        raise RuntimeError(f"Open Targets GraphQL errors: {payload['errors']}")
    return payload


def fetch_open_targets_loci(
    endpoint: str,
    study_locus_ids: list[str],
    *,
    verify_tls: bool = True,
) -> dict[str, Any]:
    payload = graphql_request(
        endpoint,
        OPEN_TARGETS_QUERY,
        {"ids": study_locus_ids},
        verify_tls=verify_tls,
    )
    validate_open_targets_snapshot(payload, study_locus_ids)
    rows = payload.get("data", {}).get("credibleSets", {}).get("rows", [])
    returned = {row["studyLocusId"] for row in rows}
    missing = set(study_locus_ids).difference(returned)
    if missing:
        raise ValueError(f"Open Targets did not return configured loci: {sorted(missing)}")
    payload["bayesv2g_snapshot"] = {
        "fetched_at_utc": datetime.now(timezone.utc).isoformat(),
        "endpoint": endpoint,
        "tls_verified": verify_tls,
        "requested_study_locus_ids": study_locus_ids,
    }
    return payload


def validate_open_targets_snapshot(
    payload: dict[str, Any],
    expected_study_locus_ids: list[str] | None = None,
    expected_release: str | None = None,
) -> None:
    """Reject stale, mismatched, or silently truncated GraphQL snapshots."""

    credible_sets = payload.get("data", {}).get("credibleSets")
    if not isinstance(credible_sets, dict):
        raise ValueError("Open Targets snapshot has no credibleSets object")
    rows = credible_sets.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError("Open Targets snapshot contains no credible-set rows")
    count = credible_sets.get("count")
    if count is not None and int(count) != len(rows):
        raise ValueError(
            "Open Targets credible-set response is truncated: "
            f"count={count}, rows={len(rows)}"
        )

    returned: list[str] = []
    for row in rows:
        locus_id = row.get("studyLocusId")
        if not isinstance(locus_id, str):
            raise ValueError("Open Targets credible-set row has no studyLocusId")
        returned.append(locus_id)
        for field in ("locus", "l2GPredictions"):
            page = row.get(field)
            if not isinstance(page, dict) or not isinstance(page.get("rows"), list):
                raise ValueError(f"Open Targets locus {locus_id} has no {field} rows")
            page_count = page.get("count")
            if page_count is not None and int(page_count) != len(page["rows"]):
                raise ValueError(
                    f"Open Targets {field} response is truncated for {locus_id}: "
                    f"count={page_count}, rows={len(page['rows'])}"
                )

    if len(returned) != len(set(returned)):
        raise ValueError("Open Targets snapshot contains duplicate studyLocusId rows")
    if expected_study_locus_ids is not None and set(returned) != set(expected_study_locus_ids):
        missing = sorted(set(expected_study_locus_ids).difference(returned))
        extra = sorted(set(returned).difference(expected_study_locus_ids))
        raise ValueError(
            f"Open Targets snapshot/config mismatch; missing={missing}, extra={extra}"
        )

    if expected_release is not None:
        version = payload.get("data", {}).get("meta", {}).get("dataVersion", {})
        actual = f"{version.get('year')}.{version.get('month')}"
        if actual != expected_release:
            raise ValueError(
                f"Open Targets snapshot release {actual!r} does not match configured "
                f"release {expected_release!r}"
            )


def write_json(payload: dict[str, Any], path: str | Path) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".part")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(output)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, destination: str | Path, *, verify_tls: bool = True) -> str:
    """Download a public artifact atomically and return its SHA-256 digest."""

    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".part")
    context = None if verify_tls else ssl._create_unverified_context()  # noqa: SLF001
    request = urllib.request.Request(url, headers={"User-Agent": "BayesV2G/0.1"})
    with urllib.request.urlopen(request, context=context, timeout=300) as response:
        with temporary.open("wb") as handle:
            shutil.copyfileobj(response, handle)
    temporary.replace(target)
    return sha256_file(target)
