"""Conservative dbSNP GRCh38 fallback resolution."""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.request

from collections.abc import Iterable
from pathlib import Path


REFSNP_URL = "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{rsid}"


def _is_grch38_placement(placement: dict[str, object]) -> bool:
    annotation = placement.get("placement_annot", {})
    if not isinstance(annotation, dict):
        return False
    traits = annotation.get("seq_id_traits_by_assembly", [])
    assembly_names = {
        str(item.get("assembly_name", "")) for item in traits if isinstance(item, dict)
    }
    return bool(annotation.get("is_chromosome")) and any(
        name.startswith("GRCh38") for name in assembly_names
    )


def parse_refsnp_grch38(
    payload: dict[str, object], *, preferred_alt: str | None = None
) -> dict[str, object] | None:
    """Extract one biallelic chr1 SNV placement from a RefSNP response."""
    primary = payload.get("primary_snapshot_data", {})
    if not isinstance(primary, dict):
        return None
    placements = primary.get("placements_with_allele", [])
    for placement in placements if isinstance(placements, list) else []:
        if not isinstance(placement, dict) or not _is_grch38_placement(placement):
            continue
        alleles = placement.get("alleles", [])
        parsed: list[tuple[str, int, str, str]] = []
        for allele in alleles if isinstance(alleles, list) else []:
            if not isinstance(allele, dict):
                continue
            spdi = allele.get("allele", {}).get("spdi", {})
            if not isinstance(spdi, dict):
                continue
            sequence = str(spdi.get("seq_id", ""))
            deleted = str(spdi.get("deleted_sequence", "")).upper()
            inserted = str(spdi.get("inserted_sequence", "")).upper()
            position = spdi.get("position")
            if (
                sequence.startswith("NC_000001.")
                and isinstance(position, int)
                and len(deleted) == 1
                and len(inserted) == 1
                and deleted in "ACGT"
                and inserted in "ACGT"
                and deleted != inserted
            ):
                parsed.append((sequence, position + 1, deleted, inserted))
        if preferred_alt:
            parsed = [item for item in parsed if item[3] == preferred_alt.upper()]
        if len(parsed) != 1:
            continue
        accession, position, reference, alternate = parsed[0]
        return {
            "chromosome": "1",
            "position": position,
            "ref": reference,
            "alt": alternate,
            "assembly_accession": accession,
            "variant": f"chr1:{position}:{reference}>{alternate}",
            "coordinate_source": "dbSNP GRCh38 chromosome placement",
        }
    return None


def resolve_rsid(
    rsid: str,
    *,
    preferred_alt: str | None = None,
    max_retries: int = 3,
) -> dict[str, object] | None:
    """Resolve one rsID using the public NCBI Variation API."""
    numeric = rsid.lower().removeprefix("rs")
    url = REFSNP_URL.format(rsid=numeric)
    for attempt in range(max_retries + 1):
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": "abca4-alpha-genome/0.1"}
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                return parse_refsnp_grch38(
                    json.load(response), preferred_alt=preferred_alt
                )
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            if exc.code not in {429, 500, 502, 503, 504}:
                raise
        except urllib.error.URLError:
            if attempt == max_retries:
                raise
        if attempt == max_retries:
            return None
        time.sleep(2**attempt)
    return None


def fill_missing_coordinates(
    records: Iterable[dict[str, object]],
    *,
    requests_per_second: float = 2.5,
    cache_path: Path | None = None,
) -> list[dict[str, object]]:
    """Resolve missing precise alleles by rsID and retain unresolved records."""
    output: list[dict[str, object]] = []
    cache: dict[str, dict[str, object] | None] = {}
    if cache_path is not None and cache_path.exists():
        loaded = json.loads(cache_path.read_text(encoding="utf-8"))
        if isinstance(loaded, dict):
            cache = loaded
    minimum_interval = 1.0 / requests_per_second
    last_request = 0.0
    for record in records:
        resolved = dict(record)
        rsid = resolved.get("rsid")
        if not resolved.get("variant") and rsid:
            key = str(rsid)
            if key not in cache:
                elapsed = time.monotonic() - last_request
                if elapsed < minimum_interval:
                    time.sleep(minimum_interval - elapsed)
                cache[key] = resolve_rsid(key)
                last_request = time.monotonic()
                if cache_path is not None:
                    cache_path.parent.mkdir(parents=True, exist_ok=True)
                    temporary_path = cache_path.with_suffix(cache_path.suffix + ".part")
                    temporary_path.write_text(
                        json.dumps(cache, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8",
                    )
                    temporary_path.replace(cache_path)
            placement = cache[key]
            if placement:
                resolved.update(placement)
        output.append(resolved)
    return output
