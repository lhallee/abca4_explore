"""Rate-limited ClinVar retrieval and variation-centric XML parsing."""

from __future__ import annotations

import gzip
import json
import random
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from dataclasses import dataclass
from pathlib import Path
from collections.abc import Iterable, Iterator, Mapping, Sequence


EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
RETRYABLE_STATUS_CODES = frozenset({429, 500, 502, 503, 504})


@dataclass(frozen=True, kw_only=True)
class AlleleLocation:
    """A precise forward-strand VCF allele location."""

    chromosome: str
    position: int
    reference: str
    alternate: str
    assembly_accession: str | None

    @property
    def variant(self) -> str:
        chromosome = self.chromosome
        if not chromosome.lower().startswith("chr"):
            chromosome = f"chr{chromosome}"
        return f"{chromosome}:{self.position}:{self.reference}>{self.alternate}"


class RateLimiter:
    """Serial request limiter based on a minimum request interval."""

    def __init__(self, requests_per_second: float) -> None:
        if requests_per_second <= 0:
            raise ValueError("requests_per_second must be positive")
        self._minimum_interval = 1.0 / requests_per_second
        self._last_request = 0.0

    def wait(self) -> None:
        elapsed = time.monotonic() - self._last_request
        remaining = self._minimum_interval - elapsed
        if remaining > 0:
            time.sleep(remaining)
        self._last_request = time.monotonic()


class ClinVarClient:
    """Small E-utilities client with bounded retries and batched requests."""

    def __init__(
        self,
        *,
        api_key: str | None = None,
        email: str | None = None,
        requests_per_second: float | None = None,
        max_retries: int = 5,
    ) -> None:
        default_rate = 8.0 if api_key else 2.5
        self._api_key = api_key
        self._email = email
        self._limiter = RateLimiter(requests_per_second or default_rate)
        self._max_retries = max_retries

    def _request(self, endpoint: str, params: Mapping[str, object]) -> bytes:
        request_params = {key: str(value) for key, value in params.items()}
        request_params["tool"] = "abca4_alpha_genome"
        if self._api_key:
            request_params["api_key"] = self._api_key
        if self._email:
            request_params["email"] = self._email

        body = urllib.parse.urlencode(request_params).encode("utf-8")
        request = urllib.request.Request(
            urllib.parse.urljoin(EUTILS_URL, endpoint),
            data=body,
            headers={"User-Agent": "abca4-alpha-genome/0.1"},
            method="POST",
        )

        for attempt in range(self._max_retries + 1):
            self._limiter.wait()
            try:
                with urllib.request.urlopen(request, timeout=120) as response:
                    return response.read()
            except urllib.error.HTTPError as exc:
                if exc.code not in RETRYABLE_STATUS_CODES:
                    raise
                retry_after = exc.headers.get("Retry-After")
                delay = float(retry_after) if retry_after else 2**attempt
            except urllib.error.URLError:
                if attempt == self._max_retries:
                    raise
                delay = 2**attempt

            if attempt == self._max_retries:
                raise RuntimeError(f"ClinVar request failed after retries: {endpoint}")
            time.sleep(delay + random.uniform(0.0, 0.25))

        raise AssertionError("unreachable")

    def search(
        self,
        query: str,
        *,
        retmax: int = 0,
        page_size: int = 500,
    ) -> dict[str, object]:
        """Return all matching ClinVar Variation IDs unless retmax is set."""
        count_payload = self._request(
            "esearch.fcgi",
            {"db": "clinvar", "term": query, "retmode": "json", "retmax": 0},
        )
        total_count = int(json.loads(count_payload)["esearchresult"].get("count", 0))
        target_count = min(retmax, total_count) if retmax else total_count
        variant_ids: list[str] = []

        for start in range(0, target_count, page_size):
            size = min(page_size, target_count - start)
            payload = self._request(
                "esearch.fcgi",
                {
                    "db": "clinvar",
                    "term": query,
                    "retmode": "json",
                    "retstart": start,
                    "retmax": size,
                },
            )
            page = json.loads(payload)["esearchresult"].get("idlist", [])
            variant_ids.extend(str(identifier) for identifier in page)

        return {
            "query": query,
            "total_count": total_count,
            "fetched_count": len(variant_ids),
            "variant_ids": variant_ids,
        }

    def summary(self, variant_ids: Sequence[str]) -> list[dict[str, object]]:
        """Fetch ESummary records in bounded batches."""
        summaries: list[dict[str, object]] = []
        for batch in batched(variant_ids, 200):
            payload = self._request(
                "esummary.fcgi",
                {"db": "clinvar", "id": ",".join(batch), "retmode": "json"},
            )
            parsed = json.loads(payload).get("result", {})
            for identifier in parsed.get("uids", []):
                record = dict(parsed.get(identifier, {}))
                record["variant_id"] = str(identifier)
                summaries.append(record)
        return summaries

    def fetch_vcv_batch(self, variant_ids: Sequence[str]) -> bytes:
        """Fetch complete variation-centric XML for one ID batch."""
        return self._request(
            "efetch.fcgi",
            {
                "db": "clinvar",
                "id": ",".join(variant_ids),
                "rettype": "vcv",
                "retmode": "xml",
                "is_variationid": "true",
                "from_esearch": "true",
            },
        )


def batched(values: Sequence[str], size: int) -> Iterator[list[str]]:
    """Yield fixed-size list batches."""
    if size < 1:
        raise ValueError("batch size must be positive")
    for start in range(0, len(values), size):
        yield list(values[start : start + size])


def _first_text(root: ET.Element, paths: Iterable[str]) -> str | None:
    for path in paths:
        node = root.find(path)
        if node is not None and node.text:
            return node.text.strip()
    return None


def _review_stars(review_status: str | None) -> int:
    text = (review_status or "").lower()
    if "practice guideline" in text:
        return 4
    if "expert panel" in text:
        return 3
    if "multiple submitters" in text and "no conflicts" in text:
        return 2
    if "criteria provided" in text:
        return 1
    return 0


def extract_allele_location(archive: ET.Element) -> AlleleLocation | None:
    """Select an allele-level GRCh38 VCF location, never a gene interval."""
    alleles = archive.findall(".//ClassifiedRecord/SimpleAllele")
    if not alleles:
        alleles = archive.findall(".//SimpleAllele")
    unique_alleles = {id(allele): allele for allele in alleles}
    if len(unique_alleles) != 1:
        return None

    allele = next(iter(unique_alleles.values()))
    candidates = [
        location
        for location in allele.findall(".//Location/SequenceLocation")
        if location.attrib.get("Assembly") == "GRCh38"
        and location.attrib.get("positionVCF")
        and location.attrib.get("referenceAlleleVCF")
        and location.attrib.get("alternateAlleleVCF")
    ]
    if not candidates:
        return None
    candidates.sort(key=lambda location: location.attrib.get("forDisplay") != "true")
    location = candidates[0]
    position = int(location.attrib["positionVCF"])
    reference = location.attrib["referenceAlleleVCF"].upper()
    alternate = location.attrib["alternateAlleleVCF"].upper()
    start = int(location.attrib.get("start", position))
    stop = int(location.attrib.get("stop", position))
    if start != stop or len(reference) != 1 or len(alternate) != 1:
        return None

    return AlleleLocation(
        chromosome=location.attrib.get("Chr", ""),
        position=position,
        reference=reference,
        alternate=alternate,
        assembly_accession=location.attrib.get("Accession"),
    )


def _extract_conditions(archive: ET.Element) -> tuple[list[str], list[str]]:
    names: set[str] = set()
    identifiers: set[str] = set()
    for condition in archive.findall(".//ConditionList/Condition"):
        name = condition.attrib.get("Name") or _first_text(
            condition, ("./Name/ElementValue", ".//ElementValue")
        )
        if name:
            names.add(name)
        for xref in condition.findall(".//XRef"):
            database = xref.attrib.get("DB")
            identifier = xref.attrib.get("ID")
            if database and identifier:
                identifiers.add(f"{database}:{identifier}")
    for trait in archive.findall(".//TraitSet/Trait"):
        name = _first_text(
            trait,
            (
                './Name/ElementValue[@Type="Preferred"]',
                "./Name/ElementValue",
            ),
        )
        if name:
            names.add(name)
        for xref in trait.findall(".//XRef"):
            database = xref.attrib.get("DB")
            identifier = xref.attrib.get("ID")
            if database and identifier:
                identifiers.add(f"{database}:{identifier}")
    return sorted(names), sorted(identifiers)


def _extract_hgvs(archive: ET.Element) -> tuple[str | None, str | None]:
    coding: str | None = None
    protein: str | None = None
    for expression in archive.findall(".//SimpleAllele//HGVS/Expression"):
        if not expression.text:
            continue
        text = expression.text.strip()
        if coding is None and ":c." in text and "NM_000350" in text:
            coding = text
        if protein is None and ":p." in text:
            protein = text
    for change in archive.findall(".//SimpleAllele//ProteinChange"):
        if protein is None and change.text:
            protein = f"p.{change.text.strip()}"
    for name in archive.findall(".//SimpleAllele//Name"):
        if not name.text:
            continue
        text = name.text.strip()
        if protein is None and (":p." in text or text.startswith("p.")):
            protein = text
    variation_name = archive.attrib.get("VariationName", "")
    if coding is None and "NM_000350" in variation_name and ":c." in variation_name:
        coding = variation_name
    if protein is None and ":p." in variation_name:
        protein = variation_name
    return coding, protein


def _extract_rsid(archive: ET.Element) -> str | None:
    for xref in archive.findall(".//SimpleAllele/XRef"):
        if xref.attrib.get("DB") == "dbSNP" and xref.attrib.get("ID"):
            identifier = xref.attrib["ID"]
            return identifier if identifier.startswith("rs") else f"rs{identifier}"
    return None


def parse_vcv_xml(payload: bytes) -> list[dict[str, object]]:
    """Parse one EFetch payload into flat variation records."""
    root = ET.fromstring(payload)
    archives = (
        [root]
        if root.tag.endswith("VariationArchive")
        else list(root.findall(".//VariationArchive"))
    )
    records: list[dict[str, object]] = []

    for archive in archives:
        variation_id = archive.attrib.get("VariationID", "")
        significance = _first_text(
            archive,
            (
                ".//Classifications/GermlineClassification/Description",
                ".//GermlineClassification/Description",
                ".//ClinicalSignificance/Description",
            ),
        )
        review_status = _first_text(
            archive,
            (
                ".//Classifications/GermlineClassification/ReviewStatus",
                ".//GermlineClassification/ReviewStatus",
                ".//ClinicalSignificance/ReviewStatus",
            ),
        )
        classification_node = archive.find(".//Classifications/GermlineClassification")
        last_evaluated = (
            classification_node.attrib.get("DateLastEvaluated")
            if classification_node is not None
            else None
        )
        location = extract_allele_location(archive)
        conditions, condition_identifiers = _extract_conditions(archive)
        hgvs_c, hgvs_p = _extract_hgvs(archive)
        consequence_nodes = [
            *archive.findall(".//FunctionalConsequence"),
            *archive.findall(".//MolecularConsequence"),
        ]
        consequences = sorted(
            {
                node.attrib.get("Value") or node.attrib.get("Type") or ""
                for node in consequence_nodes
                if node.attrib.get("Value") or node.attrib.get("Type")
            }
        )
        citations = sorted(
            {
                node.text.strip()
                for node in archive.findall('.//Citation/ID[@Source="PubMed"]')
                if node.text
            }
        )
        has_haplotype = archive.find(".//Haplotype") is not None
        has_genotype = archive.find(".//Genotype") is not None
        records.append(
            {
                "clinvar_variation_id": variation_id,
                "vcv_accession": archive.attrib.get("Accession"),
                "variation_type": archive.attrib.get("VariationType"),
                "record_status": archive.attrib.get("RecordStatus"),
                "clinical_significance_raw": significance,
                "review_status": review_status,
                "review_stars": _review_stars(review_status),
                "last_evaluated": last_evaluated,
                "conditions": "|".join(conditions),
                "condition_identifiers": "|".join(condition_identifiers),
                "molecular_consequences": "|".join(consequences),
                "hgvs_c": hgvs_c,
                "hgvs_p": hgvs_p,
                "rsid": _extract_rsid(archive),
                "pmids": "|".join(citations),
                "has_haplotype": has_haplotype,
                "has_genotype": has_genotype,
                "chromosome": location.chromosome if location else None,
                "position": location.position if location else None,
                "ref": location.reference if location else None,
                "alt": location.alternate if location else None,
                "assembly_accession": (
                    location.assembly_accession if location else None
                ),
                "variant": location.variant if location else None,
                "coordinate_source": "ClinVar allele VCF" if location else None,
            }
        )
    return records


def save_vcv_batches(
    client: ClinVarClient,
    variant_ids: Sequence[str],
    output_dir: Path,
    *,
    batch_size: int = 200,
    resume: bool = True,
) -> list[Path]:
    """Fetch compressed XML batches and resume from existing files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for batch_index, identifier_batch in enumerate(
        batched(variant_ids, batch_size), start=1
    ):
        path = output_dir / f"vcv_{batch_index:04d}.xml.gz"
        paths.append(path)
        if resume and path.exists() and path.stat().st_size > 0:
            continue
        payload = client.fetch_vcv_batch(identifier_batch)
        temporary_path = path.with_suffix(path.suffix + ".part")
        with gzip.open(temporary_path, "wb") as output_file:
            output_file.write(payload)
        temporary_path.replace(path)
    return paths


def save_summary_batches(
    client: ClinVarClient,
    variant_ids: Sequence[str],
    output_dir: Path,
    *,
    batch_size: int = 200,
    resume: bool = True,
) -> list[Path]:
    """Fetch compressed ESummary batches with per-batch checkpoints."""
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for batch_index, identifier_batch in enumerate(
        batched(variant_ids, batch_size), start=1
    ):
        path = output_dir / f"summary_{batch_index:04d}.json.gz"
        paths.append(path)
        if resume and path.exists() and path.stat().st_size > 0:
            continue
        summaries = client.summary(identifier_batch)
        temporary_path = path.with_suffix(path.suffix + ".part")
        with gzip.open(temporary_path, "wt", encoding="utf-8") as output_file:
            json.dump(summaries, output_file, sort_keys=True)
        temporary_path.replace(path)
    return paths


def parse_saved_batches(paths: Sequence[Path]) -> list[dict[str, object]]:
    """Parse compressed XML batches in deterministic order."""
    records: list[dict[str, object]] = []
    for path in sorted(paths):
        with gzip.open(path, "rb") as input_file:
            records.extend(parse_vcv_xml(input_file.read()))
    return records
