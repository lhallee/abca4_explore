"""Population-frequency control candidates for the ABCA4 analyses."""

from __future__ import annotations

import gzip
import json
import urllib.request

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd


GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET = "gnomad_r4"
ABCA4_GENE_ID = "ENSG00000198691"
ABCA4_BA1_AF = 0.163
ABCA4_BS1_STRONG_AF = 0.0163
ABCA4_BS1_SUPPORTING_AF = 0.00163
GNOMAD_V4_GRPMAX_EXCLUDED_POPULATIONS = frozenset(
    {"", "ami", "asj", "fin", "oth", "remaining"}
)
ABCA4_BS1_EXCLUSIONS = frozenset(
    {
        "c.2588G>C",
        "c.3113C>T",
        "c.4253+43G>A",
        "c.4685T>C",
        "c.5603A>T",
        "c.5882G>A",
        "c.6320G>A",
    }
)

GNOMAD_GENE_QUERY = """
query Abca4PopulationControls(
  $geneId: String!
  $dataset: DatasetId!
  $referenceGenome: ReferenceGenomeId!
) {
  gene(gene_id: $geneId, reference_genome: $referenceGenome) {
    gene_id
    symbol
    variants(dataset: $dataset) {
      variant_id
      chrom
      pos
      ref
      alt
      rsids
      transcript_consequence {
        consequence_terms
        hgvsc
        hgvsp
        is_mane_select
        major_consequence
        transcript_id
      }
      joint {
        ac
        an
        homozygote_count
        populations {
          id
          ac
          an
          homozygote_count
        }
      }
    }
  }
}
"""


def fetch_abca4_gnomad(*, timeout_seconds: int = 120) -> dict[str, Any]:
    """Fetch the compact gnomAD v4 ABCA4 gene response."""
    payload = json.dumps(
        {
            "query": GNOMAD_GENE_QUERY,
            "variables": {
                "geneId": ABCA4_GENE_ID,
                "dataset": GNOMAD_DATASET,
                "referenceGenome": "GRCh38",
            },
        }
    ).encode("utf-8")
    request = urllib.request.Request(
        GNOMAD_API_URL,
        data=payload,
        headers={
            "Content-Type": "application/json",
            "User-Agent": "abca4-alpha-genome/0.1",
        },
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
        result = json.loads(response.read())

    errors = result.get("errors")
    if errors:
        messages = "; ".join(str(error.get("message", error)) for error in errors)
        raise RuntimeError(f"gnomAD GraphQL query failed: {messages}")
    gene = result.get("data", {}).get("gene")
    if not gene:
        raise RuntimeError(f"gnomAD returned no record for {ABCA4_GENE_ID}")
    return result


def save_raw_response(response: Mapping[str, object], path: Path) -> None:
    """Persist a compressed source response for provenance."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        json.dump(response, handle, separators=(",", ":"), sort_keys=True)
        handle.write("\n")


def _coding_hgvs(value: object) -> str | None:
    if not value:
        return None
    text = str(value)
    return text.rsplit(":", maxsplit=1)[-1]


def _max_population_frequency(
    populations: Sequence[Mapping[str, object]],
) -> tuple[float | None, str | None]:
    frequencies: list[tuple[float, str]] = []
    for population in populations:
        population_id = str(population.get("id") or "")
        if population_id in GNOMAD_V4_GRPMAX_EXCLUDED_POPULATIONS:
            continue
        allele_count = population.get("ac")
        allele_number = population.get("an")
        if allele_count is None or not allele_number:
            continue
        frequencies.append(
            (float(allele_count) / float(allele_number), population_id)
        )
    if not frequencies:
        return None, None
    return max(frequencies)


def frequency_evidence_tier(
    frequency: float | None,
    *,
    excluded_by_vcep: bool,
) -> str:
    """Return a candidate tier based on ABCA4 VCEP population thresholds."""
    if excluded_by_vcep:
        return "excluded_vcep_reduced_penetrance_or_hypomorphic"
    if frequency is None:
        return "insufficient_frequency_data"
    if frequency > ABCA4_BA1_AF:
        return "ba1_candidate"
    if frequency > ABCA4_BS1_STRONG_AF:
        return "bs1_strong_candidate"
    if frequency > ABCA4_BS1_SUPPORTING_AF:
        return "bs1_supporting_candidate"
    return "below_supporting_threshold"


def normalize_gnomad_missense(response: Mapping[str, object]) -> pd.DataFrame:
    """Normalize MANE/canonical ABCA4 missense variants from a gnomAD response."""
    gene = response.get("data", {}).get("gene", {})
    rows: list[dict[str, object]] = []
    for record in gene.get("variants", []):
        consequence = record.get("transcript_consequence") or {}
        terms = set(consequence.get("consequence_terms") or [])
        if "missense_variant" not in terms:
            continue
        reference = str(record.get("ref") or "")
        alternate = str(record.get("alt") or "")
        if len(reference) != 1 or len(alternate) != 1:
            continue

        coding_hgvs = _coding_hgvs(consequence.get("hgvsc"))
        excluded_by_vcep = coding_hgvs in ABCA4_BS1_EXCLUSIONS
        joint = record.get("joint") or {}
        maximum_af, maximum_population = _max_population_frequency(
            joint.get("populations") or []
        )
        chromosome = str(record.get("chrom") or "")
        if not chromosome.lower().startswith("chr"):
            chromosome = f"chr{chromosome}"
        rows.append(
            {
                "variant": (
                    f"{chromosome}:{int(record['pos'])}:{reference}>{alternate}"
                ),
                "gnomad_variant_id": record.get("variant_id"),
                "rsids": "|".join(record.get("rsids") or []),
                "hgvsc": consequence.get("hgvsc"),
                "hgvsp": consequence.get("hgvsp"),
                "transcript_id": consequence.get("transcript_id"),
                "is_mane_select": bool(consequence.get("is_mane_select")),
                "joint_ac": joint.get("ac"),
                "joint_an": joint.get("an"),
                "joint_homozygote_count": joint.get("homozygote_count"),
                "max_population_af": maximum_af,
                "max_population": maximum_population,
                "vcep_frequency_exclusion": excluded_by_vcep,
                "population_evidence_tier": frequency_evidence_tier(
                    maximum_af,
                    excluded_by_vcep=excluded_by_vcep,
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["max_population_af", "variant"], ascending=[False, True]
    )


def join_existing_screen(
    candidates: pd.DataFrame,
    mapping: pd.DataFrame,
) -> pd.DataFrame:
    """Attach ClinVar and Atlas-PPI mapping status without pseudoreplication."""
    columns = [
        "variant",
        "clinvar_variation_id",
        "clinical_class",
        "review_stars",
        "protein_status",
        "protein_effect",
        "normalized_hgvs_p",
        "sequence_id",
        "atlas_observable",
    ]
    available = [column for column in columns if column in mapping.columns]
    unique_mapping = mapping[available].drop_duplicates(subset=["variant"])
    joined = candidates.merge(unique_mapping, on="variant", how="left")
    joined["already_in_clinvar_snapshot"] = joined["clinvar_variation_id"].notna()
    joined["already_screened_by_atlas_ppi"] = joined["sequence_id"].notna()
    joined["candidate_label_scope"] = "population-compatible; not clinical benign"
    return joined
