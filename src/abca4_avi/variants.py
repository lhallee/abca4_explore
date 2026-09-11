"""ClinVar label normalization, eligibility checks, and variant grouping."""

from __future__ import annotations

import re

from collections.abc import Iterable

import pandas as pd

from .gene_model import GeneModel, nearest_exon_distance


ABCA4_START = 93_992_834
ABCA4_END = 94_121_148
STGD_IDENTIFIERS = (
    "MedGen:C1855465",
    "MONDO:0009549",
    "OMIM:248200",
    "Orphanet:827",
)


def normalize_clinical_class(value: object) -> str:
    """Map ClinVar aggregate descriptions to stable analysis classes."""
    text = "" if pd.isna(value) else str(value)
    text = text.strip().lower().replace("_", " ")
    if not text or text in {"none", "unknown", "not provided"}:
        return "Other"
    if "conflict" in text:
        return "Conflicting"
    if "uncertain" in text or text == "vus":
        return "VUS"
    has_pathogenic = "pathogenic" in text
    has_benign = "benign" in text
    if has_pathogenic and has_benign:
        return "Conflicting"
    if has_pathogenic:
        return "P/LP"
    if has_benign:
        return "B/LB"
    return "Other"


def is_stargardt_record(conditions: object, identifiers: object) -> bool:
    """Identify records explicitly linked to STGD1."""
    condition_text = "" if pd.isna(conditions) else str(conditions).lower()
    identifier_text = "" if pd.isna(identifiers) else str(identifiers)
    aliases = ("stargardt", "fundus flavimaculatus", "stgd1")
    return any(alias in condition_text for alias in aliases) or any(
        identifier in identifier_text for identifier in STGD_IDENTIFIERS
    )


def _exclusion_reason(row: pd.Series) -> str | None:
    has_haplotype = row.get("has_haplotype")
    has_genotype = row.get("has_genotype")
    if pd.notna(has_haplotype) and bool(has_haplotype):
        return "haplotype"
    if pd.notna(has_genotype) and bool(has_genotype):
        return "genotype"
    variant = row.get("variant")
    if pd.isna(variant) or not str(variant).strip():
        raw_type = row.get("variation_type")
        variation_type = "" if pd.isna(raw_type) else str(raw_type).lower()
        if "deletion" in variation_type or "insertion" in variation_type:
            return "indel"
        if "copy number" in variation_type:
            return "copy_number_variant"
        return "missing_precise_grch38_allele"
    if str(row.get("chromosome")) not in {"1", "chr1"}:
        return "wrong_chromosome"
    position = int(row["position"])
    if not ABCA4_START <= position <= ABCA4_END:
        return "outside_abca4_locus"
    reference = str(row.get("ref") or "")
    alternate = str(row.get("alt") or "")
    if len(reference) != 1 or len(alternate) != 1:
        return "not_single_nucleotide_variant"
    if reference not in "ACGT" or alternate not in "ACGT":
        return "noncanonical_base"
    if reference == alternate:
        return "reference_equals_alternate"
    return None


def classify_variant_region(row: pd.Series, model: GeneModel | None) -> str:
    """Classify an SNV by consequence and splice-junction distance."""
    raw_consequences = row.get("molecular_consequences")
    consequences = "" if pd.isna(raw_consequences) else str(raw_consequences).lower()
    if "missense" in consequences:
        return "coding_missense"
    if pd.isna(row.get("position")) or model is None:
        return "other"
    position = int(row["position"])
    distance = nearest_exon_distance(position, model.exons)
    if distance == 0:
        return "other_exonic"
    if distance <= 2:
        return "canonical_splice"
    if distance <= 20:
        return "splice_region"
    if model.start <= position <= model.end:
        return "deep_intronic"
    return "other"


def _hgvs_offset(hgvs_c: object) -> int | None:
    text = "" if pd.isna(hgvs_c) else str(hgvs_c)
    match = re.search(r"c\.\d+([+-])(\d+)", text)
    if not match:
        return None
    return int(match.group(2))


def annotate_records(
    records: Iterable[dict[str, object]],
    model: GeneModel | None,
) -> pd.DataFrame:
    """Normalize flat ClinVar records and derive analysis fields."""
    frame = pd.DataFrame(records)
    if frame.empty:
        return frame
    frame["clinical_class"] = frame["clinical_significance_raw"].map(
        normalize_clinical_class
    )
    frame["binary_label"] = frame["clinical_class"].map({"B/LB": 0, "P/LP": 1})
    frame["is_stgd1"] = [
        is_stargardt_record(conditions, identifiers)
        for conditions, identifiers in zip(
            frame["conditions"], frame["condition_identifiers"], strict=True
        )
    ]
    frame["exclusion_reason"] = frame.apply(_exclusion_reason, axis=1)
    frame["eligible_for_avi"] = frame["exclusion_reason"].isna()
    frame["region_class"] = frame.apply(
        lambda row: classify_variant_region(row, model), axis=1
    )
    frame["hgvs_splice_offset"] = frame["hgvs_c"].map(_hgvs_offset)
    if model is not None:
        frame["nearest_exon_distance"] = frame["position"].map(
            lambda value: (
                nearest_exon_distance(int(value), model.exons)
                if pd.notna(value)
                else pd.NA
            )
        )
    else:
        frame["nearest_exon_distance"] = pd.NA
    return frame


def deduplicate_eligible_variants(frame: pd.DataFrame) -> pd.DataFrame:
    """Return one analysis row per precise genomic SNV."""
    eligible = frame[frame["eligible_for_avi"]].copy()
    if eligible.empty:
        return eligible
    eligible = eligible.sort_values(
        ["variant", "review_stars", "clinvar_variation_id"],
        ascending=[True, False, True],
    )
    rows: list[pd.Series] = []
    for _, group in eligible.groupby("variant", sort=True, dropna=False):
        row = group.iloc[0].copy()
        classes = set(group["clinical_class"].dropna()) - {"Other"}
        if "P/LP" in classes and "B/LB" in classes:
            row["clinical_class"] = "Conflicting"
            row["binary_label"] = pd.NA
        row["clinvar_variation_ids"] = "|".join(
            sorted(set(group["clinvar_variation_id"].astype(str)))
        )
        row["review_stars"] = int(group["review_stars"].max())
        row["is_stgd1"] = bool(group["is_stgd1"].any())
        rows.append(row)
    return pd.DataFrame(rows).reset_index(drop=True)
