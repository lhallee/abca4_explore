"""Resumable AlphaGenome Atlas AVI scoring."""

from __future__ import annotations

import json
import math
import random
import time

from pathlib import Path
from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd

from alphagenome.atlas import atlas
from alphagenome.data import genome

from .config import Credentials, redact_text
from .links import build_variant_url, normalize_variant


AVI_FEATURES = (
    "MERGED_SPLICING",
    "MAX_ABS_ATAC",
    "MAX_ABS_CONTACT_MAPS",
    "MAX_ABS_DNASE",
    "MAX_ABS_CHIP_TF",
    "MAX_ABS_CHIP_HISTONE",
    "MAX_ABS_CAGE",
    "MAX_ABS_PROCAP",
    "MAX_ABS_RNA_SEQ",
    "MAX_ABS_POLYADENYLATION",
    "ALPHAMISSENSE",
    "CACTUS_241_WAY",
    "PROTEIN_TERMINATION",
    "START_LOST",
    "STOP_LOST",
    "PHASTCONS_470_WAY",
    "IS_INSERTION",
    "IS_DELETION",
)
REQUESTED_SCORERS = ("AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE")
TRACK_SCORERS = (
    "SPLICE_SITES",
    "SPLICE_SITE_USAGE",
    "SPLICE_JUNCTIONS",
    "ATAC",
    "CONTACT_MAPS",
    "DNASE",
    "CHIP_TF",
    "CHIP_HISTONE",
    "CAGE",
    "PROCAP",
    "RNA_SEQ",
    "POLYADENYLATION",
)


def _append_jsonl(path: Path, record: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as output_file:
        output_file.write(json.dumps(record, sort_keys=True) + "\n")
        output_file.flush()


def read_jsonl(path: Path) -> list[dict[str, object]]:
    """Read valid JSON records from a checkpoint file."""
    if not path.exists():
        return []
    records: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as input_file:
        for line in input_file:
            if line.strip():
                records.append(json.loads(line))
    return records


def score_variant(client: atlas.AtlasClient, variant: str) -> dict[str, object]:
    """Query one SNV and return calibrated AVI fields and 18 weights."""
    normalized = normalize_variant(variant)
    scores = client.query_variant(
        genome.Variant.from_str(normalized), requested_scorers=REQUESTED_SCORERS
    )
    avi = scores["AVI_SCORE"]
    feature_data = scores["AVI_SCORE_FEATURE_IMPORTANCE"]
    raw_score = float(np.ravel(avi.X)[0])
    cdf_quantile = float(np.ravel(avi.layers["quantiles"])[0])
    tail_probability = max(1.0 - cdf_quantile, np.finfo(float).tiny)
    phred = -10.0 * math.log10(tail_probability)

    feature_values = np.ravel(feature_data.X)  # (18,)
    if "name" in feature_data.var.columns:
        feature_names = tuple(str(name) for name in feature_data.var["name"])
    else:
        feature_names = AVI_FEATURES
    if len(feature_names) != len(feature_values):
        raise ValueError("AlphaGenome returned an unexpected feature count")
    weights = {
        f"fi_{name}": float(value)
        for name, value in zip(feature_names, feature_values, strict=True)
    }
    top_index = int(np.argmax(np.abs(feature_values)))
    return {
        "variant": normalized,
        "avi_raw": raw_score,
        "avi_cdf_quantile": cdf_quantile,
        "avi_tail_probability": tail_probability,
        "avi_phred": phred,
        "top_percentile": tail_probability * 100.0,
        "top_modality": feature_names[top_index],
        "top_feature_importance": float(feature_values[top_index]),
        "atlas_url": build_variant_url(normalized),
        **weights,
    }


def query_track_details(variants: Iterable[str], *, api_key: str) -> pd.DataFrame:
    """Return maximum absolute underlying track effects for selected variants."""
    client = atlas.create(api_key)
    rows: list[dict[str, object]] = []
    for variant in variants:
        normalized = normalize_variant(str(variant))
        scores = client.query_variant(
            genome.Variant.from_str(normalized),
            requested_scorers=TRACK_SCORERS,
        )
        for scorer_name, score_data in scores.items():
            matrix = np.asarray(score_data.X)
            if matrix.size == 0:
                continue
            flat_index = int(np.nanargmax(np.abs(matrix)))
            row_index, column_index = np.unravel_index(flat_index, matrix.shape)
            track = (
                score_data.var.iloc[column_index]
                if len(score_data.var) > column_index
                else pd.Series(dtype=object)
            )
            target = (
                score_data.obs.iloc[row_index]
                if len(score_data.obs) > row_index
                else pd.Series(dtype=object)
            )
            rows.append(
                {
                    "variant": normalized,
                    "scorer": scorer_name,
                    "raw_score": float(matrix[row_index, column_index]),
                    "track_index": int(column_index),
                    "track_name": str(track.get("name", "")),
                    "biosample_name": str(track.get("biosample_name", "")),
                    "ontology_curie": str(track.get("ontology_curie", "")),
                    "gene_name": str(target.get("gene_name", target.get("name", ""))),
                    "gene_id": str(target.get("gene_id", "")),
                    "atlas_url": build_variant_url(normalized),
                }
            )
    return pd.DataFrame(rows)


def score_variants(
    variants: Iterable[str],
    *,
    api_key: str,
    checkpoint_path: Path,
    failure_path: Path,
    requests_per_second: float = 5.0,
    resume: bool = True,
    max_retries: int = 4,
    credentials: Credentials | None = None,
) -> pd.DataFrame:
    """Score variants serially with per-record checkpoints."""
    if requests_per_second <= 0:
        raise ValueError("requests_per_second must be positive")
    existing = read_jsonl(checkpoint_path) if resume else []
    completed = {str(record["variant"]) for record in existing}
    normalized_variants = [normalize_variant(str(variant)) for variant in variants]
    if all(variant in completed for variant in normalized_variants):
        return pd.DataFrame(existing)
    client = atlas.create(api_key)
    minimum_interval = 1.0 / requests_per_second
    last_request = 0.0

    for normalized in normalized_variants:
        if normalized in completed:
            continue
        for attempt in range(max_retries + 1):
            elapsed = time.monotonic() - last_request
            if elapsed < minimum_interval:
                time.sleep(minimum_interval - elapsed)
            last_request = time.monotonic()
            try:
                record = score_variant(client, normalized)
                _append_jsonl(checkpoint_path, record)
                existing.append(record)
                completed.add(normalized)
                break
            except Exception as exc:
                if attempt == max_retries:
                    _append_jsonl(
                        failure_path,
                        {
                            "variant": normalized,
                            "error_type": type(exc).__name__,
                            "error": redact_text(
                                str(exc), credentials or Credentials()
                            ),
                        },
                    )
                    break
                time.sleep((2**attempt) + random.uniform(0.0, 0.5))
    return pd.DataFrame(existing)


def annotate_table(
    input_path: Path,
    output_path: Path,
    *,
    api_key: str,
    checkpoint_path: Path,
    failure_path: Path,
    requests_per_second: float = 5.0,
    resume: bool = True,
    max_variants: int | None = None,
) -> pd.DataFrame:
    """Score a tabular callset and preserve every input column."""
    separator = "\t" if input_path.suffix.lower() in {".tsv", ".txt"} else ","
    if input_path.suffix.lower() == ".parquet":
        variants = pd.read_parquet(input_path)
    else:
        variants = pd.read_csv(input_path, sep=separator)
    if "variant" not in variants.columns:
        required = {"chromosome", "position", "ref", "alt"}
        missing = required - set(variants.columns)
        if missing:
            raise ValueError(f"Input is missing columns: {sorted(missing)}")
        variants["variant"] = [
            f"chr{str(chrom).removeprefix('chr')}:{int(pos)}:{ref}>{alt}"
            for chrom, pos, ref, alt in zip(
                variants["chromosome"],
                variants["position"],
                variants["ref"],
                variants["alt"],
                strict=True,
            )
        ]
    target_variants = variants["variant"].dropna().astype(str).tolist()
    if max_variants is not None:
        target_variants = target_variants[:max_variants]
    scores = score_variants(
        target_variants,
        api_key=api_key,
        checkpoint_path=checkpoint_path,
        failure_path=failure_path,
        requests_per_second=requests_per_second,
        resume=resume,
        credentials=Credentials(alphagenome_api_key=api_key),
    )
    annotated = variants.merge(scores, on="variant", how="inner")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.suffix.lower() == ".parquet":
        annotated.to_parquet(output_path, index=False)
    else:
        output_separator = "\t" if output_path.suffix.lower() == ".tsv" else ","
        annotated.to_csv(output_path, sep=output_separator, index=False)
    return annotated
