"""End-to-end orchestration for the ABCA4 AlphaGenome analysis."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import re

from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from .atlas import AVI_FEATURES, annotate_table, query_track_details
from .clinvar import (
    ClinVarClient,
    parse_saved_batches,
    save_summary_batches,
    save_vcv_batches,
)
from .config import ProjectPaths, load_credentials
from .dbsnp import fill_missing_coordinates
from .gene_model import load_abca4_gene_model
from .literature import prepare_functional_assays
from .plots import generate_plots
from .statistics import evaluate_all_cohorts, evaluate_functional_assays
from .variants import annotate_records, deduplicate_eligible_variants


CLINVAR_QUERY = "ABCA4[gene]"
SENTINELS = (
    "c.4539+2001G>A",
    "c.5461-10T>C",
    "c.4539+2028C>T",
    "c.769-784C>T",
)


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def _write_frame(frame: pd.DataFrame, stem: Path) -> tuple[Path, Path]:
    tsv_path = stem.with_suffix(".tsv")
    parquet_path = stem.with_suffix(".parquet")
    stem.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(tsv_path, sep="\t", index=False)
    frame.to_parquet(parquet_path, index=False)
    return tsv_path, parquet_path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as input_file:
        while chunk := input_file.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact_paths(paths: ProjectPaths) -> list[Path]:
    roots = (
        paths.raw,
        paths.processed,
        paths.reference,
        paths.tables,
        paths.plots,
        paths.logs,
    )
    artifacts = {
        path
        for root in roots
        if root.exists()
        for path in root.rglob("*")
        if path.is_file()
        and not path.name.endswith(".part")
        and path.name != "alphagenome_atlas.pdf"
    }
    for report_name in ("summary.md", "impressiveness_assessment.md"):
        report_path = paths.root / "results" / report_name
        if report_path.exists():
            artifacts.add(report_path)
    return sorted(artifacts)


def _sentinel_name(hgvs_c: object) -> str | None:
    text = "" if pd.isna(hgvs_c) else str(hgvs_c)
    return next((name for name in SENTINELS if name in text), None)


def retrieve_clinvar(
    paths: ProjectPaths,
    *,
    resume: bool,
    max_records: int | None = None,
) -> pd.DataFrame:
    """Retrieve and normalize ABCA4 ClinVar variation records."""
    credentials = load_credentials(paths.root)
    client = ClinVarClient(
        api_key=credentials.ncbi_api_key,
        email=credentials.ncbi_email,
    )
    search_path = paths.raw / "clinvar_ids.json"
    if resume and search_path.exists():
        search_result = json.loads(search_path.read_text(encoding="utf-8"))
    else:
        search_result = client.search(CLINVAR_QUERY, retmax=max_records or 0)
        _write_json(search_path, search_result)
    identifiers = [str(value) for value in search_result["variant_ids"]]
    if max_records is not None:
        identifiers = identifiers[:max_records]

    save_summary_batches(
        client,
        identifiers,
        paths.raw / "clinvar_summary_batches",
        batch_size=200,
        resume=resume,
    )

    batch_paths = save_vcv_batches(
        client,
        identifiers,
        paths.raw / "clinvar_vcv_batches",
        batch_size=200,
        resume=resume,
    )
    records = parse_saved_batches(batch_paths)
    records = fill_missing_coordinates(
        records,
        cache_path=paths.cache / "dbsnp_resolutions.json",
    )
    model = load_abca4_gene_model(paths.cache)
    frame = annotate_records(records, model)
    frame["sentinel_name"] = frame["hgvs_c"].map(_sentinel_name)
    _write_frame(frame, paths.processed / "clinvar_records")
    exclusions = frame[~frame["eligible_for_avi"]].copy()
    _write_frame(exclusions, paths.tables / "clinvar_exclusions")
    variants = deduplicate_eligible_variants(frame)
    variants["sentinel_name"] = variants["hgvs_c"].map(_sentinel_name)
    _write_frame(variants, paths.processed / "variants")
    return variants


def score_callset(
    paths: ProjectPaths,
    variants: pd.DataFrame,
    *,
    resume: bool,
    max_variants: int | None = None,
) -> pd.DataFrame:
    """Score a deduplicated callset and write audit-ready tables."""
    credentials = load_credentials(paths.root)
    if not credentials.alphagenome_api_key:
        raise RuntimeError("ALPHAGENOME_API_KEY is unavailable")
    scoring_input = variants.copy()
    scoring_input["_sentinel_priority"] = scoring_input["sentinel_name"].notna()
    scoring_input = scoring_input.sort_values(
        ["_sentinel_priority", "variant"], ascending=[False, True]
    ).drop(columns="_sentinel_priority")
    input_path = paths.processed / "variants_for_scoring.tsv"
    scoring_input.to_csv(input_path, sep="\t", index=False)
    scored = annotate_table(
        input_path,
        paths.tables / "variants_scored.tsv",
        api_key=credentials.alphagenome_api_key,
        checkpoint_path=paths.cache / "alphagenome_scores.jsonl",
        failure_path=paths.logs / "alphagenome_failures.jsonl",
        requests_per_second=5.0,
        resume=resume,
        max_variants=max_variants,
    )
    _write_frame(scored, paths.tables / "variants_scored")
    sentinels = scored[scored["sentinel_name"].notna()].copy()
    _write_frame(sentinels, paths.tables / "sentinel_variants")

    detail_variants = list(sentinels["variant"].astype(str))
    detail_variants.extend(
        scored.nlargest(20, "avi_phred")["variant"].astype(str).tolist()
    )
    unique_detail_variants = list(dict.fromkeys(detail_variants))
    detail_path = paths.tables / "selected_track_attributions.parquet"
    cached_details = (
        pd.read_parquet(detail_path)
        if resume and detail_path.exists()
        else pd.DataFrame()
    )
    cached_variants = (
        set(cached_details["variant"].astype(str))
        if "variant" in cached_details
        else set()
    )
    if set(unique_detail_variants).issubset(cached_variants):
        details = cached_details
    else:
        details = query_track_details(
            unique_detail_variants, api_key=credentials.alphagenome_api_key
        )
        _write_frame(details, paths.tables / "selected_track_attributions")
    return scored


def _load_functional_assays(paths: ProjectPaths, scored: pd.DataFrame) -> pd.DataFrame:
    assays = prepare_functional_assays(paths.reference, paths.raw / "literature")
    if assays.empty or "variant" not in assays:
        return assays
    scored = scored.copy()
    scored["coding_key"] = scored["hgvs_c"].map(_coding_key)
    scored["protein_key"] = scored["hgvs_p"].map(_protein_key)
    coding_map = (
        scored.dropna(subset=["coding_key"])
        .drop_duplicates("coding_key")
        .set_index("coding_key")["variant"]
    )
    protein_map = (
        scored.dropna(subset=["protein_key"])
        .drop_duplicates("protein_key")
        .set_index("protein_key")["variant"]
    )
    assays["coding_key"] = assays["hgvs_c"].map(_coding_key)
    assays["protein_key"] = assays["hgvs_p"].map(_protein_key)
    assays["variant"] = assays["variant"].fillna(assays["coding_key"].map(coding_map))
    assays["variant"] = assays["variant"].fillna(assays["protein_key"].map(protein_map))
    return assays.merge(
        scored.drop(columns=["coding_key", "protein_key"]),
        on="variant",
        how="left",
        suffixes=("", "_clinvar"),
    )


def _coding_key(value: object) -> str | None:
    text = "" if pd.isna(value) else str(value)
    match = re.search(r"c\.([^\s;\]]+)", text)
    return match.group(1).replace("−", "-") if match else None


def _protein_key(value: object) -> str | None:
    text = "" if pd.isna(value) else str(value)
    match = re.search(r"p\.\(?([A-Za-z]{1,3}\d+[A-Za-z]{1,3})\)?", text)
    if not match:
        return None
    change = match.group(1)
    amino_acids = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Glu": "E",
        "Gln": "Q",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
    }
    for name, letter in amino_acids.items():
        change = change.replace(name, letter)
    return change.upper()


def create_report(
    paths: ProjectPaths,
    scored: pd.DataFrame,
    *,
    bootstrap_replicates: int = 2_000,
) -> dict[str, object]:
    """Compute metrics, figures, provenance, and acceptance checks."""
    all_clinvar = pd.read_parquet(paths.processed / "clinvar_records.parquet")
    sentinels = scored[scored["sentinel_name"].notna()].copy()
    functional = _load_functional_assays(paths, scored)
    if not functional.empty:
        _write_frame(functional, paths.tables / "functional_assays_scored")
    metrics = {
        "clinical": evaluate_all_cohorts(
            scored, bootstrap_replicates=bootstrap_replicates
        ),
        "functional": evaluate_functional_assays(functional),
    }
    metrics_path = paths.tables / "metrics.json"
    _write_json(metrics_path, metrics)
    generate_plots(all_clinvar, scored, sentinels, functional, paths.plots)

    feature_columns = [f"fi_{name}" for name in AVI_FEATURES]
    primary = scored[
        (scored["review_stars"] >= 1) & scored["clinical_class"].isin(["P/LP", "B/LB"])
    ]
    checks = {
        "both_primary_classes_present": set(primary["clinical_class"])
        == {"P/LP", "B/LB"},
        "all_18_attribution_columns_present": set(feature_columns).issubset(scored),
        "all_four_sentinels_scored": set(sentinels["sentinel_name"]) == set(SENTINELS),
        "every_score_has_atlas_url": bool(scored["atlas_url"].notna().all()),
        "no_failed_scores": not (paths.logs / "alphagenome_failures.jsonl").exists()
        or (paths.logs / "alphagenome_failures.jsonl").stat().st_size == 0,
    }
    _write_json(paths.tables / "acceptance_checks.json", checks)
    primary_metrics = metrics["clinical"]["primary_at_least_one_star"]
    summary_lines = [
        "# ABCA4 AlphaGenome AVI results",
        "",
        (
            "AlphaGenome outputs are research-only molecular predictions, "
            "not clinical diagnoses."
        ),
        "",
        f"Retrieved ClinVar records: {len(all_clinvar):,}.",
        f"Unique scored GRCh38 SNVs: {len(scored):,}.",
        (
            f"Primary cohort: {int(primary_metrics['n']):,} variants "
            f"({int(primary_metrics['n_positive']):,} P/LP, "
            f"{int(primary_metrics['n_negative']):,} B/LB)."
        ),
    ]
    if primary_metrics.get("status") == "ok":
        summary_lines.extend(
            [
                f"ROC-AUC: {primary_metrics['roc_auc']:.3f}.",
                f"Average precision: {primary_metrics['average_precision']:.3f}.",
                f"Trapezoidal PR-AUC: {primary_metrics['pr_auc_trapezoidal']:.3f}.",
            ]
        )

    summary_lines.extend(
        [
            "",
            "## Clinical cohort performance",
            "",
            "| Cohort | n | P/LP | B/LB | ROC-AUC | Average precision |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for name, cohort in metrics["clinical"].items():
        roc_auc = cohort.get("roc_auc")
        average_precision = cohort.get("average_precision")
        summary_lines.append(
            "| "
            f"{name.replace('_', ' ')} | {int(cohort['n']):,} | "
            f"{int(cohort['n_positive']):,} | {int(cohort['n_negative']):,} | "
            f"{roc_auc:.3f} | {average_precision:.3f} |"
            if roc_auc is not None and average_precision is not None
            else (
                "| "
                f"{name.replace('_', ' ')} | {int(cohort['n']):,} | "
                f"{int(cohort['n_positive']):,} | "
                f"{int(cohort['n_negative']):,} | NA | NA |"
            )
        )

    summary_lines.extend(
        [
            "",
            "## Sentinel variants",
            "",
            "| Variant | GRCh38 allele | AVI Phred | Top modality | Atlas |",
            "|---|---|---:|---|---|",
        ]
    )
    for row in sentinels.sort_values("sentinel_name").itertuples(index=False):
        summary_lines.append(
            f"| {row.sentinel_name} | {row.variant} | {row.avi_phred:.3f} | "
            f"{row.top_modality} | [open]({row.atlas_url}) |"
        )

    if not functional.empty:
        summary_lines.extend(
            [
                "",
                "## Published functional assay coverage",
                "",
                "Assay scales are analyzed within each study and are not pooled.",
                "",
                "| Study | Assay rows | Variants represented | Rows with AVI |",
                "|---|---:|---:|---:|",
            ]
        )
        for study_key, group in functional.groupby("study_key", sort=True):
            identity_column = "hgvs_p" if study_key == "garces_2021" else "hgvs_c"
            summary_lines.append(
                f"| {study_key} | {len(group):,} | "
                f"{group[identity_column].nunique():,} | "
                f"{group['avi_phred'].notna().sum():,} |"
            )
    summary_lines.extend(
        [
            "",
            "## Acceptance checks",
            "",
            *(f"- {name}: {value}" for name, value in checks.items()),
            "",
            (
                "See `metrics.json` for fixed-seed 95% bootstrap intervals "
                "and threshold metrics."
            ),
        ]
    )
    summary_path = paths.root / "results" / "summary.md"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    artifacts = _artifact_paths(paths)
    manifest = {
        "created_at_utc": datetime.now(UTC).isoformat(),
        "query": CLINVAR_QUERY,
        "cohort_definitions": {
            "primary": (
                "P/LP versus B/LB, precise GRCh38 SNVs, at least one "
                "ClinVar review star"
            ),
            "all_reviews": "P/LP versus B/LB, precise GRCh38 SNVs, any review status",
            "stgd1": (
                "Primary labels restricted to explicit STGD1 identifiers or synonyms"
            ),
        },
        "sentinels": list(SENTINELS),
        "record_counts": {
            "clinvar_records": len(all_clinvar),
            "scored_unique_snvs": len(scored),
            "primary": len(primary),
        },
        "python": platform.python_version(),
        "packages": {
            name: importlib.metadata.version(name)
            for name in (
                "abca4-avi",
                "alphagenome",
                "numpy",
                "pandas",
                "pyarrow",
                "scikit-learn",
                "scipy",
            )
        },
        "sources": [
            "https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/",
            "https://deepmind.google/science/alphagenome/",
            "https://doi.org/10.3390/ijms22010185",
            "https://doi.org/10.1101/gr.226621.117",
            "https://doi.org/10.1167/iovs.65.10.2",
            "https://doi.org/10.1167/iovs.66.1.65",
        ],
        "limitations": [
            "AlphaGenome is for research use only.",
            (
                "AVI predicts molecular impact and does not establish disease "
                "causality, penetrance, or prognosis."
            ),
            (
                "ClinVar aggregate classifications are categorical clinical "
                "evidence, not functional assay measurements."
            ),
        ],
        "acceptance_checks": checks,
        "files": {
            str(path.relative_to(paths.root)): {
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
                "last_modified_utc": datetime.fromtimestamp(
                    path.stat().st_mtime, UTC
                ).isoformat(),
            }
            for path in artifacts
            if path.exists()
        },
    }
    _write_json(paths.root / "results" / "provenance_manifest.json", manifest)
    return manifest


def run_all(
    root: Path,
    *,
    resume: bool = True,
    max_records: int | None = None,
    max_variants: int | None = None,
    bootstrap_replicates: int = 2_000,
) -> dict[str, object]:
    """Run retrieval, scoring, and reporting."""
    paths = ProjectPaths.from_root(root)
    paths.create()
    variants = retrieve_clinvar(paths, resume=resume, max_records=max_records)
    scored = score_callset(paths, variants, resume=resume, max_variants=max_variants)
    return create_report(paths, scored, bootstrap_replicates=bootstrap_replicates)
