"""Find population-compatible ABCA4 missense controls in gnomAD v4."""

from __future__ import annotations

import argparse
import gzip
import json

from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from abca4_avi.population_controls import (
    ABCA4_BA1_AF,
    ABCA4_BS1_STRONG_AF,
    ABCA4_BS1_SUPPORTING_AF,
    GNOMAD_DATASET,
    fetch_abca4_gnomad,
    join_existing_screen,
    normalize_gnomad_missense,
    save_raw_response,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="ABCA4 analysis root",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Replace the cached raw gnomAD response",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    raw_path = root / "data" / "raw" / "gnomad_abca4_v4.json.gz"
    table_root = root / "results" / "tables"
    mapping_path = root / "results" / "atlas_ppi" / "variant_protein_mapping.parquet"

    if raw_path.exists() and not args.refresh:
        with gzip.open(raw_path, "rt", encoding="utf-8") as handle:
            response = json.load(handle)
    else:
        response = fetch_abca4_gnomad()
        save_raw_response(response, raw_path)
    missense = normalize_gnomad_missense(response)
    mapping = pd.read_parquet(mapping_path)
    candidates = join_existing_screen(missense, mapping)

    table_root.mkdir(parents=True, exist_ok=True)
    output_stem = table_root / "gnomad_population_control_candidates"
    candidates.to_csv(output_stem.with_suffix(".tsv"), sep="\t", index=False)
    candidates.to_parquet(output_stem.with_suffix(".parquet"), index=False)

    selected = candidates.population_evidence_tier.str.endswith("candidate")
    summary = {
        "candidate_label_scope": "population-compatible; not clinical benign",
        "dataset": GNOMAD_DATASET,
        "retrieved_at": datetime.now(UTC).isoformat(),
        "thresholds": {
            "ba1": ABCA4_BA1_AF,
            "bs1_strong": ABCA4_BS1_STRONG_AF,
            "bs1_supporting": ABCA4_BS1_SUPPORTING_AF,
        },
        "grpmax_method": (
            "maximum joint AF after excluding gnomAD v4 bottlenecked and "
            "remaining ancestry groups"
        ),
        "counts": {
            "gnomad_missense_snvs": int(len(candidates)),
            "frequency_candidates": int(selected.sum()),
            "vcep_excluded": int(candidates.vcep_frequency_exclusion.sum()),
            "already_in_clinvar_snapshot": int(
                candidates.loc[selected, "already_in_clinvar_snapshot"].sum()
            ),
            "already_screened_by_atlas_ppi": int(
                candidates.loc[selected, "already_screened_by_atlas_ppi"].sum()
            ),
            "new_sequences_requiring_embedding": int(
                (~candidates.loc[selected, "already_screened_by_atlas_ppi"]).sum()
            ),
        },
        "tiers": {
            str(key): int(value)
            for key, value in candidates.population_evidence_tier.value_counts().items()
        },
    }
    summary_path = output_stem.with_name(f"{output_stem.name}_summary.json")
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
