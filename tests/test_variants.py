from __future__ import annotations

import pandas as pd

from abca4_avi.variants import deduplicate_eligible_variants, normalize_clinical_class


def test_clinical_classification_mapping() -> None:
    assert normalize_clinical_class("Pathogenic/Likely pathogenic") == "P/LP"
    assert normalize_clinical_class("Benign/Likely benign") == "B/LB"
    assert normalize_clinical_class("Uncertain significance") == "VUS"
    assert normalize_clinical_class("Conflicting classifications") == "Conflicting"
    assert normalize_clinical_class(pd.NA) == "Other"


def test_deduplication_merges_ids_and_marks_cross_record_conflict() -> None:
    frame = pd.DataFrame(
        [
            {
                "variant": "chr1:100:A>G",
                "eligible_for_avi": True,
                "clinical_class": "P/LP",
                "binary_label": 1,
                "review_stars": 1,
                "clinvar_variation_id": "2",
                "is_stgd1": False,
            },
            {
                "variant": "chr1:100:A>G",
                "eligible_for_avi": True,
                "clinical_class": "B/LB",
                "binary_label": 0,
                "review_stars": 2,
                "clinvar_variation_id": "1",
                "is_stgd1": True,
            },
        ]
    )

    result = deduplicate_eligible_variants(frame)

    assert len(result) == 1
    assert result.loc[0, "clinical_class"] == "Conflicting"
    assert pd.isna(result.loc[0, "binary_label"])
    assert result.loc[0, "clinvar_variation_ids"] == "1|2"
    assert result.loc[0, "review_stars"] == 2
    assert bool(result.loc[0, "is_stgd1"])
