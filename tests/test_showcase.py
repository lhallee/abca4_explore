from __future__ import annotations

import pandas as pd

from abca4_avi.showcase import compute_signal_diagnostics


def test_signal_diagnostics_preserve_cohort_counts() -> None:
    rows = []
    for index, clinical_class in enumerate(("B/LB", "B/LB", "P/LP", "P/LP")):
        row = {
            "clinical_class": clinical_class,
            "review_stars": 1,
            "region_class": "deep_intronic",
            "avi_phred": float(index),
        }
        for feature in (
            "fi_MAX_ABS_ATAC",
            "fi_MAX_ABS_CAGE",
            "fi_MAX_ABS_CHIP_HISTONE",
            "fi_MAX_ABS_CHIP_TF",
            "fi_MAX_ABS_CONTACT_MAPS",
            "fi_MAX_ABS_DNASE",
            "fi_MAX_ABS_POLYADENYLATION",
            "fi_MAX_ABS_PROCAP",
            "fi_MAX_ABS_RNA_SEQ",
            "fi_MERGED_SPLICING",
            "fi_ALPHAMISSENSE",
            "fi_CACTUS_241_WAY",
            "fi_PHASTCONS_470_WAY",
        ):
            row[feature] = float(index)
        rows.append(row)

    diagnostics = compute_signal_diagnostics(pd.DataFrame(rows))
    primary = diagnostics[diagnostics["cohort"] == "all"]

    assert set(primary["n"]) == {4}
    assert set(primary["n_positive"]) == {2}
    assert set(primary["n_negative"]) == {2}
    assert set(primary["roc_auc"]) == {1.0}
