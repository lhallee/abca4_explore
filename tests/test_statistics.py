from __future__ import annotations

import pandas as pd
import pytest

from abca4_avi.statistics import evaluate_cohort


def test_metrics_separate_perfectly_ranked_classes() -> None:
    frame = pd.DataFrame(
        {
            "binary_label": [0] * 10 + [1] * 10,
            "avi_phred": list(range(10)) + list(range(20, 30)),
        }
    )

    metrics = evaluate_cohort(frame, name="test", bootstrap_replicates=20, seed=7)

    assert metrics["status"] == "ok"
    assert metrics["roc_auc"] == 1.0
    assert metrics["average_precision"] == pytest.approx(1.0)
    assert metrics["rank_biserial"] == 1.0
    assert len(metrics["confidence_intervals_95"]["roc_auc"]) == 2
