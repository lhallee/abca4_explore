"""Clinical discrimination and functional-assay statistics."""

from __future__ import annotations


import numpy as np
import pandas as pd

from scipy import stats
from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score


MIN_CLASS_SIZE = 10


def _finite_binary_data(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    subset = frame.dropna(subset=["binary_label", "avi_phred"])
    labels = subset["binary_label"].astype(int).to_numpy()  # (n,)
    scores = subset["avi_phred"].astype(float).to_numpy()  # (n,)
    return labels, scores


def _threshold_metrics(
    labels: np.ndarray, scores: np.ndarray, threshold: float
) -> dict[str, float | int | None]:
    predictions = (scores >= threshold).astype(int)  # (n,)
    tn, fp, fn, tp = confusion_matrix(labels, predictions, labels=[0, 1]).ravel()

    def divide(numerator: float, denominator: float) -> float | None:
        return float(numerator / denominator) if denominator else None

    sensitivity = divide(tp, tp + fn)
    specificity = divide(tn, tn + fp)
    precision = divide(tp, tp + fp)
    negative_predictive_value = divide(tn, tn + fn)
    f1 = divide(2 * tp, 2 * tp + fp + fn)
    balanced_accuracy = (
        (sensitivity + specificity) / 2
        if sensitivity is not None and specificity is not None
        else None
    )
    return {
        "threshold": threshold,
        "true_positive": int(tp),
        "false_positive": int(fp),
        "true_negative": int(tn),
        "false_negative": int(fn),
        "sensitivity": sensitivity,
        "specificity": specificity,
        "positive_predictive_value": precision,
        "negative_predictive_value": negative_predictive_value,
        "f1": f1,
        "balanced_accuracy": balanced_accuracy,
    }


def _point_metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    precision, recall, _ = precision_recall_curve(labels, scores)
    positives = scores[labels == 1]
    negatives = scores[labels == 0]
    mann_whitney = stats.mannwhitneyu(positives, negatives, alternative="two-sided")
    rank_biserial = (
        2.0 * float(mann_whitney.statistic) / (len(positives) * len(negatives)) - 1.0
    )
    point_biserial = stats.pointbiserialr(labels, scores)
    return {
        "roc_auc": float(roc_auc_score(labels, scores)),
        "average_precision": float(average_precision_score(labels, scores)),
        "pr_auc_trapezoidal": float(auc(recall, precision)),
        "prevalence": float(labels.mean()),
        "point_biserial_r": float(point_biserial.statistic),
        "point_biserial_p": float(point_biserial.pvalue),
        "mann_whitney_u": float(mann_whitney.statistic),
        "mann_whitney_p": float(mann_whitney.pvalue),
        "rank_biserial": rank_biserial,
    }


def _bootstrap_intervals(
    labels: np.ndarray,
    scores: np.ndarray,
    *,
    replicates: int,
    seed: int,
) -> dict[str, list[float]]:
    rng = np.random.default_rng(seed)
    positive_scores = scores[labels == 1]
    negative_scores = scores[labels == 0]
    sampled: dict[str, list[float]] = {
        "roc_auc": [],
        "average_precision": [],
        "pr_auc_trapezoidal": [],
        "point_biserial_r": [],
        "rank_biserial": [],
    }
    for _ in range(replicates):
        bootstrap_positive = rng.choice(
            positive_scores, size=len(positive_scores), replace=True
        )  # (n_positive,)
        bootstrap_negative = rng.choice(
            negative_scores, size=len(negative_scores), replace=True
        )  # (n_negative,)
        bootstrap_scores = np.concatenate(
            (bootstrap_negative, bootstrap_positive)
        )  # (n,)
        bootstrap_labels = np.concatenate(
            (
                np.zeros(len(bootstrap_negative), dtype=int),
                np.ones(len(bootstrap_positive), dtype=int),
            )
        )  # (n,)
        metrics = _point_metrics(bootstrap_labels, bootstrap_scores)
        for name in sampled:
            sampled[name].append(metrics[name])
    return {
        name: [
            float(np.percentile(values, 2.5)),
            float(np.percentile(values, 97.5)),
        ]
        for name, values in sampled.items()
    }


def evaluate_cohort(
    frame: pd.DataFrame,
    *,
    name: str,
    bootstrap_replicates: int = 2_000,
    seed: int = 20_260_911,
) -> dict[str, object]:
    """Evaluate one pathogenic-versus-benign cohort."""
    labels, scores = _finite_binary_data(frame)
    n_positive = int((labels == 1).sum())
    n_negative = int((labels == 0).sum())
    output: dict[str, object] = {
        "name": name,
        "n": int(len(labels)),
        "n_positive": n_positive,
        "n_negative": n_negative,
    }
    if n_positive < MIN_CLASS_SIZE or n_negative < MIN_CLASS_SIZE:
        output["status"] = "insufficient_class_size"
        return output

    output["status"] = "ok"
    output.update(_point_metrics(labels, scores))
    output["confidence_intervals_95"] = _bootstrap_intervals(
        labels,
        scores,
        replicates=bootstrap_replicates,
        seed=seed,
    )
    output["thresholds"] = [
        _threshold_metrics(labels, scores, threshold) for threshold in (20.0, 30.0)
    ]
    return output


def build_cohorts(frame: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Build primary, sensitivity, and region-specific binary cohorts."""
    binary = frame[frame["clinical_class"].isin(["P/LP", "B/LB"])].copy()
    cohorts = {
        "primary_at_least_one_star": binary[binary["review_stars"] >= 1],
        "all_review_statuses": binary,
        "stgd1_at_least_one_star": binary[
            (binary["review_stars"] >= 1) & binary["is_stgd1"]
        ],
    }
    for region_class in (
        "coding_missense",
        "canonical_splice",
        "splice_region",
        "deep_intronic",
    ):
        cohorts[f"primary_{region_class}"] = cohorts["primary_at_least_one_star"][
            cohorts["primary_at_least_one_star"]["region_class"] == region_class
        ]
    return cohorts


def evaluate_all_cohorts(
    frame: pd.DataFrame,
    *,
    bootstrap_replicates: int = 2_000,
) -> dict[str, object]:
    """Evaluate every prespecified analysis cohort."""
    return {
        name: evaluate_cohort(
            cohort,
            name=name,
            bootstrap_replicates=bootstrap_replicates,
        )
        for name, cohort in build_cohorts(frame).items()
    }


def evaluate_functional_assays(frame: pd.DataFrame) -> list[dict[str, object]]:
    """Correlate comparable published assay measures with AVI fields."""
    outputs: list[dict[str, object]] = []
    required = {"assay_group", "effect_value", "avi_phred"}
    if frame.empty or not required.issubset(frame.columns):
        return outputs
    group_columns = ["study_key", "assay_group", "effect_name", "effect_unit"]
    for keys, group in frame.dropna(subset=list(required)).groupby(group_columns):
        study_key, assay_group, effect_name, effect_unit = keys
        target = (
            "fi_ALPHAMISSENSE"
            if str(assay_group).startswith("protein")
            else "fi_MERGED_SPLICING"
        )
        if len(group) < 3 or target not in group:
            continue
        for predictor in ("avi_phred", target):
            subset = group.dropna(subset=[predictor, "effect_value"])
            if len(subset) < 3:
                continue
            association = stats.spearmanr(
                subset[predictor].astype(float),
                subset["effect_value"].astype(float),
            )
            outputs.append(
                {
                    "study_key": study_key,
                    "assay_group": assay_group,
                    "effect_name": effect_name,
                    "effect_unit": effect_unit,
                    "predictor": predictor,
                    "n": int(len(subset)),
                    "spearman_rho": float(association.statistic),
                    "spearman_p": float(association.pvalue),
                }
            )
    return outputs
