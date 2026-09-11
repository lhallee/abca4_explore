"""Static 300 dpi figures for the ABCA4 AVI analysis."""

from __future__ import annotations

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path

from sklearn.metrics import precision_recall_curve, roc_curve

from .atlas import AVI_FEATURES
from .statistics import build_cohorts


DPI = 300
CLASS_ORDER = ["B/LB", "VUS", "Conflicting", "P/LP", "Other"]
REGION_ORDER = [
    "coding_missense",
    "canonical_splice",
    "splice_region",
    "deep_intronic",
    "other_exonic",
    "other",
]


def _save(figure: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def configure_style() -> None:
    """Set a compact, colorblind-safe plotting style."""
    sns.set_theme(style="whitegrid", context="notebook")
    matplotlib.rcParams.update({"figure.dpi": 120, "savefig.dpi": DPI})


def plot_class_counts(frame: pd.DataFrame, path: Path) -> None:
    counts = frame["clinical_class"].value_counts().reindex(CLASS_ORDER, fill_value=0)
    figure, axis = plt.subplots(figsize=(7.2, 4.5))
    sns.barplot(x=counts.index, y=counts.values, color="#3B82F6", ax=axis)
    axis.set(xlabel="ClinVar aggregate class", ylabel="Records")
    axis.set_title("ABCA4 ClinVar records by aggregate classification")
    for index, value in enumerate(counts.values):
        axis.text(index, value, f"{value:,}", ha="center", va="bottom", fontsize=8)
    _save(figure, path)


def plot_score_distributions(frame: pd.DataFrame, path: Path) -> None:
    subset = frame.dropna(subset=["clinical_class", "avi_phred"])
    order = [name for name in CLASS_ORDER if name in set(subset["clinical_class"])]
    figure, axis = plt.subplots(figsize=(8.2, 5.0))
    sns.violinplot(
        data=subset,
        x="clinical_class",
        y="avi_phred",
        order=order,
        inner=None,
        cut=0,
        color="#93C5FD",
        ax=axis,
    )
    sns.boxplot(
        data=subset,
        x="clinical_class",
        y="avi_phred",
        order=order,
        width=0.22,
        showfliers=False,
        color="white",
        ax=axis,
    )
    axis.axhline(20, color="#F59E0B", linestyle="--", linewidth=1, label="Phred 20")
    axis.axhline(30, color="#DC2626", linestyle=":", linewidth=1, label="Phred 30")
    axis.set(xlabel="ClinVar aggregate class", ylabel="AVI Phred")
    axis.set_title("AlphaGenome AVI score distributions")
    axis.legend(frameon=False)
    _save(figure, path)


def plot_roc_pr(frame: pd.DataFrame, path: Path) -> None:
    cohort = build_cohorts(frame)["primary_at_least_one_star"].dropna(
        subset=["binary_label", "avi_phred"]
    )
    labels = cohort["binary_label"].astype(int).to_numpy()  # (n,)
    scores = cohort["avi_phred"].astype(float).to_numpy()  # (n,)
    figure, axes = plt.subplots(1, 2, figsize=(10.0, 4.4))
    if len(np.unique(labels)) == 2:
        false_positive_rate, true_positive_rate, _ = roc_curve(labels, scores)
        precision, recall, _ = precision_recall_curve(labels, scores)
        axes[0].plot(false_positive_rate, true_positive_rate, color="#2563EB")
        axes[0].plot([0, 1], [0, 1], color="#9CA3AF", linestyle="--")
        axes[1].plot(recall, precision, color="#7C3AED")
        axes[1].axhline(labels.mean(), color="#9CA3AF", linestyle="--")
    axes[0].set(
        xlabel="False-positive rate",
        ylabel="True-positive rate",
        title="ROC curve",
        xlim=(0, 1),
        ylim=(0, 1),
    )
    axes[1].set(
        xlabel="Recall",
        ylabel="Precision",
        title="Precision-recall curve",
        xlim=(0, 1),
        ylim=(0, 1),
    )
    figure.suptitle("Primary ClinVar cohort: P/LP versus B/LB")
    _save(figure, path)


def plot_region_scores(frame: pd.DataFrame, path: Path) -> None:
    subset = frame.dropna(subset=["region_class", "avi_phred"])
    order = [name for name in REGION_ORDER if name in set(subset["region_class"])]
    figure, axis = plt.subplots(figsize=(10.0, 5.0))
    sns.boxplot(
        data=subset,
        x="region_class",
        y="avi_phred",
        hue="clinical_class",
        order=order,
        showfliers=False,
        ax=axis,
    )
    axis.set(xlabel="Variant region", ylabel="AVI Phred")
    axis.set_title("AVI scores by molecular region and ClinVar class")
    axis.tick_params(axis="x", rotation=25)
    axis.legend(title="ClinVar", frameon=False, ncol=2)
    _save(figure, path)


def plot_feature_heatmap(frame: pd.DataFrame, path: Path) -> None:
    feature_columns = [f"fi_{name}" for name in AVI_FEATURES if f"fi_{name}" in frame]
    subset = frame.dropna(subset=["region_class"])
    summary = subset.groupby("region_class")[feature_columns].apply(
        lambda values: values.abs().mean()
    )
    summary = summary.reindex([name for name in REGION_ORDER if name in summary.index])
    summary.columns = [column.removeprefix("fi_") for column in summary.columns]
    figure, axis = plt.subplots(figsize=(13.0, 4.8))
    sns.heatmap(summary, cmap="mako", ax=axis)
    axis.set(xlabel="AVI feature", ylabel="Variant region")
    axis.set_title("Mean absolute AVI feature attribution")
    axis.tick_params(axis="x", rotation=55, labelsize=7)
    _save(figure, path)


def plot_sentinel_variants(frame: pd.DataFrame, path: Path) -> None:
    subset = frame.dropna(subset=["sentinel_name", "avi_phred"]).copy()
    figure, axis = plt.subplots(figsize=(8.5, 4.5))
    if not subset.empty:
        subset = subset.sort_values("avi_phred")
        sns.barplot(
            data=subset,
            y="sentinel_name",
            x="avi_phred",
            hue="top_modality",
            dodge=False,
            ax=axis,
        )
        axis.axvline(20, color="#F59E0B", linestyle="--", linewidth=1)
        axis.axvline(30, color="#DC2626", linestyle=":", linewidth=1)
        axis.legend(title="Top modality", frameon=False)
    axis.set(xlabel="AVI Phred", ylabel="Variant")
    axis.set_title("Required ABCA4 sentinel variants")
    _save(figure, path)


def plot_functional_assays(frame: pd.DataFrame, path: Path) -> None:
    subset = frame.dropna(subset=["effect_value", "avi_phred"])
    if subset.empty:
        figure, axis = plt.subplots(figsize=(7.0, 3.5))
        axis.text(0.5, 0.5, "No assay variants matched scored SNVs", ha="center")
        axis.set_axis_off()
        _save(figure, path)
        return
    assay_groups = list(dict.fromkeys(subset["assay_group"].astype(str)))
    figure, axes = plt.subplots(
        len(assay_groups),
        2,
        figsize=(10.5, max(4.5, 3.8 * len(assay_groups))),
        squeeze=False,
        layout="constrained",
    )
    for row_index, assay_group in enumerate(assay_groups):
        group = subset[subset["assay_group"].astype(str) == assay_group]
        feature = (
            "fi_ALPHAMISSENSE"
            if assay_group.startswith("protein")
            else "fi_MERGED_SPLICING"
        )
        for column_index, (axis, predictor) in enumerate(
            zip(axes[row_index], ("avi_phred", feature), strict=True)
        ):
            if predictor in group:
                sns.scatterplot(
                    data=group,
                    x=predictor,
                    y="effect_value",
                    hue="study_key",
                    ax=axis,
                )
                if column_index == 1 and axis.get_legend() is not None:
                    axis.get_legend().remove()
            axis.set(
                xlabel=predictor.removeprefix("fi_"),
                ylabel="Reported effect",
                title=(
                    f"{assay_group.replace('_', ' ')} versus "
                    f"{predictor.removeprefix('fi_')}"
                ),
            )
    _save(figure, path)


def generate_plots(
    all_clinvar: pd.DataFrame,
    scored: pd.DataFrame,
    sentinels: pd.DataFrame,
    functional: pd.DataFrame,
    output_dir: Path,
) -> list[Path]:
    """Generate every prespecified figure and return its path."""
    configure_style()
    plot_paths = [
        output_dir / "clinvar_class_counts.png",
        output_dir / "avi_by_clinvar_class.png",
        output_dir / "roc_pr_curves.png",
        output_dir / "avi_by_variant_region.png",
        output_dir / "feature_attribution_heatmap.png",
        output_dir / "sentinel_variants.png",
    ]
    plot_class_counts(all_clinvar, plot_paths[0])
    plot_score_distributions(scored, plot_paths[1])
    plot_roc_pr(scored, plot_paths[2])
    plot_region_scores(scored, plot_paths[3])
    plot_feature_heatmap(scored, plot_paths[4])
    plot_sentinel_variants(sentinels, plot_paths[5])
    if not functional.empty and "effect_value" in functional:
        functional_path = output_dir / "functional_assay_comparison.png"
        plot_functional_assays(functional, functional_path)
        plot_paths.append(functional_path)
    return plot_paths
