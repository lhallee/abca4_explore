"""Leakage-aware analysis and figures for the ABCA4 AVI results."""

from __future__ import annotations

import json

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from pathlib import Path

from scipy.stats import spearmanr
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)


DPI = 300
BLUE = "#1F4E79"
ORANGE = "#D97732"
GOLD = "#D5A021"
PINK = "#C65A7A"
GRAY = "#667085"
LIGHT_GRAY = "#D0D5DD"
CLASS_COLORS = {
    "B/LB": "#4E79A7",
    "VUS": "#D5A021",
    "Conflicting": "#B07AA1",
    "P/LP": "#D97732",
}
CLASS_STYLES = {
    "B/LB": "-",
    "VUS": "--",
    "Conflicting": ":",
    "P/LP": "-.",
}
ALPHAGENOME_FEATURES = (
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
)
SIGNAL_DEFINITIONS = {
    "Full AVI": ("avi_phred",),
    "All AlphaGenome attributions": ALPHAGENOME_FEATURES,
    "AlphaMissense attribution": ("fi_ALPHAMISSENSE",),
    "Conservation attributions": (
        "fi_CACTUS_241_WAY",
        "fi_PHASTCONS_470_WAY",
    ),
    "Splicing attribution": ("fi_MERGED_SPLICING",),
}
COHORT_LABELS = {
    "all": "All primary",
    "coding_missense": "Coding missense",
    "splice_region": "Splice region",
    "deep_intronic": "Deep intronic",
}


def _configure_style() -> None:
    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.edgecolor": "#344054",
            "axes.labelcolor": "#344054",
            "axes.titlecolor": "#101828",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "text.color": "#101828",
        }
    )


def _save(figure: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _primary_cohort(scored: pd.DataFrame) -> pd.DataFrame:
    return scored[
        (scored["review_stars"] >= 1) & scored["clinical_class"].isin(["P/LP", "B/LB"])
    ].copy()


def compute_signal_diagnostics(scored: pd.DataFrame) -> pd.DataFrame:
    """Measure ranking from full AVI and groups of signed AVI attributions."""
    primary = _primary_cohort(scored)
    primary["binary_label"] = (primary["clinical_class"] == "P/LP").astype(int)
    cohorts = {"all": primary}
    cohorts.update(
        {
            name: primary[primary["region_class"] == name]
            for name in COHORT_LABELS
            if name != "all"
        }
    )

    records: list[dict[str, object]] = []
    for cohort_name, cohort in cohorts.items():
        labels = cohort["binary_label"].to_numpy(dtype=int)  # (n,)
        if len(np.unique(labels)) != 2:
            continue
        for signal_name, columns in SIGNAL_DEFINITIONS.items():
            scores = cohort[list(columns)].sum(axis=1).to_numpy(dtype=float)  # (n,)
            records.append(
                {
                    "cohort": cohort_name,
                    "cohort_label": COHORT_LABELS[cohort_name],
                    "signal": signal_name,
                    "n": len(cohort),
                    "n_positive": int(labels.sum()),
                    "n_negative": int((1 - labels).sum()),
                    "roc_auc": roc_auc_score(labels, scores),
                    "average_precision": average_precision_score(labels, scores),
                }
            )
    return pd.DataFrame.from_records(records)


def compute_late_evaluated_metrics(
    scored: pd.DataFrame,
    *,
    cutoff: str = "2025-06-15",
) -> dict[str, float | int | str]:
    """Evaluate records last assessed after the Atlas paper's ClinVar snapshot."""
    primary = _primary_cohort(scored)
    evaluated = pd.to_datetime(primary["last_evaluated"], errors="coerce")
    cohort = primary[evaluated > pd.Timestamp(cutoff)].copy()
    labels = (cohort["clinical_class"] == "P/LP").to_numpy(dtype=int)  # (n,)
    scores = cohort["avi_phred"].to_numpy(dtype=float)  # (n,)
    return {
        "cutoff": cutoff,
        "n": len(cohort),
        "n_positive": int(labels.sum()),
        "n_negative": int((1 - labels).sum()),
        "roc_auc": roc_auc_score(labels, scores),
        "average_precision": average_precision_score(labels, scores),
    }


def _bootstrap_spearman(
    x: np.ndarray,
    y: np.ndarray,
    *,
    replicates: int,
    random_seed: int,
) -> tuple[float, float]:
    # x: (n,); y: (n,)
    rng = np.random.default_rng(random_seed)
    estimates: list[float] = []
    n = len(x)
    for _ in range(replicates):
        indices = rng.integers(0, n, size=n)  # (n,)
        estimate = spearmanr(x[indices], y[indices]).statistic
        if np.isfinite(estimate):
            estimates.append(float(estimate))
    lower, upper = np.quantile(np.asarray(estimates), [0.025, 0.975])
    return float(lower), float(upper)


def compute_functional_correlations(
    functional: pd.DataFrame,
    *,
    bootstrap_replicates: int = 2_000,
    random_seed: int = 20260911,
) -> pd.DataFrame:
    """Calculate study-specific assay correlations with bootstrap intervals."""
    predictor_by_assay = {
        "protein_basal_atpase": ("avi_phred", "fi_ALPHAMISSENSE"),
        "protein_expression": ("avi_phred", "fi_ALPHAMISSENSE"),
        "protein_f_index": ("avi_phred", "fi_ALPHAMISSENSE"),
        "protein_ret_pe_induced_atpase": (
            "avi_phred",
            "fi_ALPHAMISSENSE",
        ),
        "splicing_correct_mrna": ("avi_phred", "fi_MERGED_SPLICING"),
    }
    assay_labels = {
        "protein_basal_atpase": "Basal ATPase",
        "protein_expression": "Protein expression",
        "protein_f_index": "F-index",
        "protein_ret_pe_induced_atpase": "Retinoid-stimulated ATPase",
        "splicing_correct_mrna": "Correctly spliced mRNA",
    }
    study_labels = {
        "aslaksen_2024": "Aslaksen 2024",
        "garces_2021": "Garces 2021",
        "sangermano_2018": "Sangermano 2018",
    }
    predictor_labels = {
        "avi_phred": "Full AVI",
        "fi_ALPHAMISSENSE": "AlphaMissense attribution",
        "fi_MERGED_SPLICING": "Splicing attribution",
    }

    records: list[dict[str, object]] = []
    grouped = functional.dropna(subset=["effect_value", "avi_phred"]).groupby(
        ["study_key", "assay_group"], sort=True
    )
    for group_index, ((study_key, assay_group), group) in enumerate(grouped):
        for predictor_index, predictor in enumerate(predictor_by_assay[assay_group]):
            valid = group.dropna(subset=[predictor, "effect_value"])
            x = valid[predictor].to_numpy(dtype=float)  # (n,)
            y = valid["effect_value"].to_numpy(dtype=float)  # (n,)
            estimate = float(spearmanr(x, y).statistic)
            lower, upper = _bootstrap_spearman(
                x,
                y,
                replicates=bootstrap_replicates,
                random_seed=random_seed + group_index * 10 + predictor_index,
            )
            records.append(
                {
                    "study_key": study_key,
                    "study": study_labels[study_key],
                    "assay_group": assay_group,
                    "assay": assay_labels[assay_group],
                    "predictor": predictor,
                    "predictor_label": predictor_labels[predictor],
                    "n": len(valid),
                    "spearman_rho": estimate,
                    "ci_lower": lower,
                    "ci_upper": upper,
                }
            )
    return pd.DataFrame.from_records(records)


def plot_overview(
    scored: pd.DataFrame,
    clinical_metrics: dict[str, dict[str, object]],
    path: Path,
) -> None:
    """Plot clinical separation and subgroup performance with uncertainty."""
    _configure_style()
    primary = _primary_cohort(scored)
    labels = (primary["clinical_class"] == "P/LP").to_numpy(dtype=int)  # (n,)
    scores = primary["avi_phred"].to_numpy(dtype=float)  # (n,)
    false_positive_rate, true_positive_rate, _ = roc_curve(labels, scores)
    precision, recall, _ = precision_recall_curve(labels, scores)
    metric = clinical_metrics["primary_at_least_one_star"]
    intervals = metric["confidence_intervals_95"]

    figure, axes = plt.subplots(2, 2, figsize=(14.0, 10.5), layout="constrained")
    distribution_axis, roc_axis, pr_axis, subgroup_axis = axes.ravel()

    display = scored[
        scored["clinical_class"].isin(CLASS_COLORS) & scored["avi_phred"].notna()
    ]
    for clinical_class in CLASS_COLORS:
        values = display.loc[
            display["clinical_class"] == clinical_class, "avi_phred"
        ].to_numpy(dtype=float)  # (n_class,)
        x = np.sort(values)  # (n_class,)
        y = np.arange(1, len(x) + 1) / len(x)  # (n_class,)
        distribution_axis.plot(
            x,
            y,
            color=CLASS_COLORS[clinical_class],
            linestyle=CLASS_STYLES[clinical_class],
            linewidth=2.2,
            label=f"{clinical_class} (n={len(x):,})",
        )
    for threshold, line_style in ((20, "--"), (30, ":")):
        distribution_axis.axvline(
            threshold,
            color=GRAY,
            linestyle=line_style,
            linewidth=1.4,
        )
    distribution_axis.text(20.4, 0.04, "Top 1%", color=GRAY, fontsize=9)
    distribution_axis.text(30.4, 0.04, "Top 0.1%", color=GRAY, fontsize=9)
    distribution_axis.set(
        title="A. AVI score distributions by ClinVar class",
        xlabel="AVI Phred score",
        ylabel="Cumulative fraction",
        xlim=(0, 52),
        ylim=(0, 1),
    )
    distribution_axis.legend(frameon=False, fontsize=10, loc="lower right")

    roc_axis.plot(false_positive_rate, true_positive_rate, color=BLUE, linewidth=2.6)
    roc_axis.plot([0, 1], [0, 1], color=LIGHT_GRAY, linestyle="--")
    roc_ci = intervals["roc_auc"]
    roc_axis.set(
        title="B. ROC curve",
        xlabel="False-positive rate",
        ylabel="True-positive rate",
        xlim=(0, 1),
        ylim=(0, 1),
    )
    roc_axis.legend(
        [f"ROC-AUC {metric['roc_auc']:.3f} (95% CI {roc_ci[0]:.3f}-{roc_ci[1]:.3f})"],
        frameon=False,
        loc="lower right",
        fontsize=10,
    )

    pr_axis.plot(recall, precision, color=ORANGE, linewidth=2.6)
    pr_axis.axhline(labels.mean(), color=LIGHT_GRAY, linestyle="--")
    ap_ci = intervals["average_precision"]
    pr_axis.set(
        title="C. Precision-recall curve",
        xlabel="Recall",
        ylabel="Precision",
        xlim=(0, 1),
        ylim=(0, 1),
    )
    pr_axis.legend(
        [
            f"Average precision {metric['average_precision']:.3f} "
            f"(95% CI {ap_ci[0]:.3f}-{ap_ci[1]:.3f})"
        ],
        frameon=False,
        loc="lower left",
        fontsize=10,
    )

    cohort_order = (
        "primary_at_least_one_star",
        "stgd1_at_least_one_star",
        "primary_coding_missense",
        "primary_splice_region",
        "primary_deep_intronic",
    )
    labels_by_cohort = {
        "primary_at_least_one_star": "All primary",
        "stgd1_at_least_one_star": "STGD1-specific",
        "primary_coding_missense": "Coding missense",
        "primary_splice_region": "Splice region",
        "primary_deep_intronic": "Deep intronic",
    }
    positions = np.arange(len(cohort_order))  # (n_cohorts,)
    for offset, metric_name, color, marker, label in (
        (-0.10, "roc_auc", BLUE, "o", "ROC-AUC"),
        (0.10, "average_precision", ORANGE, "s", "Average precision"),
    ):
        estimates = np.asarray(
            [clinical_metrics[name][metric_name] for name in cohort_order]
        )  # (n_cohorts,)
        bounds = np.asarray(
            [
                clinical_metrics[name]["confidence_intervals_95"][metric_name]
                for name in cohort_order
            ]
        )  # (n_cohorts, 2)
        errors = np.vstack(
            (estimates - bounds[:, 0], bounds[:, 1] - estimates)
        )  # (2, n_cohorts)
        subgroup_axis.errorbar(
            estimates,
            positions + offset,
            xerr=errors,
            fmt=marker,
            color=color,
            capsize=3,
            linewidth=1.5,
            markersize=7,
            label=label,
        )
    subgroup_axis.axvline(0.5, color=LIGHT_GRAY, linestyle="--")
    subgroup_axis.set(
        title="D. Performance by ABCA4 cohort",
        xlabel="Metric value with 95% bootstrap CI",
        xlim=(0.45, 1.01),
        yticks=positions,
        yticklabels=[labels_by_cohort[name] for name in cohort_order],
    )
    subgroup_axis.invert_yaxis()
    subgroup_axis.legend(
        frameon=False,
        loc="upper left",
        fontsize=10,
    )

    figure.suptitle(
        "ABCA4 AlphaGenome AVI results",
        fontsize=22,
        fontweight="bold",
    )
    _save(figure, path)


def plot_signal_sources(diagnostics: pd.DataFrame, path: Path) -> None:
    """Plot which AVI attribution groups carry rank information by region."""
    _configure_style()
    signal_order = list(SIGNAL_DEFINITIONS)
    cohort_order = list(COHORT_LABELS)
    matrix = (
        diagnostics.pivot(index="signal", columns="cohort", values="roc_auc")
        .reindex(index=signal_order, columns=cohort_order)
        .rename(columns=COHORT_LABELS)
    )
    matrix.columns = [label.replace(" ", "\n") for label in matrix.columns]
    figure, axis = plt.subplots(figsize=(12.0, 6.2), layout="constrained")
    sns.heatmap(
        matrix,
        vmin=0.5,
        vmax=1.0,
        cmap=sns.light_palette(BLUE, as_cmap=True),
        annot=True,
        fmt=".3f",
        linewidths=1.0,
        linecolor="white",
        cbar_kws={"label": "ROC-AUC"},
        ax=axis,
    )
    axis.set(
        title="ROC-AUC by cohort and AVI attribution source",
        xlabel="ABCA4 primary cohort subset",
        ylabel="",
    )
    axis.tick_params(axis="x", rotation=0)
    axis.tick_params(axis="y", rotation=0)
    _save(figure, path)


def plot_functional_validation(correlations: pd.DataFrame, path: Path) -> None:
    """Plot assay correlations and bootstrap intervals by prediction source."""
    _configure_style()
    correlations = correlations.copy()
    correlations["label"] = correlations.apply(
        lambda row: f"{row['study']}: {row['assay']} (n={row['n']})", axis=1
    )
    labels = list(dict.fromkeys(correlations["label"]))
    positions = {label: index for index, label in enumerate(labels)}
    offsets = {
        "Full AVI": -0.10,
        "AlphaMissense attribution": 0.10,
        "Splicing attribution": 0.10,
    }
    colors = {
        "Full AVI": BLUE,
        "AlphaMissense attribution": ORANGE,
        "Splicing attribution": GOLD,
    }
    markers = {
        "Full AVI": "o",
        "AlphaMissense attribution": "s",
        "Splicing attribution": "D",
    }

    figure, axis = plt.subplots(figsize=(11.5, 7.5), layout="constrained")
    for predictor, group in correlations.groupby("predictor_label", sort=False):
        estimates = group["spearman_rho"].to_numpy(dtype=float)  # (n_assays,)
        lower = group["ci_lower"].to_numpy(dtype=float)  # (n_assays,)
        upper = group["ci_upper"].to_numpy(dtype=float)  # (n_assays,)
        y = np.asarray(
            [positions[label] + offsets[predictor] for label in group["label"]]
        )  # (n_assays,)
        errors = np.vstack((estimates - lower, upper - estimates))  # (2, n_assays)
        axis.errorbar(
            estimates,
            y,
            xerr=errors,
            fmt=markers[predictor],
            color=colors[predictor],
            capsize=3,
            linewidth=1.5,
            markersize=7,
            label=predictor,
        )
    axis.axvline(0, color=GRAY, linewidth=1.2)
    axis.set(
        title="ABCA4 functional assay correlations",
        xlabel=(
            "Spearman correlation with retained molecular function (95% bootstrap CI)"
        ),
        xlim=(-1.02, 1.02),
        yticks=np.arange(len(labels)),
        yticklabels=labels,
    )
    axis.invert_yaxis()
    axis.legend(frameon=False, loc="lower right", fontsize=10)
    _save(figure, path)


def _write_frame(frame: pd.DataFrame, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(stem.with_suffix(".tsv"), sep="\t", index=False)
    frame.to_parquet(stem.with_suffix(".parquet"), index=False)


def build_showcase(root: Path) -> dict[str, object]:
    """Create the leakage-aware report, source tables, and polished figures."""
    root = root.resolve()
    scored = pd.read_parquet(root / "results" / "tables" / "variants_scored.parquet")
    functional = pd.read_parquet(
        root / "results" / "tables" / "functional_assays_scored.parquet"
    )
    clinical_metrics = json.loads(
        (root / "results" / "tables" / "metrics.json").read_text(encoding="utf-8")
    )["clinical"]

    signal_diagnostics = compute_signal_diagnostics(scored)
    correlations = compute_functional_correlations(functional)
    late_metrics = compute_late_evaluated_metrics(scored)
    _write_frame(
        signal_diagnostics,
        root / "results" / "tables" / "showcase_signal_diagnostics",
    )
    _write_frame(
        correlations,
        root / "results" / "tables" / "showcase_functional_correlations",
    )

    output_dir = root / "plots" / "showcase"
    plot_overview(scored, clinical_metrics, output_dir / "abca4_avi_overview.png")
    plot_signal_sources(
        signal_diagnostics,
        output_dir / "abca4_avi_signal_sources.png",
    )
    plot_functional_validation(
        correlations,
        output_dir / "abca4_avi_functional_validation.png",
    )

    primary = clinical_metrics["primary_at_least_one_star"]
    signal_lookup = signal_diagnostics.set_index(["cohort", "signal"])
    splicing_deep = signal_lookup.loc[("deep_intronic", "Splicing attribution")]
    splicing_region = signal_lookup.loc[("splice_region", "Splicing attribution")]
    sentinels = scored[scored["sentinel_name"].notna()].sort_values("sentinel_name")
    atlas_paper = (
        "https://storage.googleapis.com/deepmind-media/DeepMind.com/Blog/"
        "alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-"
        "in-the-human-genome/alphagenome-atlas.pdf"
    )
    report_lines = [
        "# How impressive is the ABCA4 AVI result?",
        "",
        (
            "The result is impressive as a retrospective ABCA4 prioritization "
            "result, especially for splice-region and deep-intronic variants. "
            "The headline 0.990 ROC-AUC is not yet an independent estimate of "
            "prospective performance."
        ),
        "",
        "## What AlphaGenome and AVI trained on",
        "",
        (
            "The base AlphaGenome model learned to predict experimental genome "
            "tracks from human and mouse DNA sequence. It was not trained on "
            "ABCA4 ClinVar classifications."
        ),
        "",
        (
            f"AVI is a separate 18-feature ensemble. According to the "
            f"[Atlas methods paper]({atlas_paper}), it was trained on 37 million "
            "gnomAD v4.1 variants. Variants with filtering allele frequency below "
            "0.001 were proxy-impactful; more common variants were proxy-benign. "
            "Chromosome 1 was one of the AVI training chromosomes."
        ),
        "",
        (
            "Therefore AVI did not simply memorize ClinVar labels, but some ABCA4 "
            "positions may have appeared in its gnomAD-derived training data. The "
            "exact overlap cannot be measured from the public score API because "
            "the sampled 37-million-variant coordinate list is not returned."
        ),
        "",
        (
            "The official Atlas ClinVar benchmark removed every evaluation position "
            "overlapping AVI training or validation. This ABCA4 analysis cannot "
            "apply that filter, so it should be called retrospective separation, "
            "not held-out generalization. Rarity is also partly circular because "
            "population frequency contributes to clinical classification."
        ),
        "",
        "## What remains genuinely persuasive",
        "",
        (
            f"- The primary cohort contains {int(primary['n']):,} variants and "
            f"achieves ROC-AUC {primary['roc_auc']:.3f} and average precision "
            f"{primary['average_precision']:.3f}. The bootstrap intervals are narrow."
        ),
        (
            f"- Among {late_metrics['n']:,} records last evaluated after the "
            f"paper's June 15, 2025 ClinVar snapshot, ROC-AUC remains "
            f"{late_metrics['roc_auc']:.3f} and average precision remains "
            f"{late_metrics['average_precision']:.3f}. Last-evaluated date is not "
            "a first-submission date, so this is a sensitivity analysis, not a "
            "strict temporal holdout."
        ),
        (
            f"- The splicing attribution alone ranks deep-intronic P/LP versus "
            f"B/LB variants at ROC-AUC {splicing_deep['roc_auc']:.3f} and average "
            f"precision {splicing_deep['average_precision']:.3f}. In the splice-region "
            f"subset it reaches ROC-AUC {splicing_region['roc_auc']:.3f}. This is "
            "the most biologically distinctive result because AlphaMissense and "
            "conservation carry little deep-intronic signal."
        ),
        (
            "- Published ABCA4 experiments provide a smaller, more independent "
            "check. In 32 Sangermano variants, the splicing attribution correlates "
            "with correctly spliced mRNA at Spearman rho -0.429. In 35 Garces "
            "missense variants, full AVI correlations with retained protein "
            "expression and ATPase measures range from -0.377 to -0.470."
        ),
        "",
        "## What would make the claim publication-grade",
        "",
        (
            "Obtain the AVI train/validation coordinate list and exclude every "
            "overlapping ABCA4 position, then repeat the analysis in rare-only "
            "benign and pathogenic strata. A truly prospective ClinVar release or "
            "new blinded functional assays would test generalization without the "
            "rarity proxy."
        ),
        "",
        "## Sentinel variants",
        "",
        "| Variant | AVI Phred | Leading attribution | Atlas |",
        "|---|---:|---|---|",
    ]
    for row in sentinels.itertuples(index=False):
        report_lines.append(
            f"| {row.sentinel_name} | {row.avi_phred:.3f} | "
            f"{row.top_modality} | [open]({row.atlas_url}) |"
        )
    report_lines.extend(
        [
            "",
            "AlphaGenome and AVI outputs are research-only molecular predictions, "
            "not clinical diagnoses.",
        ]
    )
    report_path = root / "results" / "impressiveness_assessment.md"
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    return {
        "report": str(report_path),
        "plots": [str(path) for path in sorted(output_dir.glob("*.png"))],
        "late_evaluated": late_metrics,
    }
