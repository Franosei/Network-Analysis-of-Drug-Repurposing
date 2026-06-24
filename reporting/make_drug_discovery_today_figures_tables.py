"""
Generate the focused manuscript figure and table package for the
evidence-readiness benchmark.

The output contract follows the focused prompt:
  - 4 main PNG figures
  - 1 separate graphical abstract
  - 3 main CSV tables only
  - figure captions, table captions, and a 3-paragraph manuscript summary

Usage:
    py reporting/make_drug_discovery_today_figures_tables.py \
        --benchmark-dir outputs/evidence_readiness_benchmark_20260618_204625
"""

from __future__ import annotations

import argparse
import math
import textwrap
from datetime import datetime
from pathlib import Path
from typing import Iterable

import matplotlib
import pandas as pd

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle


REQUIRED_COLUMNS = [
    "drug_disease_pair",
    "source_panel",
    "external_evidence_maturity",
    "expected_pipeline_flag",
    "coverage_tier",
    "literature_count",
    "therapeutic_rate",
    "adverse_rate",
    "irrelevant_rate",
    "gamma_safety_overlap",
    "safety_quality_flag",
    "posterior_mean",
    "credible_interval_width",
    "evidence_readiness_score",
    "translation_evidence_level",
    "observed_quality_flag",
    "agreement_status",
    "mismatch_type",
    "requires_manual_review",
    "prioritisation_recommendation",
]

FIGURE_2_LABEL_PAIRS = [
    "Dexamethasone-COVID-19",
    "Disulfiram-SARS-CoV-2",
    "BIO101-COVID-19",
    "Aspirin-Colorectal cancer",
    "Digoxin-Medulloblastoma",
    "Decitabine-K-RAS-dependent pancreatic ductal adenocarcinoma",
]

PAIR_ABBREVIATIONS = {
    "Dexamethasone-COVID-19": "Dex-COVID-19",
    "Remdesivir-COVID-19": "Remdesivir-COVID-19",
    "Thalidomide-Multiple myeloma": "Thalidomide-MM",
    "Sildenafil-Erectile dysfunction": "Sildenafil-ED",
    "Aspirin-Colorectal cancer": "Aspirin-CRC",
    "Sorafenib-Hepatocellular carcinoma": "Sorafenib-HCC",
    "Carfilzomib-SARS-CoV-2": "Carfilzomib-SARS-CoV-2",
    "Mibefradil-Glioma": "Mibefradil-glioma",
    "BIO101-COVID-19": "BIO101-COVID-19",
    "Minoxidil-Ovarian cancer": "Minoxidil-ovarian",
    "Cimetidine-Lung adenocarcinoma": "Cimetidine-LUAD",
    "Disulfiram-SARS-CoV-2": "Disulfiram-SARS-CoV-2",
    "Digoxin-Medulloblastoma": "Digoxin-MB",
    "Niclosamide-Colorectal cancer": "Niclosamide-CRC",
    "Decitabine-K-RAS-dependent pancreatic ductal adenocarcinoma": "Decitabine-KRAS PDAC",
}

FIGURE_3_GROUPS = [
    (
        "Established / trial-supported",
        [
            "Dexamethasone-COVID-19",
            "Remdesivir-COVID-19",
            "Thalidomide-Multiple myeloma",
            "Sildenafil-Erectile dysfunction",
        ],
    ),
    (
        "Safety-sensitive",
        [
            "Aspirin-Colorectal cancer",
            "Sorafenib-Hepatocellular carcinoma",
            "Carfilzomib-SARS-CoV-2",
            "Mibefradil-Glioma",
        ],
    ),
    (
        "Sparse / emerging",
        [
            "BIO101-COVID-19",
            "Minoxidil-Ovarian cancer",
            "Cimetidine-Lung adenocarcinoma",
        ],
    ),
    (
        "Computational / preclinical",
        [
            "Disulfiram-SARS-CoV-2",
            "Digoxin-Medulloblastoma",
            "Niclosamide-Colorectal cancer",
            "Decitabine-K-RAS-dependent pancreatic ductal adenocarcinoma",
        ],
    ),
]

THEME = {
    "therapeutic": "#118A7E",
    "adverse": "#D65F00",
    "noise": "#8A93A3",
    "uncertainty": "#6B5FA7",
    "manual_review": "#E6A700",
    "blue": "#2C7FB8",
    "green": "#2CA25F",
    "red": "#C84C31",
    "ink": "#202124",
    "muted": "#5F6368",
    "panel": "#F6F7F9",
}

SAFETY_STYLES = {
    "high": {
        "label": "Safety concern",
        "marker": "^",
        "color": THEME["adverse"],
    },
    "moderate": {
        "label": "Safety-aware / adjudication",
        "marker": "s",
        "color": THEME["manual_review"],
    },
    "low": {
        "label": "No major safety overlap",
        "marker": "o",
        "color": THEME["blue"],
    },
}

SAFETY_BURDEN_COLORS = {
    "high": THEME["adverse"],
    "moderate": THEME["manual_review"],
    "low": THEME["blue"],
}


def _set_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 130,
            "savefig.dpi": 300,
            "font.family": "DejaVu Sans",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "axes.labelsize": 11,
            "axes.titlesize": 14,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "legend.title_fontsize": 9,
        }
    )


def _read_benchmark(benchmark_dir: Path) -> pd.DataFrame:
    csv_path = benchmark_dir / "benchmark_pair_outputs_full.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing benchmark file: {csv_path}")

    df = pd.read_csv(csv_path)
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Benchmark CSV is missing required columns: {missing}")

    numeric_columns = [
        "literature_count",
        "therapeutic_rate",
        "adverse_rate",
        "irrelevant_rate",
        "gamma_safety_overlap",
        "posterior_mean",
        "posterior_ci_low",
        "posterior_ci_high",
        "credible_interval_width",
        "evidence_readiness_score",
    ]
    for column in numeric_columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    if "pair_id" in df.columns:
        df["pair_id"] = pd.to_numeric(df["pair_id"], errors="coerce")
        df = df.sort_values("pair_id", kind="stable")

    for column in REQUIRED_COLUMNS:
        if column not in numeric_columns:
            df[column] = df[column].fillna("").astype(str)

    return df.reset_index(drop=True)


def _write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def _save_figure(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _wrap(text: object, width: int = 28) -> str:
    value = str(text or "").replace("_", " ")
    return "\n".join(textwrap.wrap(value, width=width, break_long_words=False)) or "missing"


def _compact_label(text: object, width: int = 34) -> str:
    value = str(text or "")
    return "\n".join(textwrap.wrap(value, width=width, break_long_words=False)) or "missing"


def _format_number(value: object, decimals: int) -> str:
    try:
        if pd.isna(value):
            return ""
        return f"{float(value):.{decimals}f}"
    except (TypeError, ValueError):
        return ""


def _boolish(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _first_seen(values: Iterable[object]) -> list[str]:
    seen: list[str] = []
    for value in values:
        label = str(value or "")
        if label and label not in seen:
            seen.append(label)
    return seen


def _short_pair(pair: object) -> str:
    value = str(pair or "")
    return PAIR_ABBREVIATIONS.get(value, value)


def _safety_burden(row: pd.Series) -> str:
    flag = str(row.get("safety_quality_flag", "")).strip().lower()
    gamma = row.get("gamma_safety_overlap")
    if flag == "high_safety_overlap":
        return "high"
    if flag == "safety_aware":
        return "moderate"
    try:
        if not pd.isna(gamma):
            gamma_value = float(gamma)
            if gamma_value >= 0.70:
                return "high"
            if gamma_value >= 0.30:
                return "moderate"
    except (TypeError, ValueError):
        pass
    return "low"


def _maturity_overcalled(row: pd.Series) -> bool:
    expected = str(row.get("expected_pipeline_flag", "")).lower()
    observed = str(row.get("observed_quality_flag", "")).lower()
    immature_terms = [
        "preclinical",
        "computational",
        "validation_needed",
        "not_established",
        "sparse",
        "emerging",
        "translation_limited",
        "incomplete",
    ]
    return "high evidence quality" in observed and any(term in expected for term in immature_terms)


def _driver_counts(df: pd.DataFrame) -> dict[str, int]:
    mismatch = df["mismatch_type"].fillna("").astype(str).str.lower()
    expected = df["expected_pipeline_flag"].fillna("").astype(str).str.lower()
    observed = df["observed_quality_flag"].fillna("").astype(str).str.lower()

    overcalled = int(
        (mismatch == "overcalled_maturity").sum()
        + df.apply(_maturity_overcalled, axis=1).sum()
    )
    missed_safety = int(
        (mismatch == "missed_safety_signal").sum()
        + ((expected.str.contains("safety")) & (~observed.str.contains("safety"))).sum()
    )
    translation_limited = (
        int(df["translation_limitation_flag"].map(_boolish).sum())
        if "translation_limitation_flag" in df.columns
        else int(
            df["translation_evidence_level"]
            .isin(["preclinical", "computational_only", "combination_context", "screening_supported"])
            .sum()
        )
    )
    context_issue = int(
        (
            (mismatch == "other")
            | expected.str.contains("context|reviewable|uncertainty|incomplete|not_established")
        ).sum()
    )
    manual_review = int(df["requires_manual_review"].map(_boolish).sum())
    return {
        "Manual review retained": manual_review,
        "Translation limitation": translation_limited,
        "Context-specific evidence gap": context_issue,
        "Overcalled maturity": overcalled,
        "Missed safety burden": missed_safety,
    }


def _create_output_dir(output_dir: Path | None) -> Path:
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = Path("outputs") / f"drug_discovery_today_figures_tables_{stamp}"
    out.mkdir(parents=True, exist_ok=False)
    return out


def make_table_1(df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    columns = [
        "drug_disease_pair",
        "source_panel",
        "external_evidence_maturity",
        "expected_pipeline_flag",
    ]
    table = df[columns].copy()
    _write_csv(table, output_dir / "Table_1_benchmark_design.csv")
    return table


def make_table_2(df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    table = df[
        [
            "drug_disease_pair",
            "coverage_tier",
            "literature_count",
            "therapeutic_rate",
            "adverse_rate",
            "gamma_safety_overlap",
            "posterior_mean",
            "credible_interval_width",
            "evidence_readiness_score",
            "translation_evidence_level",
            "observed_quality_flag",
            "prioritisation_recommendation",
        ]
    ].copy()
    table = table.sort_values("evidence_readiness_score", ascending=False, kind="stable")
    table["literature_count"] = table["literature_count"].fillna(0).round(0).astype(int)
    table["therapeutic_rate"] = table["therapeutic_rate"].map(lambda value: _format_number(value, 2))
    table["adverse_rate"] = table["adverse_rate"].map(lambda value: _format_number(value, 2))
    table["gamma_safety_overlap"] = table["gamma_safety_overlap"].map(
        lambda value: _format_number(value, 2)
    )
    table["posterior_mean"] = table["posterior_mean"].map(lambda value: _format_number(value, 3))
    table["credible_interval_width"] = table["credible_interval_width"].map(
        lambda value: _format_number(value, 3)
    )
    table["evidence_readiness_score"] = table["evidence_readiness_score"].map(
        lambda value: _format_number(value, 1)
    )
    _write_csv(table, output_dir / "Table_2_core_evidence_readiness_results.csv")
    return table


def make_table_3(df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    table = df[
        [
            "drug_disease_pair",
            "expected_pipeline_flag",
            "observed_quality_flag",
            "agreement_status",
            "mismatch_type",
            "requires_manual_review",
        ]
    ].copy()
    table["requires_manual_review"] = table["requires_manual_review"].map(_boolish)
    agreement_order = {"mismatch": 0, "partial_match": 1, "match": 2}
    table["_agreement_sort"] = table["agreement_status"].map(agreement_order).fillna(99)
    table["_review_sort"] = table["requires_manual_review"].astype(int)
    table = table.sort_values(
        ["_agreement_sort", "_review_sort", "drug_disease_pair"],
        ascending=[True, False, True],
        kind="stable",
    ).drop(columns=["_agreement_sort", "_review_sort"])
    _write_csv(table, output_dir / "Table_3_external_adjudication.csv")
    return table


def make_figure_1(output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(14.2, 6.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    zones = [
        {
            "x": 0.04,
            "w": 0.25,
            "title": "Candidate signal",
            "subtitle": "Computational candidate generation",
            "items": [
                "drug-disease pair",
                "literature",
                "trial records",
                "safety data",
                "network/omics data",
            ],
            "face": "#F6F7F9",
            "edge": "#7A869A",
        },
        {
            "x": 0.365,
            "w": 0.27,
            "title": "Evidence-readiness audit",
            "subtitle": "Decision-quality layer",
            "items": [
                "Terminology standardisation",
                "Literature composition",
                "Safety burden",
                "External support",
                "Bayesian uncertainty and maturity",
            ],
            "face": "#EDF7F5",
            "edge": THEME["therapeutic"],
        },
        {
            "x": 0.71,
            "w": 0.25,
            "title": "Decision-ready interpretation",
            "subtitle": "Prioritisation with safeguards",
            "items": [
                "prioritise for review",
                "safety adjudication",
                "preclinical validation",
                "watchlist",
                "hold back",
            ],
            "face": "#FFF7E0",
            "edge": THEME["manual_review"],
        },
    ]

    for zone in zones:
        x = zone["x"]
        w = zone["w"]
        box = FancyBboxPatch(
            (x, 0.24),
            w,
            0.50,
            boxstyle="round,pad=0.018,rounding_size=0.025",
            linewidth=1.5,
            edgecolor=zone["edge"],
            facecolor=zone["face"],
        )
        ax.add_patch(box)
        ax.text(
            x + w / 2,
            0.68,
            zone["title"],
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
            color=THEME["ink"],
        )
        ax.text(
            x + w / 2,
            0.62,
            zone["subtitle"],
            ha="center",
            va="center",
            fontsize=9,
            color=THEME["muted"],
        )
        for index, item in enumerate(zone["items"]):
            y = 0.55 - index * 0.065
            tile = FancyBboxPatch(
                (x + 0.03, y - 0.025),
                w - 0.06,
                0.043,
                boxstyle="round,pad=0.007,rounding_size=0.012",
                linewidth=0.8,
                edgecolor="#FFFFFF",
                facecolor="#FFFFFF",
                alpha=0.95,
            )
            ax.add_patch(tile)
            ax.text(
                x + w / 2,
                y - 0.003,
                item,
                ha="center",
                va="center",
                fontsize=8.8,
                color=THEME["ink"],
            )

    for start, end in [((0.305, 0.49), (0.345, 0.49)), ((0.655, 0.49), (0.695, 0.49))]:
        ax.annotate(
            "",
            xy=end,
            xytext=start,
            arrowprops=dict(arrowstyle="-|>", lw=2.0, color=THEME["ink"]),
        )

    ax.text(
        0.5,
        0.10,
        "The audit evaluates evidence quality before prioritisation. It does not infer clinical efficacy.",
        ha="center",
        va="center",
        fontsize=10.5,
        color=THEME["ink"],
    )
    _save_figure(fig, output_dir / "Figure_1_workflow.png")


def make_figure_2(df: pd.DataFrame, output_dir: Path) -> None:
    plot_df = df.copy()
    plot_df["_safety_burden"] = plot_df.apply(_safety_burden, axis=1)

    max_literature = max(float(plot_df["literature_count"].max()), 1.0)
    plot_df["_point_size"] = 70 + 360 * (
        plot_df["literature_count"].fillna(0).map(lambda value: math.sqrt(max(value, 0)))
        / math.sqrt(max_literature)
    )

    median_score = float(plot_df["evidence_readiness_score"].median())
    median_width = float(plot_df["credible_interval_width"].median())
    x_min = max(0.0, float(plot_df["evidence_readiness_score"].min()) - 4)
    x_max = min(100.0, float(plot_df["evidence_readiness_score"].max()) + 5)
    y_min = max(0.0, float(plot_df["credible_interval_width"].min()) - 0.025)
    y_max = float(plot_df["credible_interval_width"].max()) + 0.065

    fig, ax = plt.subplots(figsize=(11.2, 7.0))
    zones = [
        (x_min, median_score, median_width, y_max, "#F1EFF8", "Immature\nor sparse"),
        (median_score, x_max, median_width, y_max, "#FFF3CF", "Promising\nbut unstable"),
        (x_min, median_score, y_min, median_width, "#F0F2F4", "Consistently\nweak"),
        (median_score, x_max, y_min, median_width, "#E7F4F1", "Reviewable\nsignal"),
    ]
    for x0, x1, y0, y1, color, label in zones:
        ax.add_patch(
            Rectangle(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                facecolor=color,
                edgecolor="none",
                zorder=0,
            )
        )
        if label == "Reviewable\nsignal":
            label_x = x0 + (x1 - x0) * 0.66
            label_y = y0 + (y1 - y0) * 0.25
        elif label == "Consistently\nweak":
            label_x = x0 + (x1 - x0) * 0.06
            label_y = y0 + (y1 - y0) * 0.25
        else:
            label_x = x0 + (x1 - x0) * 0.06
            label_y = y1 - (y1 - y0) * 0.16
        ax.text(
            label_x,
            label_y,
            label,
            ha="left",
            va="center" if label in {"Reviewable\nsignal", "Consistently\nweak"} else "top",
            fontsize=11,
            fontweight="bold",
            color=THEME["ink"],
            alpha=0.86,
        )

    for burden, style in SAFETY_STYLES.items():
        subset = plot_df[plot_df["_safety_burden"] == burden]
        if subset.empty:
            continue
        ax.scatter(
            subset["evidence_readiness_score"],
            subset["credible_interval_width"],
            s=subset["_point_size"],
            c=style["color"],
            marker=style["marker"],
            edgecolors=THEME["ink"],
            linewidths=0.85,
            alpha=0.86,
            label=style["label"],
            zorder=3,
        )

    ax.axvline(median_score, color=THEME["ink"], lw=1.0, ls="--", alpha=0.65)
    ax.axhline(median_width, color=THEME["ink"], lw=1.0, ls="--", alpha=0.65)

    label_offsets = {
        "Dexamethasone-COVID-19": (18, 20),
        "Disulfiram-SARS-CoV-2": (18, -24),
        "BIO101-COVID-19": (12, 14),
        "Aspirin-Colorectal cancer": (-110, 20),
        "Digoxin-Medulloblastoma": (-116, -30),
        "Decitabine-K-RAS-dependent pancreatic ductal adenocarcinoma": (14, 12),
    }
    for pair in FIGURE_2_LABEL_PAIRS:
        match = plot_df[plot_df["drug_disease_pair"] == pair]
        if match.empty:
            continue
        row = match.iloc[0]
        dx, dy = label_offsets.get(pair, (8, 8))
        ax.annotate(
            _short_pair(pair),
            xy=(row["evidence_readiness_score"], row["credible_interval_width"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=8.2,
            arrowprops=dict(arrowstyle="-", color=THEME["muted"], lw=0.65),
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.88),
            zorder=4,
        )

    ax.set_xlabel("Evidence-readiness score")
    ax.set_ylabel("Posterior uncertainty (credible interval width)")
    ax.grid(axis="both", color="#FFFFFF", lw=1.0, alpha=0.95)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    safety_handles = [
        Line2D(
            [0],
            [0],
            marker=style["marker"],
            linestyle="",
            markerfacecolor=style["color"],
            markeredgecolor=THEME["ink"],
            label=style["label"],
            markersize=8,
        )
        for burden, style in SAFETY_STYLES.items()
        if burden in set(plot_df["_safety_burden"])
    ]
    ax.legend(
        handles=safety_handles,
        title="Safety status",
        loc="upper right",
        frameon=True,
        facecolor="white",
        edgecolor="#D0D4D9",
    )
    _save_figure(fig, output_dir / "Figure_2_readiness_uncertainty_landscape.png")


def make_figure_3(df: pd.DataFrame, output_dir: Path) -> None:
    records = []
    order = 0
    for group, pairs in FIGURE_3_GROUPS:
        for pair in pairs:
            match = df[df["drug_disease_pair"] == pair]
            if match.empty:
                continue
            row = match.iloc[0].copy()
            row["_group"] = group
            row["_order"] = order
            row["_safety_burden"] = _safety_burden(row)
            records.append(row)
            order += 1
    plot_df = pd.DataFrame(records).sort_values("_order").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(11.8, 7.4))
    fig.subplots_adjust(left=0.28, right=0.88, bottom=0.15, top=0.88)
    n_rows = len(plot_df)
    y_positions = range(n_rows)

    therapeutic = plot_df["therapeutic_rate"].fillna(0)
    adverse = plot_df["adverse_rate"].fillna(0)
    irrelevant = plot_df["irrelevant_rate"].fillna(0)

    ax.barh(y_positions, therapeutic, color=THEME["therapeutic"], label="Therapeutic")
    ax.barh(y_positions, adverse, left=therapeutic, color=THEME["adverse"], label="Adverse")
    ax.barh(
        y_positions,
        irrelevant,
        left=therapeutic + adverse,
        color=THEME["noise"],
        label="Irrelevant / noisy",
    )

    for y, row in plot_df.iterrows():
        burden = row["_safety_burden"]
        ax.add_patch(
            Rectangle(
                (-0.07, y - 0.35),
                0.035,
                0.70,
                facecolor=SAFETY_BURDEN_COLORS[burden],
                edgecolor="white",
                lw=0.8,
                clip_on=False,
            )
        )
        literature_count = int(row["literature_count"]) if not pd.isna(row["literature_count"]) else 0
        ax.text(
            1.015,
            y,
            f"n={literature_count}",
            va="center",
            ha="left",
            fontsize=8.2,
            color=THEME["muted"],
        )

    for _, group in plot_df.groupby("_group", sort=False):
        first_index = int(group.index.min())
        ax.axhline(first_index - 0.52, color="#D9DDE3", lw=0.9)
        ax.text(
            -0.09,
            first_index - 0.42,
            str(group["_group"].iloc[0]),
            va="bottom",
            ha="right",
            fontsize=8.3,
            color=THEME["muted"],
            fontweight="bold",
            transform=ax.get_yaxis_transform(),
        )

    ax.set_yticks(list(y_positions))
    ax.set_yticklabels([_short_pair(pair) for pair in plot_df["drug_disease_pair"]])
    ax.invert_yaxis()
    ax.set_xlim(-0.08, 1.10)
    ax.set_xlabel("Share of classified literature")
    ax.grid(axis="x", color="#D8D8D8", lw=0.7, alpha=0.65)
    ax.text(-0.073, -0.95, "Safety\nburden", ha="left", va="center", fontsize=8.0, color=THEME["muted"])
    safety_handles = [
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            markerfacecolor=color,
            markeredgecolor=color,
            label=label,
            markersize=7,
        )
        for label, color in [
            ("low", SAFETY_BURDEN_COLORS["low"]),
            ("moderate", SAFETY_BURDEN_COLORS["moderate"]),
            ("high", SAFETY_BURDEN_COLORS["high"]),
        ]
    ]
    ax.legend(
        handles=ax.get_legend_handles_labels()[0] + safety_handles,
        labels=ax.get_legend_handles_labels()[1] + ["low safety", "moderate safety", "high safety"],
        loc="lower center",
        bbox_to_anchor=(0.52, -0.16),
        ncol=6,
        frameon=False,
        fontsize=8,
    )
    _save_figure(fig, output_dir / "Figure_3_literature_composition_safety.png")


def make_figure_4(df: pd.DataFrame, output_dir: Path) -> None:
    counts = df["agreement_status"].value_counts()
    matches = int(counts.get("match", 0))
    partial = int(counts.get("partial_match", 0))
    mismatches = int(counts.get("mismatch", 0))
    total = max(len(df), 1)
    agreement_rate = (matches + 0.5 * partial) / total

    fig, (ax_a, ax_b) = plt.subplots(
        1,
        2,
        figsize=(12.4, 5.8),
        gridspec_kw={"width_ratios": [1.05, 1.25], "wspace": 0.34},
    )
    agreement_labels = ["Match", "Partial match", "Mismatch"]
    agreement_values = [matches, partial, mismatches]
    agreement_colors = [THEME["therapeutic"], THEME["manual_review"], THEME["adverse"]]
    left = 0
    for label, value, color in zip(agreement_labels, agreement_values, agreement_colors):
        ax_a.barh([0], [value], left=left, height=0.36, color=color, label=label)
        ax_a.text(
            left + value / 2,
            0,
            f"{value}",
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
            color="white" if value >= 3 else THEME["ink"],
        )
        left += value
    ax_a.set_xlim(0, total)
    ax_a.set_ylim(-0.65, 0.85)
    ax_a.set_yticks([])
    ax_a.set_xlabel("Benchmark pairs")
    ax_a.set_title("A  Agreement summary", loc="left", pad=12)
    ax_a.legend(loc="lower center", bbox_to_anchor=(0.5, -0.28), ncol=3, frameon=False)
    ax_a.text(
        total / 2,
        0.54,
        f"Weighted agreement = {agreement_rate:.3f}",
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
        color=THEME["ink"],
    )
    drivers = _driver_counts(df)
    driver_series = pd.Series(drivers).sort_values(ascending=True)
    colors = [
        THEME["manual_review"] if "Manual" in label else THEME["uncertainty"]
        for label in driver_series.index
    ]
    ax_b.barh(driver_series.index, driver_series.values, color=colors, height=0.58)
    for y, value in enumerate(driver_series.values):
        ax_b.text(value + 0.5, y, str(int(value)), va="center", fontsize=9, color=THEME["ink"])
    ax_b.set_xlim(0, max(driver_series.max() + 4, 10))
    ax_b.set_xlabel("Pairs flagged")
    ax_b.set_title("B  Main triage signals", loc="left", pad=12)
    ax_b.grid(axis="x", color="#D8D8D8", lw=0.7, alpha=0.7)
    _save_figure(fig, output_dir / "Figure_4_expected_vs_observed_maturity.png")


def make_graphical_abstract(output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.8, 6.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    boxes = [
        (0.22, 0.74, 0.56, 0.12, "Computational candidate", "#F6F7F9", "#7A869A"),
        (0.17, 0.43, 0.66, 0.18, "Evidence-readiness audit", "#EDF7F5", THEME["therapeutic"]),
        (0.22, 0.17, 0.56, 0.12, "Decision-ready interpretation", "#FFF7E0", THEME["manual_review"]),
    ]
    for x, y, w, h, label, face, edge in boxes:
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.018,rounding_size=0.026",
                facecolor=face,
                edgecolor=edge,
                linewidth=1.7,
            )
        )
        ax.text(
            x + w / 2,
            y + h - 0.045 if label == "Evidence-readiness audit" else y + h / 2,
            label,
            ha="center",
            va="center",
            fontsize=13 if label == "Evidence-readiness audit" else 14,
            fontweight="bold",
            color=THEME["ink"],
        )

    for start, end in [((0.50, 0.725), (0.50, 0.63)), ((0.50, 0.405), (0.50, 0.31))]:
        ax.annotate(
            "",
            xy=end,
            xytext=start,
            arrowprops=dict(arrowstyle="-|>", lw=2.0, color=THEME["ink"]),
        )

    audit_labels = [
        ("terminology", THEME["blue"]),
        ("literature\ncomposition", THEME["therapeutic"]),
        ("safety\nburden", THEME["adverse"]),
        ("external\nsupport", THEME["manual_review"]),
        ("uncertainty /\nmaturity", THEME["uncertainty"]),
    ]
    x_positions = [0.24, 0.37, 0.50, 0.63, 0.76]
    for x, (label, color) in zip(x_positions, audit_labels):
        ax.add_patch(
            FancyBboxPatch(
                (x - 0.055, 0.47),
                0.11,
                0.07,
                boxstyle="round,pad=0.01,rounding_size=0.018",
                facecolor="white",
                edgecolor=color,
                linewidth=1.2,
            )
        )
        ax.text(x, 0.505, label, ha="center", va="center", fontsize=8.4, color=THEME["ink"])

    ax.text(
        0.5,
        0.08,
        "From ranked candidates to interpretable evidence profiles",
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
        color=THEME["ink"],
    )
    _save_figure(fig, output_dir / "Graphical_Abstract.png")


def write_figure_captions(output_dir: Path) -> None:
    text = """# Figure Captions

**Figure 1. Evidence-readiness audit workflow.** The diagram shows the sequential audit from drug-disease input through terminology standardisation, literature classification, safety review, graph support, Bayesian uncertainty, and prioritisation. It matters because the framework evaluates evidence quality and traceability, not clinical efficacy. This supports the data-quality-first argument by placing curation, safety, and uncertainty before candidate prioritisation.

**Figure 2. Evidence-readiness versus uncertainty landscape.** The scatter plot shows evidence-readiness score against posterior credible interval width, with point size scaled by literature volume, colour indicating translation evidence level, and marker style indicating safety status. It matters because high readiness and low uncertainty are distinct properties. This supports the data-quality-first argument by highlighting candidates that need safety review, uncertainty review, or manual triage before follow-up.

**Figure 3. Literature composition and safety burden.** The stacked horizontal bars show the therapeutic, adverse, and irrelevant shares of classified literature for each benchmark pair, annotated by literature count, safety overlap, and safety flag. It matters because literature volume alone can hide adverse or noisy evidence. This supports the data-quality-first argument by separating supportive evidence from safety burden and irrelevant retrieval.

**Figure 4. Expected versus observed evidence maturity.** The matrix compares expected pipeline flags with observed quality flags and reports match, partial-match, mismatch, and weighted agreement summaries. It matters because an evidence-readiness audit should avoid overcalling immature, sparse, computational, or preclinical examples. This supports the data-quality-first argument by making disagreement and manual-review cases visible.
"""
    (output_dir / "figure_captions.md").write_text(text, encoding="utf-8")


def write_table_captions(output_dir: Path) -> None:
    text = """# Table Captions

**Table 1. Benchmark design.** This table lists the selected drug-disease pairs, benchmark panel, external evidence maturity, and expected pipeline flag. It matters because the benchmark intentionally spans established, safety-sensitive, sparse, computational, preclinical, and translation-limited examples. This supports the data-quality-first argument by showing that the audit was tested across diverse evidence states.

**Table 2. Core evidence-readiness results.** This table reports concise quantitative outputs for each pair, including coverage, literature composition, safety overlap, posterior mean, uncertainty width, evidence-readiness score, translation level, quality flag, and prioritisation recommendation. It matters because the same score is interpreted alongside safety, uncertainty, and translation maturity. This supports the data-quality-first argument by keeping prioritisation linked to evidence limitations.

**Table 3. External adjudication and manual-review cases.** This table compares expected pipeline flags with observed quality flags and highlights match status, mismatch type, and manual-review requirements. It matters because disagreements reveal overcalled maturity, missed safety signals, translation limitations, or cases needing expert judgement. This supports the data-quality-first argument by treating mismatch as information for curation rather than as a simple failure.
"""
    (output_dir / "table_captions.md").write_text(text, encoding="utf-8")


def write_results_summary(df: pd.DataFrame, output_dir: Path) -> None:
    panel_counts = df["source_panel"].value_counts(sort=False)
    panel_summary = "; ".join(f"{panel}: {count}" for panel, count in panel_counts.items())
    total_pairs = len(df)
    full_audit = int((df["coverage_tier"] == "full_bayesian_audit").sum())
    graph_only = int((df["coverage_tier"] == "graph_only").sum())
    safety_quality_count = int((df["observed_quality_flag"] == "Safety-concerning").sum())
    safety_signal_count = int(
        df["safety_quality_flag"].isin(["high_safety_overlap", "safety_aware"]).sum()
    )
    if "translation_limitation_flag" in df.columns:
        translation_limited = int(df["translation_limitation_flag"].map(_boolish).sum())
    else:
        translation_limited = int(
            df["translation_evidence_level"]
            .isin(["preclinical", "computational_only", "combination_context", "screening_supported"])
            .sum()
        )

    counts = df["agreement_status"].value_counts()
    matches = int(counts.get("match", 0))
    partial = int(counts.get("partial_match", 0))
    mismatches = int(counts.get("mismatch", 0))
    agreement_rate = (matches + 0.5 * partial) / max(total_pairs, 1)
    manual_review = int(df["requires_manual_review"].map(_boolish).sum())

    paragraphs = [
        (
            f"The benchmark included {total_pairs} drug-disease pairs across {len(panel_counts)} "
            f"source panels ({panel_summary}). The selected pairs were included to cover "
            "established clinical successes, COVID-19 and SARS-CoV-2 examples, oncology "
            "repurposing cases, safety-sensitive examples, sparse or emerging candidates, "
            "preclinical and computational signals, and translation-limited contexts."
        ),
        (
            f"The pipeline processed all {total_pairs} pairs, with {full_audit} reaching full "
            f"Bayesian audit coverage and {graph_only} retained as graph-only evidence. The "
            f"outputs separated evidence volume from evidence quality: {safety_signal_count} "
            "pairs had safety-aware or high-overlap safety signals, and "
            f"{safety_quality_count} received the Safety-concerning observed quality flag. "
            f"{translation_limited} pairs were marked by translation limitations. The pipeline "
            "audits evidence readiness; it does not predict clinical efficacy and does not "
            "replace experimental or clinical validation."
        ),
        (
            f"External adjudication showed {matches} matches, {partial} partial matches, and "
            f"{mismatches} mismatches, giving a weighted agreement rate of {agreement_rate:.3f}. "
            f"{manual_review} pairs were marked for manual review, demonstrating that the output "
            "supports prioritisation, manual review, and evidence-quality triage. These findings "
            "argue for data-quality-first review before computational repurposing candidates are "
            "advanced to experimental or clinical follow-up."
        ),
    ]
    (output_dir / "results_summary_for_manuscript.md").write_text(
        "\n\n".join(paragraphs) + "\n",
        encoding="utf-8",
    )


def generate_package(
    benchmark_dir: Path,
    output_dir: Path | None = None,
    *,
    figures_only: bool = False,
) -> Path:
    _set_style()
    out = _create_output_dir(output_dir)
    df = _read_benchmark(benchmark_dir)

    if not figures_only:
        make_table_1(df, out)
        make_table_2(df, out)
        make_table_3(df, out)

    make_figure_1(out)
    make_figure_2(df, out)
    make_figure_3(df, out)
    make_figure_4(df, out)
    make_graphical_abstract(out)

    if not figures_only:
        write_figure_captions(out)
        write_table_captions(out)
        write_results_summary(df, out)
    return out


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate focused manuscript figures and tables from benchmark outputs.",
    )
    parser.add_argument(
        "--benchmark-dir",
        default="outputs/evidence_readiness_benchmark_20260618_204625",
        help="Directory containing benchmark_pair_outputs_full.csv.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional output directory. Defaults to outputs/drug_discovery_today_figures_tables_<timestamp>/.",
    )
    parser.add_argument(
        "--figures-only",
        action="store_true",
        help="Only regenerate figure PNGs and the graphical abstract.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    out = generate_package(
        benchmark_dir=Path(args.benchmark_dir),
        output_dir=Path(args.output_dir) if args.output_dir else None,
        figures_only=args.figures_only,
    )
    files = [
        "Graphical_Abstract.png",
        "Figure_1_workflow.png",
        "Figure_2_readiness_uncertainty_landscape.png",
        "Figure_3_literature_composition_safety.png",
        "Figure_4_expected_vs_observed_maturity.png",
    ]
    if not args.figures_only:
        files.extend(
            [
                "Table_1_benchmark_design.csv",
                "Table_2_core_evidence_readiness_results.csv",
                "Table_3_external_adjudication.csv",
                "figure_captions.md",
                "table_captions.md",
                "results_summary_for_manuscript.md",
            ]
        )
    print(f"Generated manuscript figure/table package: {out}")
    for file_name in files:
        print(f"  {file_name}")


if __name__ == "__main__":
    main()
