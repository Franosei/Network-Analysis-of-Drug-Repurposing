"""
Figure 1 (case study): "What disappears when different evidence sources become one row."

Minimal three-part editorial figure from real pipeline outputs for one drug-disease
pair (hydroxychloroquine / covid-19):
  (A) four evidence sources that measure different things for the same pair,
  (B) a small number of meaningful integration / judgement steps,
  (C) the final, conflicted evidence profile for this pair.

All numeric content is read at run time from the pipeline's own audit artifacts:
  outputs/20260923_hcq_covid_full_pipeline/manuscript_tables/pair_level_evidence_quality.csv
  outputs/20260923_hcq_covid_full_pipeline/ledgers/full_evidence_quality_ledger.csv

Run provenance (run id, timestamp, source file paths) is intentionally kept out of
the image and belongs in the figure caption / source note instead.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = PROJECT_ROOT / "outputs/20260923_hcq_covid_full_pipeline"
PAIR_CSV = RUN_ROOT / "manuscript_tables/pair_level_evidence_quality.csv"
LEDGER_CSV = RUN_ROOT / "ledgers/full_evidence_quality_ledger.csv"
OUT_DIR = RUN_ROOT / "manuscript_figures"

DRUG = "hydroxychloroquine"
DISEASE = "covid-19"

# Palette consistent with visualisation/make_all_publication_figures.py
C_TRIALS = "#4C72B0"     # blue
C_LIT = "#55A868"        # green
C_SAFETY = "#C44E52"     # red
C_NETWORK = "#8172B2"    # purple
C_COMPOSITE = "#946E3D"  # amber/bronze for the composite/final row
C_MIXED = "#C9A227"      # amber for "mixed / moderate" profile rows
C_INK = "#1F1F1F"
C_MUTED = "#5B5B5B"
C_NEUTRAL_BORDER = "#B7A98C"
C_NEUTRAL_FILL = "#F7F5F2"


def load_data() -> dict:
    pair_df = pd.read_csv(PAIR_CSV)
    row = pair_df[
        (pair_df["drug"].str.lower() == DRUG) & (pair_df["disease"].str.lower() == DISEASE)
    ].iloc[0]

    ledger_df = pd.read_csv(LEDGER_CSV)
    lrow = ledger_df[
        (ledger_df["drug"].str.lower() == DRUG) & (ledger_df["disease"].str.lower() == DISEASE)
    ].iloc[0]
    run_path = next(RUN_ROOT.glob("runs/run_hydroxychloroquine_covid-19_*.json"))
    run_payload = json.loads(run_path.read_text(encoding="utf-8"))
    graph_df = pd.read_csv(RUN_ROOT / "graph/graph_features_known.csv")
    graph_row = graph_df[
        (graph_df["Drug"].str.lower() == DRUG) & (graph_df["Disease"].str.lower() == DISEASE)
    ].iloc[0]

    return {
        "trial_count": int(lrow["trial_count"]),
        "records_scanned": int(run_payload["components"]["records_retrieved"]),
        "records_retrieved": int(row["records_retrieved"]),
        "therapeutic": int(row["therapeutic_count"]),
        "adverse": int(row["adverse_count"]),
        "irrelevant": int(row["irrelevant_count"]),
        "therapeutic_rate": float(row["therapeutic_relevance_ratio"]),
        "gamma": float(row["gamma_safety_overlap"]),
        "n_safety_terms": int(row["safety_overlap_term_count"]),
        "structural_consistency": float(row["structural_consistency_score"]),
        "graph_probability": float(graph_row["GraphProbability"]),
        "alternative_paths": int(graph_row["AlternativePathCountLength3"]),
        "entity_mapping": float(row["entity_mapping_quality_score"]),
        "posterior_mean": float(row["posterior_mean"]),
        "ci_low": float(row["posterior_ci_low"]),
        "ci_high": float(row["posterior_ci_high"]),
        "ci_width": float(row["credible_interval_width"]),
        "evidence_readiness_score": float(row["evidence_readiness_score"]),
        "quality_flag": str(row["quality_flag"]),
    }


def rbox(ax, x, y, w, h, fc, ec, lw=1.1, alpha=1.0, rounding=0.02, z=2):
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={rounding}",
        linewidth=lw, edgecolor=ec, facecolor=fc, alpha=alpha, zorder=z,
        mutation_aspect=1,
    )
    ax.add_patch(box)
    return box


def down_arrow(ax, x, y0, y1, color=C_MUTED, lw=1.3):
    ax.add_patch(FancyArrowPatch(
        (x, y0), (x, y1),
        arrowstyle="-|>", mutation_scale=10, linewidth=lw,
        color=color, zorder=3, shrinkA=0, shrinkB=0,
    ))


def right_arrow(ax, x0, x1, y, color=C_NEUTRAL_BORDER, lw=1.0):
    ax.add_patch(FancyArrowPatch(
        (x0, y), (x1, y),
        arrowstyle="-|>", mutation_scale=9, linewidth=lw,
        color=color, zorder=3, shrinkA=0, shrinkB=0,
    ))


def build_figure(d: dict, out_dir: Path) -> None:
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig = plt.figure(figsize=(16.5, 21.0), dpi=300)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 137)
    ax.invert_yaxis()
    ax.axis("off")

    # ---------------------------------------------------------------- title
    ax.text(50, 4.0, "What disappears when different evidence sources become one row",
            ha="center", va="center", fontsize=18.5, fontweight="bold", color=C_INK)
    ax.text(50, 7.6,
            "Hydroxychloroquine → COVID-19 — one drug–disease pair, four evidence sources",
            ha="center", va="center", fontsize=11.5, color=C_MUTED, style="italic")

    # =================================================================
    # PANEL A — four sources
    # =================================================================
    ay0 = 13.0
    ax.text(2, ay0, "A", fontsize=15, fontweight="bold", color=C_INK, va="top")
    ax.text(6.5, ay0, "Four sources, four different meanings",
            fontsize=12.5, fontweight="bold", color=C_INK, va="top")

    box_top = 17.0
    box_h = 16.0
    box_w = 22.3
    gap = 2.27
    xs = [2 + i * (box_w + gap) for i in range(4)]

    sources = [
        dict(
            color=C_TRIALS, title="ClinicalTrials.gov",
            big=f"{d['trial_count']} canonical trial records",
            extra=None,
            sentence="Measures investigation activity — not treatment effectiveness.",
        ),
        dict(
            color=C_LIT, title="PubMed / PMC",
            big=f"{d['records_scanned']:,} retrieved records",
            extra=f"{d['records_retrieved']:,} exact-pair records",
            sentence="Publication volume reflects attention — not therapeutic support.",
        ),
        dict(
            color=C_SAFETY, title="openFDA (FAERS)",
            big=f"γ = {d['gamma']:.2f} safety overlap",
            extra=f"{d['n_safety_terms']} overlapping terms",
            sentence="Reported events — not a causal treatment effect.",
        ),
        dict(
            color=C_NETWORK, title="Drug–disease network",
            big=f"graph probability {d['graph_probability']:.3f}",
            extra=f"{d['alternative_paths']:,} alternative length-3 paths",
            sentence="Structural proximity — not clinical efficacy.",
        ),
    ]

    for x, s in zip(xs, sources):
        rbox(ax, x, box_top, box_w, box_h, fc="white", ec=s["color"], lw=1.6, rounding=0.6)
        rbox(ax, x, box_top, box_w, 4.4, fc=s["color"], ec=s["color"], rounding=0.6, z=3)
        ax.text(x + box_w / 2, box_top + 2.2, s["title"], ha="center", va="center",
                fontsize=10.6, fontweight="bold", color="white", zorder=4)

        yy = box_top + 8.6
        ax.text(x + box_w / 2, yy, s["big"], ha="center", va="center",
                fontsize=12.4, fontweight="bold", color=C_INK)
        if s["extra"]:
            ax.text(x + box_w / 2, yy + 3.1, s["extra"], ha="center", va="center",
                    fontsize=7.6, color=C_MUTED)

        ax.plot([x + 1.2, x + box_w - 1.2], [box_top + box_h - 5.6, box_top + box_h - 5.6],
                color="#E3E3E3", lw=0.8)
        ax.text(x + box_w / 2, box_top + box_h - 2.9, s["sentence"], ha="center", va="center",
                fontsize=7.9, color=s["color"], style="italic", fontweight="bold")

    a_bottom = box_top + box_h
    for x in xs:
        down_arrow(ax, x + box_w / 2, a_bottom + 0.5, a_bottom + 3.7)

    # =================================================================
    # PANEL B — integration
    # =================================================================
    by0 = a_bottom + 6.6
    ax.text(2, by0, "B", fontsize=15, fontweight="bold", color=C_INK, va="top")
    ax.text(6.5, by0, "Integration: a few judgement points, not a simple join",
            fontsize=12.5, fontweight="bold", color=C_INK, va="top")
    ax.text(6.5, by0 + 3.1,
            "Each step below is a place where meaning is decided, not just looked up.",
            fontsize=8.8, color=C_MUTED, style="italic")

    node_top = by0 + 6.8
    node_h = 16.5
    node_w = 22.3
    node_gap = 2.27
    nxs = [2 + i * (node_w + node_gap) for i in range(4)]

    steps = [
        dict(title="Entity resolution", body=[
            "hydroxychloroquine / HCQ /",
            "hydroxychloroquine sulfate",
            "→ one drug entity",
        ]),
        dict(title="Disease harmonisation", body=[
            "COVID-19 / SARS-CoV-2 infection /",
            "coronavirus disease 2019",
            "→ one disease entity",
        ]),
        dict(title="Evidence classification", body=[
            f"{d['therapeutic']} therapeutic / {d['adverse']} adverse /",
            f"{d['irrelevant']} irrelevant",
        ]),
        dict(title="Safety-conflict flag", body=[
            "adverse-event overlap conflicts",
            "with the therapeutic literature",
            "→ flagged for review",
        ]),
    ]

    for x, s in zip(nxs, steps):
        rbox(ax, x, node_top, node_w, node_h, fc=C_NEUTRAL_FILL, ec=C_NEUTRAL_BORDER, lw=1.1, rounding=0.5)
        ax.text(x + node_w / 2, node_top + 2.8, s["title"], ha="center", va="center",
                fontsize=9.0, fontweight="bold", color=C_INK)
        ax.plot([x + 1.2, x + node_w - 1.2], [node_top + 4.6, node_top + 4.6], color="#D8CDB4", lw=0.7)
        yy = node_top + 7.2
        for ln in s["body"]:
            ax.text(x + node_w / 2, yy, ln, ha="center", va="center", fontsize=7.6, color=C_MUTED)
            yy += 2.35

    for i in range(3):
        right_arrow(ax, nxs[i] + node_w, nxs[i + 1], node_top + node_h / 2)

    node_bottom = node_top + node_h
    ax.text(50, node_bottom + 3.0,
            "Uncertainty about drug identity, evidence relevance and safety overlap carries forward into the final estimate.",
            ha="center", va="center", fontsize=8.0, color=C_MUTED, style="italic")

    merge_y = node_bottom + 6.6
    down_arrow(ax, 50, node_bottom + 4.6, merge_y - 0.6, color=C_NEUTRAL_BORDER)
    ax.text(50, merge_y + 1.6, "→ one evidence record: hydroxychloroquine × COVID-19",
            ha="center", va="center", fontsize=9.6, fontweight="bold", color=C_COMPOSITE)

    # =================================================================
    # PANEL C — final, conflicted evidence profile
    # =================================================================
    cy0 = merge_y + 8.5
    ax.text(2, cy0, "C", fontsize=15, fontweight="bold", color=C_INK, va="top")
    ax.text(6.5, cy0, "The evidence profile for this pair — and where it conflicts",
            fontsize=12.5, fontweight="bold", color=C_INK, va="top")

    profile_rows = [
        dict(label="Drug & disease mapping", value="exact canonical mapping", dot=C_LIT),
        dict(label="Literature coverage", value=f"{d['records_retrieved']:,} exact-pair records classified", dot=C_LIT),
        dict(label="Literature signal", value=f"mixed — {d['therapeutic']} therapeutic, {d['adverse']} adverse, {d['irrelevant']} irrelevant", dot=C_MIXED),
        dict(label="Safety overlap with disease symptoms", value=f"high, conflicting (γ = {d['gamma']:.2f})", dot=C_SAFETY),
        dict(label="Network structural support", value=f"held-out graph probability {d['graph_probability']:.3f}", dot=C_LIT),
        dict(label="Certainty of the final estimate", value=f"95% CrI {d['ci_low']:.3f}–{d['ci_high']:.3f}", dot=C_LIT),
    ]

    row_top = cy0 + 5.5
    row_h = 4.6
    label_x = 6.5
    value_x = 46
    dot_x = 4.2

    for i, r in enumerate(profile_rows):
        y = row_top + i * row_h
        ax.add_patch(Circle((dot_x, y), 0.75, facecolor=r["dot"], edgecolor="none", zorder=3))
        ax.text(label_x, y, r["label"], ha="left", va="center", fontsize=9.0,
                fontweight="bold", color=C_INK)
        ax.text(value_x, y, r["value"], ha="left", va="center", fontsize=9.0, color=C_MUTED)

    rows_bottom = row_top + len(profile_rows) * row_h
    divider_y = rows_bottom + 1.6
    ax.plot([2, 98], [divider_y, divider_y], color="#CCCCCC", lw=1.0)

    score_y = divider_y + 4.0
    ax.text(label_x, score_y, "Evidence-readiness score", ha="left", va="center",
            fontsize=10.8, fontweight="bold", color=C_INK)
    ax.text(value_x, score_y, f"{d['evidence_readiness_score']:.1f} / 100", ha="left", va="center",
            fontsize=13.5, fontweight="bold", color=C_COMPOSITE)

    ax.text(50, score_y + 4.6,
            f"This is not an efficacy probability — estimated probability of therapeutic benefit: "
            f"{d['posterior_mean']:.4f} (95% CrI {d['ci_low']:.3f}–{d['ci_high']:.3f}), reflecting the safety conflict above.",
            ha="center", va="center", fontsize=8.6, color=C_SAFETY, style="italic")

    fig_bottom = score_y + 9.0

    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(out_dir / f"Figure_case_study_hcq_covid_data_quality.{ext}",
                    dpi=300 if ext == "png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved figure to {out_dir} (content extends to y={fig_bottom:.1f} of ylim=128)")


if __name__ == "__main__":
    data = load_data()
    build_figure(data, OUT_DIR)
