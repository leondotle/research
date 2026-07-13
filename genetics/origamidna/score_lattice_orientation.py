#!/usr/bin/env python3
"""Score aptamer orientation on mapped DNA-origami lattice sites.

Beginner picture:
The lattice mapper tells us where each aptamer can attach.
This script asks which way each aptamer probably points.

Why this matters:
An aptamer is useful only if it can reach upward toward the EV. If it points
sideways or partly into the tile, it may still exist chemically, but it is less
useful for capture.

This is a simplified orientation model, not a full molecular simulation. It
adds an honest penalty for attachment sites that are not predicted to face the
EV well.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import score_ev_capture_geometry as scg
from score_ev_capture_clinical_73nm import enrich_sparse_metrics
from score_ev_capture_geometry import load_linker_models, score_layout

ROOT = Path(__file__).resolve().parent
IN_MAPPED = ROOT / "origami_lattice_mapped_layouts.csv"
OUT_ANCHORS = ROOT / "lattice_orientation_anchor_scores.csv"
OUT_SUMMARY = ROOT / "lattice_orientation_summary.csv"
OUT_JSON = ROOT / "lattice_orientation_summary.json"
OUT_PLOT = ROOT / "lattice_orientation_scores.png"
OUT_LAYOUT_PLOT = ROOT / "lattice_orientation_layouts.png"

EV_DIAMETER_NM = 73.0
RECEPTOR_COUNTS = (2, 5, 10)
RNG_SEED = 20260712

# Simplified orientation assumptions.
# The phase model says neighboring helix rows and neighboring attachment
# registers do not all point in exactly the same direction.
REGISTER_COUNT = 3
HELIX_STAGGER_RADIANS = math.pi / 3.0
GOOD_UP_SCORE = 0.75
WEAK_UP_SCORE = 0.45

Anchor = dict[str, float | str]
LinkerModels = dict[str, tuple[np.ndarray, np.ndarray]]


def read_mapped_layouts() -> dict[str, list[Anchor]]:
    layouts: dict[str, list[Anchor]] = defaultdict(list)
    with open(IN_MAPPED, newline="", encoding="ascii") as f:
        for row in csv.DictReader(f):
            layouts[row["mapped_layout"]].append(
                {
                    "x_nm": float(row["snapped_x_nm"]),
                    "y_nm": float(row["snapped_y_nm"]),
                    "linker_reach_nm": 15.0,
                    "linker_construct": row.get("linker_construct", "polyT30"),
                    "anchor_id": int(row["anchor_id"]),
                    "helix_id": int(row["helix_id"]),
                    "base_index": int(row["base_index"]),
                    "site_id": int(row["site_id"]),
                }
            )
    return dict(layouts)


def orientation_metrics(anchor: Anchor) -> dict[str, float | str]:
    """Estimate whether the site faces upward toward the EV.

    up_score ranges from 0 to 1.
    1 means strongly upward-facing.
    0 means strongly downward/blocked-facing.
    """
    helix_id = int(anchor["helix_id"])
    base_index = int(anchor["base_index"])
    phase = (
        2.0 * math.pi * ((base_index % REGISTER_COUNT) / REGISTER_COUNT)
        + HELIX_STAGGER_RADIANS * (helix_id % 2)
    )
    up_score = 0.5 + 0.5 * math.cos(phase)
    if up_score >= GOOD_UP_SCORE:
        class_name = "up_facing"
    elif up_score >= WEAK_UP_SCORE:
        class_name = "side_facing"
    else:
        class_name = "poor_facing"

    # Side-facing aptamers are not useless, but their practical reach is lower.
    reach_multiplier = 0.55 + 0.45 * up_score
    binding_multiplier = 0.35 + 0.65 * up_score
    return {
        "phase_radians": phase,
        "up_score": up_score,
        "orientation_class": class_name,
        "reach_multiplier": reach_multiplier,
        "binding_multiplier": binding_multiplier,
    }


def oriented_linker_models(
    anchors: list[Anchor],
    base_models: LinkerModels,
    layout_name: str,
) -> tuple[list[Anchor], LinkerModels, list[dict[str, str]]]:
    models = dict(base_models)
    oriented: list[Anchor] = []
    rows: list[dict[str, str]] = []

    for anchor in anchors:
        metrics = orientation_metrics(anchor)
        base_construct = str(anchor["linker_construct"])
        xs, ys = base_models[base_construct]
        construct_name = f"{layout_name}_anchor_{int(anchor['anchor_id'])}"
        models[construct_name] = (
            xs * float(metrics["reach_multiplier"]),
            ys * float(metrics["binding_multiplier"]),
        )
        oriented_anchor = {
            **anchor,
            "linker_construct": construct_name,
            "linker_reach_nm": float(anchor["linker_reach_nm"]) * float(metrics["reach_multiplier"]),
        }
        oriented.append(oriented_anchor)
        rows.append(
            {
                "layout": layout_name,
                "anchor_id": str(int(anchor["anchor_id"])),
                "x_nm": f"{float(anchor['x_nm']):.3f}",
                "y_nm": f"{float(anchor['y_nm']):.3f}",
                "helix_id": str(int(anchor["helix_id"])),
                "base_index": str(int(anchor["base_index"])),
                "site_id": str(int(anchor["site_id"])),
                "orientation_class": str(metrics["orientation_class"]),
                "up_score": f"{float(metrics['up_score']):.4f}",
                "reach_multiplier": f"{float(metrics['reach_multiplier']):.4f}",
                "binding_multiplier": f"{float(metrics['binding_multiplier']):.4f}",
            }
        )
    return oriented, models, rows


def sparse_score_average(layout: list[Anchor], linker_models: LinkerModels) -> dict[str, float]:
    scg.LATERAL_STEP_NM = 8.0
    scg.N_RECEPTOR_REALIZATIONS = 3
    scg.N_BINDING_TRIALS = 6

    rng = np.random.default_rng(RNG_SEED)
    scores = []
    mean_contacts = []
    p1_values = []
    p2_values = []
    p3_values = []
    for receptor_count in RECEPTOR_COUNTS:
        metrics = score_layout(
            layout,
            linker_models,
            EV_DIAMETER_NM,
            0.0,
            rng,
            fixed_receptor_count=receptor_count,
        )
        enriched = enrich_sparse_metrics(metrics)
        scores.append(enriched["clinical_sparse_score"])
        mean_contacts.append(enriched["mean_contacts"])
        p1_values.append(enriched["p_at_least_1_contact"])
        p2_values.append(enriched["p_at_least_2_contacts"])
        p3_values.append(enriched["p_at_least_3_contacts"])

    return {
        "mean_clinical_sparse_score": float(np.mean(scores)),
        "mean_contacts": float(np.mean(mean_contacts)),
        "mean_p_at_least_1_contact": float(np.mean(p1_values)),
        "mean_p_at_least_2_contacts": float(np.mean(p2_values)),
        "mean_p_at_least_3_contacts": float(np.mean(p3_values)),
    }


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with open(path, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_summary(rows: list[dict[str, str]]) -> None:
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    labels = [row["layout"] for row in rows]
    x = np.arange(len(rows))
    ax.bar(x - 0.18, [float(row["unoriented_sparse_score"]) for row in rows], width=0.36, label="before orientation penalty", color="#187a72")
    ax.bar(x + 0.18, [float(row["oriented_sparse_score"]) for row in rows], width=0.36, label="after orientation penalty", color="#c45a2c")
    ax.set_xticks(x, labels, rotation=20, ha="right")
    ax.set_ylabel("mean clinical sparse score")
    ax.set_ylim(0, max(0.25, max(float(row["unoriented_sparse_score"]) for row in rows) * 1.15))
    ax.set_title("Effect of aptamer orientation on mapped lattice layouts")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)


def plot_layouts(layouts: dict[str, list[Anchor]], anchor_rows: list[dict[str, str]]) -> None:
    class_colors = {
        "up_facing": "#187a72",
        "side_facing": "#d69f22",
        "poor_facing": "#b23b3b",
    }
    rows_by_layout: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in anchor_rows:
        rows_by_layout[row["layout"]].append(row)

    fig, axes = plt.subplots(1, len(layouts), figsize=(12, 4), constrained_layout=True)
    if len(layouts) == 1:
        axes = [axes]
    for ax, (layout_name, anchors) in zip(axes, layouts.items()):
        ax.add_patch(plt.Rectangle((-45, -30), 90, 60, fill=False, color="0.3"))
        for row in rows_by_layout[layout_name]:
            ax.scatter(
                float(row["x_nm"]),
                float(row["y_nm"]),
                s=52,
                color=class_colors[row["orientation_class"]],
                edgecolor="white",
                linewidth=0.5,
            )
        ax.set_title(layout_name)
        ax.set_xlim(-50, 50)
        ax.set_ylim(-35, 35)
        ax.set_aspect("equal")
        ax.grid(alpha=0.2)
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, label=label, markersize=8)
        for label, color in class_colors.items()
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3)
    fig.savefig(OUT_LAYOUT_PLOT, dpi=220)
    plt.close(fig)


def main() -> None:
    base_models = load_linker_models()
    layouts = read_mapped_layouts()
    anchor_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []
    oriented_layouts: dict[str, list[Anchor]] = {}

    for layout_name, anchors in layouts.items():
        unoriented = sparse_score_average(anchors, base_models)
        oriented, oriented_models, rows = oriented_linker_models(anchors, base_models, layout_name)
        oriented_layouts[layout_name] = oriented
        anchor_rows.extend(rows)
        oriented_score = sparse_score_average(oriented, oriented_models)

        up_scores = np.asarray([float(row["up_score"]) for row in rows])
        class_counts = {
            class_name: sum(row["orientation_class"] == class_name for row in rows)
            for class_name in ("up_facing", "side_facing", "poor_facing")
        }
        summary_rows.append(
            {
                "layout": layout_name,
                "aptamer_count": str(len(anchors)),
                "up_facing_aptamers": str(class_counts["up_facing"]),
                "side_facing_aptamers": str(class_counts["side_facing"]),
                "poor_facing_aptamers": str(class_counts["poor_facing"]),
                "mean_up_score": f"{float(np.mean(up_scores)):.4f}",
                "min_up_score": f"{float(np.min(up_scores)):.4f}",
                "unoriented_sparse_score": f"{unoriented['mean_clinical_sparse_score']:.4f}",
                "oriented_sparse_score": f"{oriented_score['mean_clinical_sparse_score']:.4f}",
                "orientation_score_change": f"{oriented_score['mean_clinical_sparse_score'] - unoriented['mean_clinical_sparse_score']:.4f}",
                "oriented_mean_contacts": f"{oriented_score['mean_contacts']:.4f}",
                "oriented_p_at_least_1_contact": f"{oriented_score['mean_p_at_least_1_contact']:.4f}",
                "oriented_p_at_least_2_contacts": f"{oriented_score['mean_p_at_least_2_contacts']:.4f}",
                "oriented_p_at_least_3_contacts": f"{oriented_score['mean_p_at_least_3_contacts']:.4f}",
            }
        )

    summary_rows.sort(
        key=lambda row: (
            float(row["oriented_sparse_score"]),
            float(row["mean_up_score"]),
        ),
        reverse=True,
    )
    write_csv(OUT_ANCHORS, anchor_rows)
    write_csv(OUT_SUMMARY, summary_rows)
    plot_summary(summary_rows)
    plot_layouts({row["layout"]: layouts[row["layout"]] for row in summary_rows}, anchor_rows)

    summary = {
        "model": "simplified aptamer orientation screen for lattice-mapped DNA origami",
        "plain_language_summary": [
            "The lattice mapper decides where aptamers attach.",
            "This orientation screen estimates whether those attachment sites point upward toward the EV.",
            "Poor-facing aptamers keep existing as anchors, but the model reduces their effective reach and binding probability.",
        ],
        "orientation_assumptions": {
            "register_count": REGISTER_COUNT,
            "helix_stagger_radians": HELIX_STAGGER_RADIANS,
            "good_up_score_threshold": GOOD_UP_SCORE,
            "weak_up_score_threshold": WEAK_UP_SCORE,
            "reach_multiplier_range": "0.55 to 1.00",
            "binding_multiplier_range": "0.35 to 1.00",
        },
        "best_oriented_layout": summary_rows[0],
        "ranked_layouts": summary_rows,
        "limitations": [
            "This is not a full nucleotide-level orientation simulation.",
            "It assumes a simplified repeating orientation pattern across helix rows and attachment registers.",
            "A real caDNAno/oxDNA model is needed to confirm exact aptamer display direction.",
        ],
        "outputs": {
            "anchor_scores": OUT_ANCHORS.name,
            "layout_summary": OUT_SUMMARY.name,
            "summary_plot": OUT_PLOT.name,
            "layout_plot": OUT_LAYOUT_PLOT.name,
        },
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")

    best = summary_rows[0]
    print(f"Wrote {OUT_ANCHORS.name}")
    print(f"Wrote {OUT_SUMMARY.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_PLOT.name}")
    print(f"Wrote {OUT_LAYOUT_PLOT.name}")
    print(
        f"Best oriented layout: {best['layout']} "
        f"oriented_sparse={best['oriented_sparse_score']} "
        f"up={best['up_facing_aptamers']}/{best['aptamer_count']}"
    )


if __name__ == "__main__":
    main()
