#!/usr/bin/env python3
"""Score whether capture layouts are practical DNA-origami patterns.

Beginner picture:
The other scripts mostly ask, "Can this pattern catch the EV?"
This script asks a second question: "Could this pattern be a reasonable
DNA-origami tile?"

Definitions used here:
* anchor: one spot on the DNA origami where an aptamer can be attached.
* nearest-neighbor spacing: distance from one anchor to the closest anchor.
* crowding: anchors are so close that flexible linkers or aptamers may bump.
* edge margin: distance from an anchor to the edge of the DNA-origami tile.

The constants below are screening assumptions, not hard experimental laws.
They make the design search more honest by penalizing patterns that are too
packed, too edge-heavy, or too overloaded with aptamers.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from advanced_capture_design_search import candidate_layouts

ROOT = Path(__file__).resolve().parent
IN_ADVANCED = ROOT / "advanced_capture_design_summary.csv"
IN_WASH = ROOT / "capture_wash_protocol_best.csv"
OUT_CSV = ROOT / "origami_buildability_scores.csv"
OUT_JSON = ROOT / "origami_buildability_summary.json"
OUT_PLOT = ROOT / "origami_buildability_scores.png"
OUT_LAYOUTS = ROOT / "origami_buildability_layouts.png"

TILE_WIDTH_NM = 90.0
TILE_HEIGHT_NM = 60.0

# Screening assumptions. These values are deliberately simple and explainable.
HARD_MIN_SPACING_NM = 6.0
COMFORTABLE_SPACING_NM = 9.0
COMFORTABLE_EDGE_MARGIN_NM = 5.0
LOCAL_CROWDING_RADIUS_NM = 12.0
MAX_COMFORTABLE_LOCAL_NEIGHBORS = 4
MAX_COMFORTABLE_APTAMERS = 24


def read_best_capture_scores() -> dict[tuple[str, str], dict[str, str]]:
    """Prefer final wash-optimized capture results when available."""
    source = IN_WASH if IN_WASH.exists() else IN_ADVANCED
    with open(source, newline="", encoding="ascii") as f:
        rows = list(csv.DictReader(f))
    return {(row["layout"], row["formulation"]): row for row in rows}


def pairwise_distances(xy: np.ndarray) -> np.ndarray:
    delta = xy[:, None, :] - xy[None, :, :]
    distances = np.sqrt(np.sum(delta * delta, axis=2))
    np.fill_diagonal(distances, np.inf)
    return distances


def clamp01(value: float) -> float:
    return max(0.0, min(1.0, value))


def buildability_metrics(layout: list[dict[str, float | str]]) -> dict[str, float]:
    xy = np.asarray([[float(anchor["x_nm"]), float(anchor["y_nm"])] for anchor in layout])
    distances = pairwise_distances(xy)
    nearest = np.min(distances, axis=1)

    edge_margins = np.minimum(
        TILE_WIDTH_NM / 2.0 - np.abs(xy[:, 0]),
        TILE_HEIGHT_NM / 2.0 - np.abs(xy[:, 1]),
    )
    local_neighbor_counts = np.sum(distances < LOCAL_CROWDING_RADIUS_NM, axis=1)

    hard_spacing_violations = int(np.sum(nearest < HARD_MIN_SPACING_NM))
    crowding_violations = int(np.sum(nearest < COMFORTABLE_SPACING_NM))
    edge_violations = int(np.sum(edge_margins < COMFORTABLE_EDGE_MARGIN_NM))
    max_local_neighbors = int(np.max(local_neighbor_counts))

    min_spacing = float(np.min(nearest))
    median_spacing = float(np.median(nearest))
    min_edge_margin = float(np.min(edge_margins))
    aptamer_count = len(layout)

    spacing_score = clamp01((median_spacing - HARD_MIN_SPACING_NM) / (COMFORTABLE_SPACING_NM - HARD_MIN_SPACING_NM))
    if hard_spacing_violations:
        spacing_score *= 0.65 ** hard_spacing_violations

    edge_score = clamp01(min_edge_margin / COMFORTABLE_EDGE_MARGIN_NM)
    local_score = clamp01(
        1.0 - max(0, max_local_neighbors - MAX_COMFORTABLE_LOCAL_NEIGHBORS)
        / MAX_COMFORTABLE_LOCAL_NEIGHBORS
    )
    count_score = clamp01(1.0 - max(0, aptamer_count - MAX_COMFORTABLE_APTAMERS) / 12.0)

    # Weighted average: spacing matters most because crowding is the most direct
    # way to make aptamers interfere with each other on a small tile.
    origami_score = (
        0.42 * spacing_score
        + 0.22 * local_score
        + 0.18 * edge_score
        + 0.18 * count_score
    )

    return {
        "aptamer_count": float(aptamer_count),
        "min_spacing_nm": min_spacing,
        "median_nearest_spacing_nm": median_spacing,
        "hard_spacing_violations": float(hard_spacing_violations),
        "crowding_violations": float(crowding_violations),
        "min_edge_margin_nm": min_edge_margin,
        "edge_violations": float(edge_violations),
        "max_local_neighbors": float(max_local_neighbors),
        "spacing_score": spacing_score,
        "local_crowding_score": local_score,
        "edge_score": edge_score,
        "aptamer_count_score": count_score,
        "origami_buildability_score": origami_score,
    }


def fmt(value: float) -> str:
    if math.isclose(value, round(value)):
        return str(int(round(value)))
    return f"{value:.4f}"


def main() -> None:
    layouts = candidate_layouts()
    capture = read_best_capture_scores()
    rows: list[dict[str, str]] = []

    for layout_name, layout in layouts.items():
        metrics = buildability_metrics(layout)
        capture_rows = [
            row for (name, _formulation), row in capture.items()
            if name == layout_name
        ]
        if capture_rows:
            best_capture = max(
                capture_rows,
                key=lambda row: float(row.get("p10_useful_capture", row.get("median_useful_capture", 0.0))),
            )
            formulation = best_capture["formulation"]
            capture_score = float(best_capture.get("p10_useful_capture", best_capture.get("median_useful_capture", 0.0)))
            reliable_90 = best_capture.get("reliable_90", best_capture.get("reliable_90_useful", "no"))
        else:
            formulation = "not_tested"
            capture_score = 0.0
            reliable_90 = "no"

        # This combined score keeps capture as the main goal but rewards layouts
        # that are easier to explain, fabricate, and troubleshoot experimentally.
        practical_score = 0.72 * capture_score + 0.28 * metrics["origami_buildability_score"]

        row = {
            "layout": layout_name,
            "best_formulation_or_protocol": formulation,
            "capture_score_used": f"{capture_score:.4f}",
            "reliable_90_capture": reliable_90,
            **{key: fmt(value) for key, value in metrics.items()},
            "origami_adjusted_practical_score": f"{practical_score:.4f}",
        }
        rows.append(row)

    rows.sort(key=lambda row: float(row["origami_adjusted_practical_score"]), reverse=True)
    with open(OUT_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(10, 5.8), constrained_layout=True)
    top = rows[:10][::-1]
    labels = [row["layout"] for row in top]
    capture_values = [float(row["capture_score_used"]) for row in top]
    build_values = [float(row["origami_buildability_score"]) for row in top]
    y = np.arange(len(top))
    ax.barh(y - 0.18, capture_values, height=0.34, label="capture after protocol", color="#187a72")
    ax.barh(y + 0.18, build_values, height=0.34, label="origami buildability", color="#5b6ee1")
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 1)
    ax.set_xlabel("score from 0 to 1")
    ax.set_title("Capture performance compared with DNA-origami buildability")
    ax.grid(axis="x", alpha=0.25)
    ax.legend()
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(11, 7), constrained_layout=True)
    for ax, row in zip(axes.flat, rows[:6]):
        layout = layouts[row["layout"]]
        xy = np.asarray([[float(anchor["x_nm"]), float(anchor["y_nm"])] for anchor in layout])
        ax.add_patch(plt.Rectangle((-45, -30), 90, 60, fill=False, color="0.3"))
        ax.scatter(xy[:, 0], xy[:, 1], s=38, color="#187a72")
        ax.set_title(
            f"{row['layout']}\n"
            f"capture {row['capture_score_used']}, origami {row['origami_buildability_score']}"
        )
        ax.set_xlim(-50, 50)
        ax.set_ylim(-35, 35)
        ax.set_aspect("equal")
        ax.grid(alpha=0.2)
    fig.savefig(OUT_LAYOUTS, dpi=220)
    plt.close(fig)

    best = rows[0]
    summary = {
        "model": "DNA-origami buildability screen for EV capture layouts",
        "plain_language_summary": [
            "A high capture score means the simulated EV is likely to be caught.",
            "A high origami buildability score means the aptamer anchors are spaced in a less crowded, more tile-friendly way.",
            "The practical score combines both, with capture still weighted more strongly than buildability.",
        ],
        "screening_assumptions": {
            "tile_width_nm": TILE_WIDTH_NM,
            "tile_height_nm": TILE_HEIGHT_NM,
            "hard_min_spacing_nm": HARD_MIN_SPACING_NM,
            "comfortable_spacing_nm": COMFORTABLE_SPACING_NM,
            "comfortable_edge_margin_nm": COMFORTABLE_EDGE_MARGIN_NM,
            "local_crowding_radius_nm": LOCAL_CROWDING_RADIUS_NM,
            "max_comfortable_local_neighbors": MAX_COMFORTABLE_LOCAL_NEIGHBORS,
            "max_comfortable_aptamers": MAX_COMFORTABLE_APTAMERS,
        },
        "best_practical_layout": best,
        "top_three": rows[:3],
        "limitations": [
            "This is a coarse buildability screen, not a caDNAno or oxDNA full-tile folding simulation.",
            "It does not assign exact staple sequences or check scaffold routing.",
            "The spacing thresholds should be replaced by lab-specific attachment chemistry constraints when available.",
        ],
        "outputs": {
            "ranked_scores": OUT_CSV.name,
            "score_plot": OUT_PLOT.name,
            "top_layout_plot": OUT_LAYOUTS.name,
        },
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
    print(f"Wrote {OUT_CSV.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_PLOT.name}")
    print(f"Wrote {OUT_LAYOUTS.name}")
    print(
        f"Best practical layout: {best['layout']} "
        f"capture={best['capture_score_used']} "
        f"origami={best['origami_buildability_score']} "
        f"combined={best['origami_adjusted_practical_score']}"
    )


if __name__ == "__main__":
    main()
