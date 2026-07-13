#!/usr/bin/env python3
"""Move mapped aptamers to nearby upward-facing lattice sites.

Beginner picture:
The first lattice mapping picked the closest DNA-origami site.
This script asks a better question:

"Can we move each aptamer a tiny distance to a nearby site that points upward?"

That is like moving a hook from one pegboard hole to the next hole so the hook
faces the right way, while keeping the same overall broad-grid shape.
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

from score_ev_capture_geometry import load_linker_models
from score_lattice_orientation import (
    GOOD_UP_SCORE,
    orientation_metrics,
    oriented_linker_models,
    sparse_score_average,
)

ROOT = Path(__file__).resolve().parent
IN_SITES = ROOT / "origami_lattice_sites.csv"
IN_MAPPED = ROOT / "origami_lattice_mapped_layouts.csv"
OUT_LAYOUTS = ROOT / "orientation_optimized_mapped_layouts.csv"
OUT_SUMMARY = ROOT / "orientation_optimized_summary.csv"
OUT_JSON = ROOT / "orientation_optimized_summary.json"
OUT_PLOT = ROOT / "orientation_optimized_layouts.png"

MAX_SHIFT_FROM_ORIGINAL_NM = 4.0
SHIFT_PENALTY = 0.08

Anchor = dict[str, float | str]


def read_sites() -> list[dict[str, float | int]]:
    with open(IN_SITES, newline="", encoding="ascii") as f:
        return [
            {
                "site_id": int(row["site_id"]),
                "helix_id": int(row["helix_id"]),
                "base_index": int(row["base_index"]),
                "x_nm": float(row["x_nm"]),
                "y_nm": float(row["y_nm"]),
            }
            for row in csv.DictReader(f)
        ]


def read_original_mapped_rows() -> dict[str, list[dict[str, str]]]:
    rows_by_layout: dict[str, list[dict[str, str]]] = defaultdict(list)
    with open(IN_MAPPED, newline="", encoding="ascii") as f:
        for row in csv.DictReader(f):
            rows_by_layout[row["mapped_layout"]].append(row)
    return dict(rows_by_layout)


def candidate_score(site: dict[str, float | int], original_x: float, original_y: float) -> tuple[float, float]:
    shift = math.hypot(float(site["x_nm"]) - original_x, float(site["y_nm"]) - original_y)
    probe: Anchor = {
        "x_nm": float(site["x_nm"]),
        "y_nm": float(site["y_nm"]),
        "helix_id": int(site["helix_id"]),
        "base_index": int(site["base_index"]),
        "site_id": int(site["site_id"]),
        "linker_reach_nm": 15.0,
        "linker_construct": "polyT30",
    }
    up_score = float(orientation_metrics(probe)["up_score"])
    return up_score - SHIFT_PENALTY * shift, shift


def choose_site(
    original_x: float,
    original_y: float,
    sites: list[dict[str, float | int]],
    used_site_ids: set[int],
) -> tuple[dict[str, float | int], float]:
    candidates = []
    for site in sites:
        if int(site["site_id"]) in used_site_ids:
            continue
        shift = math.hypot(float(site["x_nm"]) - original_x, float(site["y_nm"]) - original_y)
        if shift <= MAX_SHIFT_FROM_ORIGINAL_NM:
            score, _ = candidate_score(site, original_x, original_y)
            candidates.append((score, -shift, site, shift))

    if not candidates:
        available = [site for site in sites if int(site["site_id"]) not in used_site_ids]
        site = min(
            available,
            key=lambda item: (float(item["x_nm"]) - original_x) ** 2
            + (float(item["y_nm"]) - original_y) ** 2,
        )
        shift = math.hypot(float(site["x_nm"]) - original_x, float(site["y_nm"]) - original_y)
        return site, shift

    _score, _negative_shift, site, shift = max(candidates, key=lambda item: (item[0], item[1]))
    return site, shift


def optimize_layout(
    layout_name: str,
    rows: list[dict[str, str]],
    sites: list[dict[str, float | int]],
) -> tuple[list[Anchor], list[dict[str, str]]]:
    used_site_ids: set[int] = set()
    anchors: list[Anchor] = []
    out_rows: list[dict[str, str]] = []

    for row in sorted(rows, key=lambda item: int(item["anchor_id"])):
        original_x = float(row["original_x_nm"])
        original_y = float(row["original_y_nm"])
        old_x = float(row["snapped_x_nm"])
        old_y = float(row["snapped_y_nm"])
        site, shift = choose_site(original_x, original_y, sites, used_site_ids)
        used_site_ids.add(int(site["site_id"]))
        anchor: Anchor = {
            "x_nm": float(site["x_nm"]),
            "y_nm": float(site["y_nm"]),
            "linker_reach_nm": 15.0,
            "linker_construct": row.get("linker_construct", "polyT30"),
            "anchor_id": int(row["anchor_id"]),
            "helix_id": int(site["helix_id"]),
            "base_index": int(site["base_index"]),
            "site_id": int(site["site_id"]),
        }
        metrics = orientation_metrics(anchor)
        anchors.append(anchor)
        out_rows.append(
            {
                "source_layout": layout_name,
                "optimized_layout": layout_name.replace("_lattice", "_orientation_optimized"),
                "anchor_id": row["anchor_id"],
                "original_x_nm": f"{original_x:.3f}",
                "original_y_nm": f"{original_y:.3f}",
                "old_lattice_x_nm": f"{old_x:.3f}",
                "old_lattice_y_nm": f"{old_y:.3f}",
                "optimized_x_nm": f"{float(site['x_nm']):.3f}",
                "optimized_y_nm": f"{float(site['y_nm']):.3f}",
                "shift_from_original_nm": f"{shift:.3f}",
                "shift_from_old_lattice_nm": f"{math.hypot(float(site['x_nm']) - old_x, float(site['y_nm']) - old_y):.3f}",
                "helix_id": str(int(site["helix_id"])),
                "base_index": str(int(site["base_index"])),
                "site_id": str(int(site["site_id"])),
                "orientation_class": str(metrics["orientation_class"]),
                "up_score": f"{float(metrics['up_score']):.4f}",
            }
        )
    return anchors, out_rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with open(path, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_layouts(rows: list[dict[str, str]]) -> None:
    by_layout: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_layout[row["optimized_layout"]].append(row)

    colors = {"up_facing": "#187a72", "side_facing": "#d69f22", "poor_facing": "#b23b3b"}
    fig, axes = plt.subplots(1, len(by_layout), figsize=(12, 4), constrained_layout=True)
    if len(by_layout) == 1:
        axes = [axes]
    for ax, (layout_name, layout_rows) in zip(axes, by_layout.items()):
        ax.add_patch(plt.Rectangle((-45, -30), 90, 60, fill=False, color="0.3"))
        for row in layout_rows:
            ax.scatter(
                float(row["optimized_x_nm"]),
                float(row["optimized_y_nm"]),
                s=52,
                color=colors[row["orientation_class"]],
                edgecolor="white",
                linewidth=0.5,
            )
        ax.set_title(layout_name)
        ax.set_xlim(-50, 50)
        ax.set_ylim(-35, 35)
        ax.set_aspect("equal")
        ax.grid(alpha=0.2)
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)


def main() -> None:
    sites = read_sites()
    source = read_original_mapped_rows()
    base_models = load_linker_models()
    all_layout_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []

    for layout_name, rows in source.items():
        anchors, layout_rows = optimize_layout(layout_name, rows, sites)
        all_layout_rows.extend(layout_rows)
        oriented, models, anchor_scores = oriented_linker_models(anchors, base_models, layout_name)
        score = sparse_score_average(oriented, models)

        up_count = sum(row["orientation_class"] == "up_facing" for row in layout_rows)
        side_count = sum(row["orientation_class"] == "side_facing" for row in layout_rows)
        poor_count = sum(row["orientation_class"] == "poor_facing" for row in layout_rows)
        shifts = np.asarray([float(row["shift_from_original_nm"]) for row in layout_rows])
        up_scores = np.asarray([float(row["up_score"]) for row in layout_rows])
        summary_rows.append(
            {
                "source_layout": layout_name,
                "optimized_layout": layout_name.replace("_lattice", "_orientation_optimized"),
                "aptamer_count": str(len(anchors)),
                "up_facing_aptamers": str(up_count),
                "side_facing_aptamers": str(side_count),
                "poor_facing_aptamers": str(poor_count),
                "mean_up_score": f"{float(np.mean(up_scores)):.4f}",
                "mean_shift_from_original_nm": f"{float(np.mean(shifts)):.3f}",
                "max_shift_from_original_nm": f"{float(np.max(shifts)):.3f}",
                "oriented_sparse_score": f"{score['mean_clinical_sparse_score']:.4f}",
                "oriented_mean_contacts": f"{score['mean_contacts']:.4f}",
                "oriented_p_at_least_1_contact": f"{score['mean_p_at_least_1_contact']:.4f}",
                "oriented_p_at_least_2_contacts": f"{score['mean_p_at_least_2_contacts']:.4f}",
                "oriented_p_at_least_3_contacts": f"{score['mean_p_at_least_3_contacts']:.4f}",
            }
        )

    summary_rows.sort(
        key=lambda row: (
            float(row["oriented_sparse_score"]),
            float(row["mean_up_score"]),
            -float(row["mean_shift_from_original_nm"]),
        ),
        reverse=True,
    )
    write_csv(OUT_LAYOUTS, all_layout_rows)
    write_csv(OUT_SUMMARY, summary_rows)
    plot_layouts(all_layout_rows)

    report = {
        "model": "orientation-aware lattice remapper",
        "plain_language_summary": [
            "The first lattice map chose the nearest attachment site.",
            "This optimizer chooses a nearby site that faces upward when possible.",
            "The goal is to keep the broad layout shape while improving aptamer display direction.",
        ],
        "settings": {
            "maximum_shift_from_original_nm": MAX_SHIFT_FROM_ORIGINAL_NM,
            "shift_penalty": SHIFT_PENALTY,
            "good_up_score_threshold": GOOD_UP_SCORE,
        },
        "best_optimized_layout": summary_rows[0],
        "ranked_layouts": summary_rows,
        "limitations": [
            "This is still a simplified register/orientation model.",
            "It does not replace caDNAno routing or full oxDNA validation.",
            "The optimized sites should be checked against a real tile design before ordering strands.",
        ],
        "outputs": {
            "optimized_layout_coordinates": OUT_LAYOUTS.name,
            "optimized_summary": OUT_SUMMARY.name,
            "optimized_layout_plot": OUT_PLOT.name,
        },
    }
    OUT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="ascii")

    best = summary_rows[0]
    print(f"Wrote {OUT_LAYOUTS.name}")
    print(f"Wrote {OUT_SUMMARY.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_PLOT.name}")
    print(
        f"Best orientation-optimized layout: {best['optimized_layout']} "
        f"oriented_sparse={best['oriented_sparse_score']} "
        f"up={best['up_facing_aptamers']}/{best['aptamer_count']} "
        f"mean_shift={best['mean_shift_from_original_nm']} nm"
    )


if __name__ == "__main__":
    main()
