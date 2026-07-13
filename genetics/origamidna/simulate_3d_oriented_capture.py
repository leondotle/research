#!/usr/bin/env python3
"""Simplified 3D capture simulation for orientation-optimized origami layouts.

Beginner picture:
The older model mostly asked, "Is the receptor close enough in x/y space?"
This script asks a more 3D question:

"Is the receptor inside the aptamer's 3D reach cone?"

Definitions:
* 3D position: x, y, and z location.
* direction vector: which way an aptamer points in 3D.
* reach cone: the 3D region an aptamer can plausibly explore.
* receptor: a CD133 point on the lower half of a spherical EV.

This is still a coarse model. It is not a full oxDNA molecular simulation, but
it is a real 3D geometry screen: receptors must be close enough and in the
right direction from the aptamer.
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

from ev_population_generator import receptor_points
from score_ev_capture_geometry import (
    EV_SURFACE_CLEARANCE_NM,
    RNG_SEED,
    load_linker_models,
    max_bipartite_matches,
    reach_probability,
)
from score_lattice_orientation import REGISTER_COUNT, HELIX_STAGGER_RADIANS

ROOT = Path(__file__).resolve().parent
IN_ORIENTED = ROOT / "orientation_optimized_mapped_layouts.csv"
OUT_CSV = ROOT / "capture_3d_oriented_scores.csv"
OUT_JSON = ROOT / "capture_3d_oriented_summary.json"
OUT_PLOT = ROOT / "capture_3d_oriented_scores.png"
OUT_LAYOUT_PLOT = ROOT / "capture_3d_oriented_layouts.png"

EV_DIAMETER_NM = 73.0
RECEPTOR_COUNTS = (2, 5, 10)
RECEPTOR_PATTERNS = ("random", "single_cluster", "two_cluster", "bottom_cap")
PATTERN_PROBABILITIES = (0.40, 0.25, 0.20, 0.15)

LATERAL_SCAN_NM = 24.0
LATERAL_STEP_NM = 8.0
N_RECEPTOR_REALIZATIONS = 8
N_BINDING_TRIALS = 12

APTAMER_BASE_HEIGHT_NM = 1.1
REACH_CONE_HALF_ANGLE_DEG = 70.0
SOFT_OUTER_ANGLE_DEG = 90.0
MIN_DIRECTIONAL_PROBABILITY = 0.08
SEED = RNG_SEED + 909

Anchor = dict[str, float | str]
LinkerModels = dict[str, tuple[np.ndarray, np.ndarray]]


def read_orientation_optimized_layouts() -> dict[str, list[Anchor]]:
    layouts: dict[str, list[Anchor]] = defaultdict(list)
    with open(IN_ORIENTED, newline="", encoding="ascii") as f:
        for row in csv.DictReader(f):
            layouts[row["optimized_layout"]].append(
                {
                    "x_nm": float(row["optimized_x_nm"]),
                    "y_nm": float(row["optimized_y_nm"]),
                    "z_nm": APTAMER_BASE_HEIGHT_NM,
                    "linker_reach_nm": 15.0,
                    "linker_construct": "polyT30",
                    "anchor_id": int(row["anchor_id"]),
                    "helix_id": int(row["helix_id"]),
                    "base_index": int(row["base_index"]),
                    "orientation_class": row["orientation_class"],
                    "up_score": float(row["up_score"]),
                }
            )
    return dict(layouts)


def direction_vector(anchor: Anchor) -> np.ndarray:
    """Return a simple 3D direction vector for an aptamer attachment site."""
    helix_id = int(anchor["helix_id"])
    base_index = int(anchor["base_index"])
    phase = (
        2.0 * math.pi * ((base_index % REGISTER_COUNT) / REGISTER_COUNT)
        + HELIX_STAGGER_RADIANS * (helix_id % 2)
    )
    # Treat each helix as running along x; the surface normal rotates in y/z.
    vector = np.array([0.0, math.sin(phase), math.cos(phase)], dtype=float)
    if vector[2] < 0:
        # Linkers are flexible, so a poor register can bend upward somewhat,
        # but it starts at a disadvantage.
        vector[2] *= 0.35
    vector /= np.linalg.norm(vector)
    return vector


def angular_factor(anchor_to_receptor: np.ndarray, direction: np.ndarray) -> float:
    distance = float(np.linalg.norm(anchor_to_receptor))
    if distance <= 1e-9:
        return 1.0
    cos_angle = float(np.dot(anchor_to_receptor / distance, direction))
    cos_inner = math.cos(math.radians(REACH_CONE_HALF_ANGLE_DEG))
    cos_outer = math.cos(math.radians(SOFT_OUTER_ANGLE_DEG))
    if cos_angle >= cos_inner:
        return 1.0
    if cos_angle <= cos_outer:
        return MIN_DIRECTIONAL_PROBABILITY
    return MIN_DIRECTIONAL_PROBABILITY + (1.0 - MIN_DIRECTIONAL_PROBABILITY) * (
        (cos_angle - cos_outer) / (cos_inner - cos_outer)
    )


def probability_matrix(
    anchors: list[Anchor],
    receptors: np.ndarray,
    linker_models: LinkerModels,
) -> np.ndarray:
    probabilities = np.zeros((len(anchors), len(receptors)), dtype=float)
    for i, anchor in enumerate(anchors):
        anchor_xyz = np.array(
            [float(anchor["x_nm"]), float(anchor["y_nm"]), float(anchor["z_nm"])],
            dtype=float,
        )
        direction = direction_vector(anchor)
        for j, receptor_xyz in enumerate(receptors):
            vector = receptor_xyz - anchor_xyz
            distance = float(np.linalg.norm(vector))
            reach = reach_probability(linker_models, str(anchor["linker_construct"]), distance)
            probabilities[i, j] = reach * angular_factor(vector, direction)
    return probabilities


def contact_metrics_3d(
    anchors: list[Anchor],
    receptor_body: np.ndarray,
    ev_radius_nm: float,
    offset_x_nm: float,
    offset_y_nm: float,
    linker_models: LinkerModels,
    rng: np.random.Generator,
) -> tuple[float, float, float, float, float]:
    center = np.array(
        [offset_x_nm, offset_y_nm, ev_radius_nm + EV_SURFACE_CLEARANCE_NM],
        dtype=float,
    )
    receptors = receptor_body + center
    probabilities = probability_matrix(anchors, receptors, linker_models)
    possible_contacts = float(max_bipartite_matches(probabilities > 0.05))
    sampled = np.empty(N_BINDING_TRIALS, dtype=float)
    for trial in range(N_BINDING_TRIALS):
        sampled[trial] = max_bipartite_matches(rng.random(probabilities.shape) < probabilities)
    return (
        float(np.mean(sampled)),
        possible_contacts,
        float(np.mean(sampled >= 1)),
        float(np.mean(sampled >= 2)),
        float(np.mean(sampled >= 3)),
    )


def score_layout_3d(
    anchors: list[Anchor],
    receptor_count: int,
    linker_models: LinkerModels,
    rng: np.random.Generator,
) -> dict[str, float]:
    ev_radius = EV_DIAMETER_NM / 2.0
    offsets = np.arange(-LATERAL_SCAN_NM, LATERAL_SCAN_NM + 0.001, LATERAL_STEP_NM)
    mean_contacts = []
    possible_contacts = []
    p1_values = []
    p2_values = []
    p3_values = []

    for ox in offsets:
        for oy in offsets:
            for _ in range(N_RECEPTOR_REALIZATIONS):
                pattern = str(rng.choice(RECEPTOR_PATTERNS, p=PATTERN_PROBABILITIES))
                receptor_body = receptor_points(pattern, receptor_count, ev_radius, rng)
                mean_c, possible_c, p1, p2, p3 = contact_metrics_3d(
                    anchors,
                    receptor_body,
                    ev_radius,
                    float(ox),
                    float(oy),
                    linker_models,
                    rng,
                )
                mean_contacts.append(mean_c)
                possible_contacts.append(possible_c)
                p1_values.append(p1)
                p2_values.append(p2)
                p3_values.append(p3)

    mean_contact_value = float(np.mean(mean_contacts))
    p1 = float(np.mean(p1_values))
    p2 = float(np.mean(p2_values))
    p3 = float(np.mean(p3_values))
    sparse_score = 0.40 * p1 + 0.30 * p2 + 0.20 * p3 + 0.10 * min(mean_contact_value / 3.0, 1.0)
    return {
        "mean_contacts": mean_contact_value,
        "max_possible_contacts": float(np.max(possible_contacts)),
        "p_at_least_1_contact": p1,
        "p_at_least_2_contacts": p2,
        "p_at_least_3_contacts": p3,
        "clinical_sparse_3d_score": sparse_score,
    }


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with open(path, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_scores(rows: list[dict[str, str]]) -> None:
    layouts = sorted({row["layout"] for row in rows})
    x = np.arange(len(RECEPTOR_COUNTS))
    width = 0.24
    fig, ax = plt.subplots(figsize=(9.5, 5.2), constrained_layout=True)
    for i, layout in enumerate(layouts):
        subset = {int(row["receptor_count"]): row for row in rows if row["layout"] == layout}
        ax.bar(
            x + (i - (len(layouts) - 1) / 2.0) * width,
            [float(subset[count]["clinical_sparse_3d_score"]) for count in RECEPTOR_COUNTS],
            width=width,
            label=layout.replace("_orientation_optimized", ""),
        )
    ax.set_xticks(x, [str(count) for count in RECEPTOR_COUNTS])
    ax.set_xlabel("CD133 receptors per 73 nm EV")
    ax.set_ylabel("3D clinical sparse score")
    ax.set_ylim(0, 1)
    ax.set_title("Simplified 3D oriented capture screen")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)


def plot_layouts(layouts: dict[str, list[Anchor]]) -> None:
    fig, axes = plt.subplots(1, len(layouts), figsize=(12, 4), constrained_layout=True)
    if len(layouts) == 1:
        axes = [axes]
    for ax, (layout_name, anchors) in zip(axes, layouts.items()):
        xy = np.asarray([[float(a["x_nm"]), float(a["y_nm"])] for a in anchors])
        uv = np.asarray([direction_vector(a)[1:] for a in anchors])
        ax.add_patch(plt.Rectangle((-45, -30), 90, 60, fill=False, color="0.3"))
        ax.scatter(xy[:, 0], xy[:, 1], s=36, color="#187a72")
        # Draw y/z direction as a tiny in-plane cue: upward-facing sites are dots
        # with short marks, while side-biased sites show longer lateral marks.
        ax.quiver(
            xy[:, 0],
            xy[:, 1],
            np.zeros(len(anchors)),
            uv[:, 0],
            angles="xy",
            scale_units="xy",
            scale=0.35,
            width=0.004,
            color="#c45a2c",
            alpha=0.75,
        )
        ax.set_title(layout_name.replace("_orientation_optimized", ""))
        ax.set_xlim(-50, 50)
        ax.set_ylim(-35, 35)
        ax.set_aspect("equal")
        ax.grid(alpha=0.2)
    fig.savefig(OUT_LAYOUT_PLOT, dpi=220)
    plt.close(fig)


def main() -> None:
    rng = np.random.default_rng(SEED)
    linker_models = load_linker_models()
    layouts = read_orientation_optimized_layouts()
    rows: list[dict[str, str]] = []

    for layout_name, anchors in layouts.items():
        up_count = sum(str(anchor["orientation_class"]) == "up_facing" for anchor in anchors)
        for receptor_count in RECEPTOR_COUNTS:
            metrics = score_layout_3d(anchors, receptor_count, linker_models, rng)
            rows.append(
                {
                    "layout": layout_name,
                    "aptamer_count": str(len(anchors)),
                    "up_facing_aptamers": str(up_count),
                    "ev_diameter_nm": f"{EV_DIAMETER_NM:.0f}",
                    "receptor_count": str(receptor_count),
                    "mean_contacts": f"{metrics['mean_contacts']:.4f}",
                    "max_possible_contacts": f"{metrics['max_possible_contacts']:.0f}",
                    "p_at_least_1_contact": f"{metrics['p_at_least_1_contact']:.4f}",
                    "p_at_least_2_contacts": f"{metrics['p_at_least_2_contacts']:.4f}",
                    "p_at_least_3_contacts": f"{metrics['p_at_least_3_contacts']:.4f}",
                    "clinical_sparse_3d_score": f"{metrics['clinical_sparse_3d_score']:.4f}",
                }
            )

    write_csv(OUT_CSV, rows)
    plot_scores(rows)
    plot_layouts(layouts)

    best_by_count = {}
    for count in RECEPTOR_COUNTS:
        subset = [row for row in rows if row["receptor_count"] == str(count)]
        best_by_count[str(count)] = max(
            subset,
            key=lambda row: float(row["clinical_sparse_3d_score"]),
        )

    layout_summary = []
    for layout_name in sorted(layouts):
        subset = [row for row in rows if row["layout"] == layout_name]
        scores = np.asarray([float(row["clinical_sparse_3d_score"]) for row in subset])
        p1 = np.asarray([float(row["p_at_least_1_contact"]) for row in subset])
        layout_summary.append(
            {
                "layout": layout_name,
                "mean_3d_sparse_score": f"{float(np.mean(scores)):.4f}",
                "mean_p_at_least_1_contact": f"{float(np.mean(p1)):.4f}",
                "up_facing_aptamers": subset[0]["up_facing_aptamers"],
                "aptamer_count": subset[0]["aptamer_count"],
            }
        )
    layout_summary.sort(key=lambda row: float(row["mean_3d_sparse_score"]), reverse=True)

    report = {
        "model": "simplified 3D oriented EV capture screen",
        "plain_language_summary": [
            "Each aptamer has an x/y/z base position and a 3D direction vector.",
            "Each receptor is a point on the lower half of a 73 nm spherical EV.",
            "A contact only counts strongly when the receptor is close enough and inside the aptamer reach cone.",
        ],
        "geometry_assumptions": {
            "ev_diameter_nm": EV_DIAMETER_NM,
            "aptamer_base_height_nm": APTAMER_BASE_HEIGHT_NM,
            "reach_cone_half_angle_degrees": REACH_CONE_HALF_ANGLE_DEG,
            "soft_outer_angle_degrees": SOFT_OUTER_ANGLE_DEG,
            "receptor_counts": RECEPTOR_COUNTS,
            "lateral_scan_nm": LATERAL_SCAN_NM,
            "lateral_step_nm": LATERAL_STEP_NM,
            "receptor_realizations_per_offset": N_RECEPTOR_REALIZATIONS,
            "binding_trials_per_realization": N_BINDING_TRIALS,
        },
        "best_by_receptor_count": best_by_count,
        "layout_summary": layout_summary,
        "best_overall_layout": layout_summary[0],
        "limitations": [
            "This is a coarse 3D reach-cone model, not full oxDNA molecular dynamics.",
            "The aptamer is represented as a flexible cone, not an explicit folded sequence.",
            "The EV is treated as a rigid sphere with receptor points on the lower hemisphere.",
        ],
        "outputs": {
            "scores": OUT_CSV.name,
            "summary": OUT_JSON.name,
            "score_plot": OUT_PLOT.name,
            "layout_plot": OUT_LAYOUT_PLOT.name,
        },
    }
    OUT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="ascii")

    print(f"Wrote {OUT_CSV.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_PLOT.name}")
    print(f"Wrote {OUT_LAYOUT_PLOT.name}")
    print(
        f"Best 3D layout: {layout_summary[0]['layout']} "
        f"mean_3d_sparse={layout_summary[0]['mean_3d_sparse_score']} "
        f"p1={layout_summary[0]['mean_p_at_least_1_contact']}"
    )


if __name__ == "__main__":
    main()
