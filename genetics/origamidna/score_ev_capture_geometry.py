#!/usr/bin/env python3
"""Score probabilistic multivalent CD133+ EV capture by origami aptamer layouts."""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Union

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
IN_CSV = ROOT / "ev_origami_aptamer_layouts.csv"
LINKER_MODEL_CSV = ROOT / "linker_reach_models.csv"
OUT_CSV = ROOT / "ev_capture_scores.csv"
OUT_JSON = ROOT / "ev_capture_summary.json"
OUT_HEATMAP = ROOT / "ev_capture_score_heatmap.png"
OUT_PROB = ROOT / "ev_capture_multivalent_probability.png"

EV_DIAMETERS_NM = [50.0, 100.0, 150.0]
CD133_DENSITIES = {
    # Receptors are modeled as finite stochastic sites on the lower EV hemisphere.
    # Values are nominal CD133 sites per 1000 nm^2 of lower-hemisphere surface.
    "low": 1.5,
    "medium": 4.0,
    "high": 8.0,
}

CONTACT_THRESHOLD = 3
STRONG_CONTACT_THRESHOLD = 6
EV_SURFACE_CLEARANCE_NM = 2.0
LATERAL_SCAN_NM = 24.0
LATERAL_STEP_NM = 4.0
N_RECEPTOR_REALIZATIONS = 8
N_BINDING_TRIALS = 16
REACH_EDGE_THRESHOLD = 0.05
RNG_SEED = 20260525

Anchor = dict[str, Union[float, str]]
LinkerModels = dict[str, tuple[np.ndarray, np.ndarray]]


def load_layouts() -> dict[str, list[Anchor]]:
    layouts: dict[str, list[Anchor]] = defaultdict(list)
    with open(IN_CSV, newline="", encoding="ascii") as f:
        for row in csv.DictReader(f):
            layouts[row["layout"]].append(
                {
                    "x_nm": float(row["x_nm"]),
                    "y_nm": float(row["y_nm"]),
                    "linker_reach_nm": float(row["linker_reach_nm"]),
                    "linker_construct": row.get("linker_construct", "polyT20"),
                }
            )
    return dict(layouts)


def load_linker_models() -> LinkerModels:
    if not LINKER_MODEL_CSV.exists():
        raise FileNotFoundError(
            f"{LINKER_MODEL_CSV.name} is missing. Run python3 calibrate_linker_reach.py first."
        )

    models: LinkerModels = {}
    with open(LINKER_MODEL_CSV, newline="", encoding="ascii") as f:
        for row in csv.DictReader(f):
            xs = []
            ys = []
            for key, value in row.items():
                if key.startswith("p_reach_") and key.endswith("_nm"):
                    xs.append(float(key.removeprefix("p_reach_").removesuffix("_nm")))
                    ys.append(float(value))
            order = np.argsort(xs)
            models[row["construct"]] = (np.asarray(xs)[order], np.asarray(ys)[order])
    return models


def reach_probability(models: LinkerModels, construct: str, distance_nm: float) -> float:
    xs, ys = models[construct]
    if distance_nm <= 0:
        return 1.0
    if distance_nm > xs[-1]:
        return 0.0
    return float(np.interp(distance_nm, xs, ys))


def receptor_count(ev_radius_nm: float, density_per_1000_nm2: float) -> int:
    lower_hemisphere_area_nm2 = 2.0 * math.pi * ev_radius_nm**2
    return max(1, int(round(lower_hemisphere_area_nm2 * density_per_1000_nm2 / 1000.0)))


def random_lower_hemisphere(n: int, radius_nm: float, rng: np.random.Generator) -> np.ndarray:
    theta = rng.uniform(0.0, 2.0 * math.pi, size=n)
    z_unit = -rng.uniform(0.0, 1.0, size=n)
    xy_unit = np.sqrt(np.clip(1.0 - z_unit * z_unit, 0.0, None))
    return np.column_stack(
        (
            radius_nm * xy_unit * np.cos(theta),
            radius_nm * xy_unit * np.sin(theta),
            radius_nm * z_unit,
        )
    )


def max_bipartite_matches(edges: np.ndarray) -> int:
    if edges.size == 0:
        return 0

    n_anchors, n_receptors = edges.shape
    match_to_anchor = np.full(n_receptors, -1, dtype=int)
    adjacency = [np.flatnonzero(edges[i]) for i in range(n_anchors)]

    def assign(anchor_index: int, seen: np.ndarray) -> bool:
        for receptor_index in adjacency[anchor_index]:
            if seen[receptor_index]:
                continue
            seen[receptor_index] = True
            if match_to_anchor[receptor_index] == -1 or assign(match_to_anchor[receptor_index], seen):
                match_to_anchor[receptor_index] = anchor_index
                return True
        return False

    matches = 0
    # Anchors with fewer candidate receptors are harder to place, so try them first.
    for anchor_index in sorted(range(n_anchors), key=lambda i: len(adjacency[i])):
        if assign(anchor_index, np.zeros(n_receptors, dtype=bool)):
            matches += 1
    return matches


def contact_metrics(
    anchors: list[Anchor],
    linker_models: LinkerModels,
    receptor_points: np.ndarray,
    ev_radius_nm: float,
    offset_x_nm: float,
    offset_y_nm: float,
    rng: np.random.Generator,
) -> tuple[float, float, float, float, float, float]:
    center = np.array([offset_x_nm, offset_y_nm, ev_radius_nm + EV_SURFACE_CLEARANCE_NM])
    receptors = receptor_points + center

    probabilities = np.zeros((len(anchors), len(receptors)), dtype=float)
    for anchor_index, anchor in enumerate(anchors):
        anchor_xyz = np.array([float(anchor["x_nm"]), float(anchor["y_nm"]), 0.0])
        distances = np.linalg.norm(receptors - anchor_xyz, axis=1)
        probabilities[anchor_index, :] = [
            reach_probability(linker_models, str(anchor["linker_construct"]), float(distance))
            for distance in distances
        ]

    possible_contacts = float(max_bipartite_matches(probabilities > REACH_EDGE_THRESHOLD))
    sampled_contacts = np.empty(N_BINDING_TRIALS, dtype=float)
    for trial in range(N_BINDING_TRIALS):
        sampled_contacts[trial] = max_bipartite_matches(rng.random(probabilities.shape) < probabilities)

    expected_contacts = float(np.mean(sampled_contacts))
    p_at_least_1 = float(np.mean(sampled_contacts >= 1))
    p_at_least_2 = float(np.mean(sampled_contacts >= 2))
    p_at_least_3 = float(np.mean(sampled_contacts >= CONTACT_THRESHOLD))
    p_at_least_6 = float(np.mean(sampled_contacts >= STRONG_CONTACT_THRESHOLD))
    return expected_contacts, possible_contacts, p_at_least_1, p_at_least_2, p_at_least_3, p_at_least_6


def score_layout(
    anchors: list[Anchor],
    linker_models: LinkerModels,
    ev_diameter_nm: float,
    density_per_1000_nm2: float,
    rng: np.random.Generator,
    fixed_receptor_count: int | None = None,
) -> dict[str, float]:
    ev_radius_nm = ev_diameter_nm / 2.0
    n_receptors = (
        max(1, int(fixed_receptor_count))
        if fixed_receptor_count is not None
        else receptor_count(ev_radius_nm, density_per_1000_nm2)
    )

    offsets = np.arange(-LATERAL_SCAN_NM, LATERAL_SCAN_NM + 0.001, LATERAL_STEP_NM)
    expected_contacts = []
    possible_contacts = []
    p1_values = []
    p2_values = []
    p3_values = []
    p6_values = []
    for ox in offsets:
        for oy in offsets:
            for _ in range(N_RECEPTOR_REALIZATIONS):
                receptor_points = random_lower_hemisphere(n_receptors, ev_radius_nm, rng)
                expected, possible, p1, p2, p3, p6 = contact_metrics(
                    anchors, linker_models, receptor_points, ev_radius_nm, ox, oy, rng
                )
                expected_contacts.append(expected)
                possible_contacts.append(possible)
                p1_values.append(p1)
                p2_values.append(p2)
                p3_values.append(p3)
                p6_values.append(p6)

    expected_arr = np.asarray(expected_contacts, dtype=float)
    possible_arr = np.asarray(possible_contacts, dtype=float)
    p3 = float(np.mean(p3_values))
    p6 = float(np.mean(p6_values))
    p1 = float(np.mean(p1_values))
    p2 = float(np.mean(p2_values))
    mean_contacts = float(np.mean(expected_arr))
    max_contacts = float(np.max(possible_arr))
    normalized_contacts = min(mean_contacts / 8.0, 1.0)
    capture_score = 0.45 * p3 + 0.35 * p6 + 0.20 * normalized_contacts

    return {
        "receptor_count": float(n_receptors),
        "mean_contacts": mean_contacts,
        "max_contacts": max_contacts,
        "p_at_least_1_contact": p1,
        "p_at_least_2_contacts": p2,
        "p_at_least_3_contacts": p3,
        "p_at_least_6_contacts": p6,
        "capture_score": capture_score,
    }


def write_scores(rows: list[dict[str, str]]) -> None:
    with open(OUT_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "layout",
                "aptamer_count",
                "linker_constructs",
                "ev_diameter_nm",
                "cd133_density",
                "receptor_count",
                "mean_contacts",
                "max_contacts",
                "p_at_least_3_contacts",
                "p_at_least_6_contacts",
                "capture_score",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def plot_heatmap(rows: list[dict[str, str]]) -> None:
    layouts = sorted({r["layout"] for r in rows})
    densities = ["low", "medium", "high"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True, sharey=True)
    for ax, diameter in zip(axes, EV_DIAMETERS_NM):
        mat = np.zeros((len(layouts), len(densities)))
        for i, layout in enumerate(layouts):
            for j, density in enumerate(densities):
                match = [
                    r
                    for r in rows
                    if r["layout"] == layout
                    and float(r["ev_diameter_nm"]) == diameter
                    and r["cd133_density"] == density
                ][0]
                mat[i, j] = float(match["capture_score"])
        im = ax.imshow(mat, vmin=0, vmax=1, cmap="viridis")
        ax.set_title(f"EV diameter: {diameter:.0f} nm")
        ax.set_xticks(range(len(densities)))
        ax.set_xticklabels(densities)
        ax.set_yticks(range(len(layouts)))
        ax.set_yticklabels(layouts)
        for i in range(len(layouts)):
            for j in range(len(densities)):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", color="white", fontsize=8)
    fig.colorbar(im, ax=axes, label="capture score")
    fig.suptitle("CD133+ EV capture score with probabilistic tether reach")
    fig.savefig(OUT_HEATMAP, dpi=220)
    plt.close(fig)


def plot_probability(rows: list[dict[str, str]]) -> None:
    layouts = sorted({r["layout"] for r in rows})
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for layout in layouts:
        xs = []
        ys = []
        for diameter in EV_DIAMETERS_NM:
            match = [
                r
                for r in rows
                if r["layout"] == layout
                and float(r["ev_diameter_nm"]) == diameter
                and r["cd133_density"] == "medium"
            ][0]
            xs.append(diameter)
            ys.append(float(match["p_at_least_6_contacts"]))
        ax.plot(xs, ys, marker="o", label=layout)
    ax.set_xlabel("EV diameter (nm)")
    ax.set_ylabel("P(at least 6 simultaneous contacts)")
    ax.set_ylim(-0.03, 1.03)
    ax.grid(alpha=0.25)
    ax.legend(title="layout")
    ax.set_title("Medium CD133 density strong multivalent-capture probability")
    fig.savefig(OUT_PROB, dpi=220)
    plt.close(fig)


def main() -> None:
    layouts = load_layouts()
    linker_models = load_linker_models()
    rng = np.random.default_rng(RNG_SEED)
    rows: list[dict[str, str]] = []
    for layout, anchors in sorted(layouts.items()):
        for diameter in EV_DIAMETERS_NM:
            for density_name, density_per_1000_nm2 in CD133_DENSITIES.items():
                metrics = score_layout(anchors, linker_models, diameter, density_per_1000_nm2, rng)
                rows.append(
                    {
                        "layout": layout,
                        "aptamer_count": str(len(anchors)),
                        "linker_constructs": ";".join(sorted({str(a["linker_construct"]) for a in anchors})),
                        "ev_diameter_nm": f"{diameter:.0f}",
                        "cd133_density": density_name,
                        "receptor_count": f"{metrics['receptor_count']:.0f}",
                        "mean_contacts": f"{metrics['mean_contacts']:.3f}",
                        "max_contacts": f"{metrics['max_contacts']:.0f}",
                        "p_at_least_3_contacts": f"{metrics['p_at_least_3_contacts']:.4f}",
                        "p_at_least_6_contacts": f"{metrics['p_at_least_6_contacts']:.4f}",
                        "capture_score": f"{metrics['capture_score']:.4f}",
                    }
                )

    write_scores(rows)
    plot_heatmap(rows)
    plot_probability(rows)

    best = max(rows, key=lambda r: float(r["capture_score"]))
    summary = {
        "model": "coarse EV surface plus probabilistic aptamer-tether reach",
        "linker_model_file": LINKER_MODEL_CSV.name,
        "contact_threshold": CONTACT_THRESHOLD,
        "strong_contact_threshold": STRONG_CONTACT_THRESHOLD,
        "surface_clearance_nm": EV_SURFACE_CLEARANCE_NM,
        "cd133_density_sites_per_1000_nm2": CD133_DENSITIES,
        "receptor_realizations": N_RECEPTOR_REALIZATIONS,
        "binding_trials_per_realization": N_BINDING_TRIALS,
        "receptor_occupancy": "finite one-to-one aptamer/CD133 matching",
        "rng_seed": RNG_SEED,
        "best_overall": best,
        "notes": [
            "Scores are not binding free energies.",
            "CD133 sites are randomly sampled on the lower EV hemisphere and cannot be reused by multiple aptamers in a single binding trial.",
            "Current reach model uses oxDNA trajectories when present, otherwise WLC fallback.",
            "Use LAMMPS or HOOMD-blue next for dynamic EV capture.",
        ],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")

    print(f"Wrote {OUT_CSV.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_HEATMAP.name}")
    print(f"Wrote {OUT_PROB.name}")
    print("Best overall:")
    print(
        f"  {best['layout']} EV={best['ev_diameter_nm']} nm "
        f"density={best['cd133_density']} score={best['capture_score']}"
    )


if __name__ == "__main__":
    main()
