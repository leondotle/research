#!/usr/bin/env python3
"""Simulate repeated DNA-origami capture tiles on a surface.

Beginner picture:
One sticky tile is one trap. A capture field is many reasonable traps repeated
on a surface. This tests whether repeated chances can beat one over-crowded
80-aptamer tile.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Union

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from clinical_layouts import clinical_candidate_layouts
from ev_population_generator import generate_population
from score_ev_capture_geometry import load_layouts, load_linker_models, reach_probability

ROOT = Path(__file__).resolve().parent
POPULATION_FILE = ROOT / "ev_population_optimization_run.npz"
OUT_SUMMARY_CSV = ROOT / "multitile_capture_field_summary.csv"
OUT_LAYOUTS_CSV = ROOT / "multitile_capture_field_layouts.csv"
OUT_JSON = ROOT / "multitile_capture_field_summary.json"
OUT_PLOT = ROOT / "multitile_capture_field_scores.png"
OUT_LAYOUT_PNG = ROOT / "multitile_capture_field_layouts.png"

RNG_SEED = 20260619
LINKER = "polyT30"
N_POPULATION_EVS = 160
N_STEPS = 300
DT_SECONDS = 0.05
CAPTURE_DWELL_SECONDS = 1.0
CAPTURE_THRESHOLD = 1
STRONG_THRESHOLD = 2
SURFACE_CLEARANCE_NM = 2.0
INITIAL_GAP_NM = 8.0
MAX_GAP_NM = 24.0
D_FREE_NM2_PER_S = 35.0
D_BOUND_FLOOR_FRACTION = 0.10
K_ON_PER_STEP = 0.18
K_OFF_PER_S = 0.12

Anchor = dict[str, Union[float, str]]


def with_linker(anchors: list[Anchor]) -> list[Anchor]:
    return [
        {
            **anchor,
            "linker_construct": LINKER,
            "linker_reach_nm": float(anchor.get("linker_reach_nm", 15.0)),
        }
        for anchor in anchors
    ]


def base_layouts() -> dict[str, list[Anchor]]:
    layouts = load_layouts()
    layouts.update(clinical_candidate_layouts(LINKER))
    return {
        "clinical_grid_12": with_linker(layouts["clinical_grid_12"]),
        "clinical_grid_18": with_linker(layouts["clinical_grid_18"]),
        "clinical_dense_24": with_linker(layouts["clinical_dense_24"]),
    }


def build_field(
    tile_anchors: list[Anchor],
    cols: int,
    rows: int,
    spacing_x_nm: float,
    spacing_y_nm: float,
) -> list[Anchor]:
    anchors: list[Anchor] = []
    x_offsets = np.linspace(-spacing_x_nm * (cols - 1) / 2.0, spacing_x_nm * (cols - 1) / 2.0, cols)
    y_offsets = np.linspace(-spacing_y_nm * (rows - 1) / 2.0, spacing_y_nm * (rows - 1) / 2.0, rows)
    for ox in x_offsets:
        for oy in y_offsets:
            for anchor in tile_anchors:
                anchors.append(
                    {
                        "x_nm": float(anchor["x_nm"]) + float(ox),
                        "y_nm": float(anchor["y_nm"]) + float(oy),
                        "linker_reach_nm": 15.0,
                        "linker_construct": LINKER,
                    }
                )
    return anchors


def anchor_array(anchors: list[Anchor]) -> np.ndarray:
    return np.asarray(
        [[float(anchor["x_nm"]), float(anchor["y_nm"]), 0.0] for anchor in anchors],
        dtype=float,
    )


def field_bounds(cols: int, rows: int, spacing_x_nm: float, spacing_y_nm: float) -> tuple[float, float]:
    width = spacing_x_nm * max(cols - 1, 0) + 90.0
    height = spacing_y_nm * max(rows - 1, 0) + 60.0
    return width, height


def break_bonds(
    distances: np.ndarray,
    anchor_to_receptor: np.ndarray,
    receptor_to_anchor: np.ndarray,
    rng: np.random.Generator,
) -> None:
    p_off = 1.0 - math.exp(-K_OFF_PER_S * DT_SECONDS)
    for anchor_index, receptor_index in enumerate(anchor_to_receptor):
        if receptor_index == -1:
            continue
        if distances[anchor_index, receptor_index] > 14.0 or rng.random() < p_off:
            anchor_to_receptor[anchor_index] = -1
            receptor_to_anchor[receptor_index] = -1


def try_form_bonds(
    distances: np.ndarray,
    anchors: list[Anchor],
    linker_models,
    anchor_to_receptor: np.ndarray,
    receptor_to_anchor: np.ndarray,
    rng: np.random.Generator,
) -> None:
    free_anchor_indices = np.flatnonzero(anchor_to_receptor == -1)
    free_receptor_indices = np.flatnonzero(receptor_to_anchor == -1)
    if len(free_anchor_indices) == 0 or len(free_receptor_indices) == 0:
        return
    sub_distances = distances[np.ix_(free_anchor_indices, free_receptor_indices)]
    probabilities = np.zeros_like(sub_distances)
    for local_i, anchor_index in enumerate(free_anchor_indices):
        construct = str(anchors[int(anchor_index)]["linker_construct"])
        probabilities[local_i, :] = [
            K_ON_PER_STEP * reach_probability(linker_models, construct, float(distance))
            for distance in sub_distances[local_i, :]
        ]
    candidate_pairs = np.argwhere(rng.random(probabilities.shape) < probabilities)
    if len(candidate_pairs) == 0:
        return
    rng.shuffle(candidate_pairs)
    for local_anchor, local_receptor in candidate_pairs:
        anchor_index = int(free_anchor_indices[local_anchor])
        receptor_index = int(free_receptor_indices[local_receptor])
        if anchor_to_receptor[anchor_index] == -1 and receptor_to_anchor[receptor_index] == -1:
            anchor_to_receptor[anchor_index] = receptor_index
            receptor_to_anchor[receptor_index] = anchor_index


def simulate_one_ev(
    anchors: list[Anchor],
    linker_models,
    receptor_body: np.ndarray,
    diameter_nm: float,
    field_width_nm: float,
    field_height_nm: float,
    rng: np.random.Generator,
) -> dict[str, float]:
    radius = diameter_nm / 2.0
    anchor_xyz = anchor_array(anchors)
    # Start somewhere over the field. This models an EV that has entered the
    # capture region, not the whole journey from a cell.
    center = np.array(
        [
            rng.uniform(-field_width_nm / 2.0, field_width_nm / 2.0),
            rng.uniform(-field_height_nm / 2.0, field_height_nm / 2.0),
            radius + SURFACE_CLEARANCE_NM + INITIAL_GAP_NM,
        ],
        dtype=float,
    )
    min_z = radius + SURFACE_CLEARANCE_NM
    max_z = radius + SURFACE_CLEARANCE_NM + MAX_GAP_NM
    anchor_to_receptor = np.full(len(anchors), -1, dtype=int)
    receptor_to_anchor = np.full(len(receptor_body), -1, dtype=int)
    required_steps = max(1, int(round(CAPTURE_DWELL_SECONDS / DT_SECONDS)))
    consecutive = 0
    ever_captured = False
    contacts = []

    for _ in range(N_STEPS):
        n_contacts = int(np.count_nonzero(anchor_to_receptor != -1))
        mobility = max(D_BOUND_FLOOR_FRACTION, 1.0 / (1.0 + 0.75 * n_contacts))
        sigma = math.sqrt(2.0 * D_FREE_NM2_PER_S * mobility * DT_SECONDS)
        center[:2] += rng.normal(0.0, sigma, size=2)
        center[2] += rng.normal(0.0, sigma * 0.45)
        center[2] -= 0.035 * n_contacts
        center[2] = float(np.clip(center[2], min_z, max_z))
        center[0] = float(np.clip(center[0], -field_width_nm / 2.0, field_width_nm / 2.0))
        center[1] = float(np.clip(center[1], -field_height_nm / 2.0, field_height_nm / 2.0))

        receptors = receptor_body + center
        distances = np.linalg.norm(anchor_xyz[:, None, :] - receptors[None, :, :], axis=2)
        break_bonds(distances, anchor_to_receptor, receptor_to_anchor, rng)
        try_form_bonds(distances, anchors, linker_models, anchor_to_receptor, receptor_to_anchor, rng)

        n_contacts = int(np.count_nonzero(anchor_to_receptor != -1))
        contacts.append(n_contacts)
        consecutive = consecutive + 1 if n_contacts >= CAPTURE_THRESHOLD else 0
        ever_captured = ever_captured or consecutive >= required_steps

    contacts_arr = np.asarray(contacts, dtype=float)
    return {
        "ever_captured": float(ever_captured),
        "mean_contacts": float(np.mean(contacts_arr)),
        "max_contacts": float(np.max(contacts_arr)),
        "strong_fraction": float(np.mean(contacts_arr >= STRONG_THRESHOLD)),
    }


def ensure_population() -> np.lib.npyio.NpzFile:
    if not POPULATION_FILE.exists():
        records, receptor_array = generate_population(N_POPULATION_EVS, RNG_SEED)
        np.savez_compressed(
            POPULATION_FILE,
            receptor_points=receptor_array,
            diameter_nm=np.asarray([r.diameter_nm for r in records], dtype=float),
            receptor_count=np.asarray([r.receptor_count for r in records], dtype=int),
            pattern=np.asarray([r.pattern for r in records]),
        )
    return np.load(POPULATION_FILE)


def field_configs() -> list[dict[str, object]]:
    configs = []
    for unit in ("clinical_grid_18", "clinical_dense_24"):
        for cols, rows in ((1, 1), (3, 3), (4, 3), (5, 3), (5, 4)):
            for spacing in (70.0, 90.0):
                configs.append(
                    {
                        "unit_layout": unit,
                        "cols": cols,
                        "rows": rows,
                        "spacing_nm": spacing,
                    }
                )
    return configs


def write_layout_preview(configs: list[dict[str, object]], layouts: dict[str, list[Anchor]]) -> None:
    top = configs[:6]
    with open(OUT_LAYOUTS_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["field", "anchor_id", "x_nm", "y_nm", "linker_construct"],
        )
        writer.writeheader()
        for config in top:
            name = str(config["field"])
            for i, anchor in enumerate(layouts[name], start=1):
                writer.writerow(
                    {
                        "field": name,
                        "anchor_id": f"A{i:03d}",
                        "x_nm": f"{float(anchor['x_nm']):.3f}",
                        "y_nm": f"{float(anchor['y_nm']):.3f}",
                        "linker_construct": anchor["linker_construct"],
                    }
                )


def plot_scores(rows: list[dict[str, str]]) -> None:
    top = rows[:14]
    labels = [row["field"] for row in top][::-1]
    values = [float(row["capture_probability"]) for row in top][::-1]
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    ax.barh(labels, values, color="#1565c0")
    ax.axvline(0.90, color="crimson", ls="--", lw=1.2, label="90% target")
    ax.set_xlim(0, 1)
    ax.set_xlabel("dynamic capture probability")
    ax.set_title("Repeated-tile capture fields")
    ax.grid(axis="x", alpha=0.25)
    ax.legend()
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)


def plot_layouts(configs: list[dict[str, object]], layouts: dict[str, list[Anchor]]) -> None:
    top = configs[:6]
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for ax, config in zip(axes.flat, top):
        name = str(config["field"])
        anchors = layouts[name]
        xy = np.asarray([[float(a["x_nm"]), float(a["y_nm"])] for a in anchors], dtype=float)
        ax.scatter(xy[:, 0], xy[:, 1], s=12, color="#1565c0")
        ax.set_title(name, fontsize=8)
        ax.set_aspect("equal")
        ax.grid(alpha=0.2)
    fig.savefig(OUT_LAYOUT_PNG, dpi=220)
    plt.close(fig)


def main() -> None:
    rng = np.random.default_rng(RNG_SEED)
    population = ensure_population()
    linker_models = load_linker_models()
    units = base_layouts()

    rows: list[dict[str, str]] = []
    field_layouts: dict[str, list[Anchor]] = {}
    for index, config in enumerate(field_configs()):
        unit_name = str(config["unit_layout"])
        cols = int(config["cols"])
        rows_count = int(config["rows"])
        spacing = float(config["spacing_nm"])
        field = build_field(units[unit_name], cols, rows_count, spacing, spacing * 0.75)
        width, height = field_bounds(cols, rows_count, spacing, spacing * 0.75)
        field_name = f"{unit_name}_{cols}x{rows_count}_{int(spacing)}nm"
        field_layouts[field_name] = field

        per_ev = []
        config_rng = np.random.default_rng(RNG_SEED + index)
        for ev_index, diameter in enumerate(population["diameter_nm"]):
            receptor_count = int(population["receptor_count"][ev_index])
            receptor_body = population["receptor_points"][ev_index, :receptor_count, :]
            per_ev.append(
                simulate_one_ev(field, linker_models, receptor_body, float(diameter), width, height, config_rng)
            )
        rows.append(
            {
                "field": field_name,
                "unit_layout": unit_name,
                "tile_copies": str(cols * rows_count),
                "aptamers_per_tile": str(len(units[unit_name])),
                "total_aptamers": str(len(field)),
                "spacing_nm": f"{spacing:.0f}",
                "field_width_nm": f"{width:.0f}",
                "field_height_nm": f"{height:.0f}",
                "population_evs": str(len(per_ev)),
                "capture_probability": f"{np.mean([r['ever_captured'] for r in per_ev]):.4f}",
                "mean_contacts": f"{np.mean([r['mean_contacts'] for r in per_ev]):.4f}",
                "mean_max_contacts": f"{np.mean([r['max_contacts'] for r in per_ev]):.4f}",
                "strong_fraction_p_ge_2": f"{np.mean([r['strong_fraction'] for r in per_ev]):.4f}",
            }
        )
        print(f"{field_name}: capture={rows[-1]['capture_probability']}", flush=True)

    rows.sort(key=lambda row: float(row["capture_probability"]), reverse=True)
    with open(OUT_SUMMARY_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    top_configs = []
    for row in rows[:6]:
        top_configs.append({"field": row["field"]})
    write_layout_preview(top_configs, field_layouts)
    plot_scores(rows)
    plot_layouts(top_configs, field_layouts)

    summary = {
        "model": "repeated moderate-density DNA-origami capture tile field",
        "rng_seed": RNG_SEED,
        "population_file": POPULATION_FILE.name,
        "capture_definition": "at least one aptamer/CD133 contact sustained for 1 second",
        "best_field": rows[0],
        "top_10_fields": rows[:10],
        "spacing_interpretation": [
            "Increasing field size helps only when the field has enough tile copies to avoid big gaps.",
            "Wider spacing lowers aptamer density and can create miss zones.",
            "Closer spacing improves capture in this model, but real fabrication may need spacing to avoid steric and electrostatic crowding.",
        ],
        "outputs": {
            "summary_csv": OUT_SUMMARY_CSV.name,
            "top_layout_coordinates": OUT_LAYOUTS_CSV.name,
            "score_plot": OUT_PLOT.name,
            "layout_plot": OUT_LAYOUT_PNG.name,
        },
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
    print(f"Wrote {OUT_SUMMARY_CSV.name}")
    print(f"Wrote {OUT_LAYOUTS_CSV.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Best field: {rows[0]['field']} capture={rows[0]['capture_probability']}")


if __name__ == "__main__":
    main()
