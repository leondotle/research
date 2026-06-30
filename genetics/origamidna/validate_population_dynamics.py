#!/usr/bin/env python3
"""Brownian validation for top population-optimized layouts."""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Union

import numpy as np

from score_ev_capture_geometry import load_linker_models, reach_probability

ROOT = Path(__file__).resolve().parent
IN_POPULATION = ROOT / "ev_population_optimization_run.npz"
IN_LAYOUTS = ROOT / "population_optimized_layouts.csv"
IN_SCORES = ROOT / "population_layout_scores.csv"
OUT_CSV = ROOT / "population_layout_dynamics_validation.csv"
OUT_JSON = ROOT / "population_layout_dynamics_validation_summary.json"

RNG_SEED = 20260617
TOP_LAYOUTS = 6
N_STEPS = 360
DT_SECONDS = 0.05
CAPTURE_DWELL_SECONDS = 1.0
CAPTURE_THRESHOLD = 1
STRONG_THRESHOLD = 2
SURFACE_CLEARANCE_NM = 2.0
INITIAL_GAP_NM = 8.0
MAX_GAP_NM = 22.0
LATERAL_START_NM = 22.0
LATERAL_ESCAPE_NM = 70.0
D_FREE_NM2_PER_S = 35.0
D_BOUND_FLOOR_FRACTION = 0.10
K_ON_PER_STEP = 0.18
K_OFF_PER_S = 0.12

Anchor = dict[str, Union[float, str]]


def load_top_layout_names() -> list[str]:
    with open(IN_SCORES, newline="", encoding="ascii") as f:
        rows = list(csv.DictReader(f))
    return [row["layout"] for row in rows[:TOP_LAYOUTS]]


def load_layouts(names: list[str]) -> dict[str, list[Anchor]]:
    layouts: dict[str, list[Anchor]] = defaultdict(list)
    with open(IN_LAYOUTS, newline="", encoding="ascii") as f:
        for row in csv.DictReader(f):
            if row["layout"] not in names:
                continue
            layouts[row["layout"]].append(
                {
                    "x_nm": float(row["x_nm"]),
                    "y_nm": float(row["y_nm"]),
                    "linker_construct": row["linker_construct"],
                    "linker_reach_nm": 15.0,
                }
            )
    return dict(layouts)


def anchor_array(anchors: list[Anchor]) -> np.ndarray:
    return np.asarray(
        [[float(anchor["x_nm"]), float(anchor["y_nm"]), 0.0] for anchor in anchors],
        dtype=float,
    )


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


def simulate_one(
    anchors: list[Anchor],
    linker_models,
    receptor_body: np.ndarray,
    diameter_nm: float,
    rng: np.random.Generator,
) -> dict[str, float]:
    radius = diameter_nm / 2.0
    anchor_xyz = anchor_array(anchors)
    center = np.array(
        [
            rng.uniform(-LATERAL_START_NM, LATERAL_START_NM),
            rng.uniform(-LATERAL_START_NM, LATERAL_START_NM),
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
        lateral_norm = np.linalg.norm(center[:2])
        if lateral_norm > LATERAL_ESCAPE_NM:
            center[:2] *= LATERAL_ESCAPE_NM / lateral_norm
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
        "mean_max_contacts": float(np.max(contacts_arr)),
        "strong_fraction": float(np.mean(contacts_arr >= STRONG_THRESHOLD)),
    }


def main() -> None:
    names = load_top_layout_names()
    layouts = load_layouts(names)
    population = np.load(IN_POPULATION)
    linker_models = load_linker_models()
    rows = []
    for name in names:
        rng = np.random.default_rng(RNG_SEED + len(rows))
        anchors = layouts[name]
        per_ev = []
        for ev_index, diameter in enumerate(population["diameter_nm"]):
            receptor_count = int(population["receptor_count"][ev_index])
            receptor_body = population["receptor_points"][ev_index, :receptor_count, :]
            per_ev.append(simulate_one(anchors, linker_models, receptor_body, float(diameter), rng))
        rows.append(
            {
                "layout": name,
                "aptamer_count": str(len(anchors)),
                "population_evs": str(len(per_ev)),
                "capture_probability_p_ge_1_dwell": f"{np.mean([r['ever_captured'] for r in per_ev]):.4f}",
                "mean_contacts": f"{np.mean([r['mean_contacts'] for r in per_ev]):.4f}",
                "mean_max_contacts": f"{np.mean([r['mean_max_contacts'] for r in per_ev]):.4f}",
                "strong_fraction_p_ge_2": f"{np.mean([r['strong_fraction'] for r in per_ev]):.4f}",
            }
        )
        print(f"Validated {name}: capture={rows[-1]['capture_probability_p_ge_1_dwell']}", flush=True)
    rows.sort(key=lambda row: float(row["capture_probability_p_ge_1_dwell"]), reverse=True)
    with open(OUT_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "model": "Brownian validation of top population-optimized layouts",
        "rng_seed": RNG_SEED,
        "top_layouts_validated": TOP_LAYOUTS,
        "steps": N_STEPS,
        "dt_seconds": DT_SECONDS,
        "capture_definition": "at least one contact for 1 second",
        "best_layout": rows[0],
        "summary_rows": rows,
        "outputs": {"summary_csv": OUT_CSV.name},
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
    print(f"Wrote {OUT_CSV.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(
        f"Best dynamic validation: {rows[0]['layout']} "
        f"capture={rows[0]['capture_probability_p_ge_1_dwell']}"
    )


if __name__ == "__main__":
    main()
