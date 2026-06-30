#!/usr/bin/env python3
"""Score matched random-layout controls for the EV capture geometry screen."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Union

import numpy as np

from score_ev_capture_geometry import (
    CD133_DENSITIES,
    EV_DIAMETERS_NM,
    RNG_SEED,
    load_linker_models,
    load_layouts,
    score_layout,
)

ROOT = Path(__file__).resolve().parent
OUT_CSV = ROOT / "ev_capture_random_layout_controls.csv"
OUT_JSON = ROOT / "ev_capture_random_layout_controls_summary.json"

TILE_WIDTH_NM = 90.0
TILE_HEIGHT_NM = 60.0
CONTROL_SEED = RNG_SEED + 101
CONTROL_REPLICATES = 4


def random_control_layout(
    aptamer_count: int, rng: np.random.Generator
) -> list[dict[str, Union[float, str]]]:
    xs = rng.uniform(-TILE_WIDTH_NM / 2.0, TILE_WIDTH_NM / 2.0, size=aptamer_count)
    ys = rng.uniform(-TILE_HEIGHT_NM / 2.0, TILE_HEIGHT_NM / 2.0, size=aptamer_count)
    return [
        {
            "x_nm": float(x),
            "y_nm": float(y),
            "linker_reach_nm": 15.0,
            "linker_construct": "polyT20",
        }
        for x, y in zip(xs, ys)
    ]


def matched_designed_scores() -> dict[tuple[int, float, str], dict[str, str]]:
    with open(ROOT / "ev_capture_scores.csv", newline="", encoding="ascii") as f:
        rows = list(csv.DictReader(f))

    best_by_key: dict[tuple[int, float, str], dict[str, str]] = {}
    for row in rows:
        key = (int(row["aptamer_count"]), float(row["ev_diameter_nm"]), row["cd133_density"])
        if key not in best_by_key or float(row["capture_score"]) > float(best_by_key[key]["capture_score"]):
            best_by_key[key] = row
    return best_by_key


def summarize_replicates(rows: list[dict[str, float]]) -> dict[str, float]:
    return {
        "receptor_count": rows[0]["receptor_count"],
        "mean_contacts": float(np.mean([row["mean_contacts"] for row in rows])),
        "max_contacts": float(np.max([row["max_contacts"] for row in rows])),
        "p_at_least_3_contacts": float(np.mean([row["p_at_least_3_contacts"] for row in rows])),
        "p_at_least_6_contacts": float(np.mean([row["p_at_least_6_contacts"] for row in rows])),
        "capture_score": float(np.mean([row["capture_score"] for row in rows])),
        "capture_score_sd": float(np.std([row["capture_score"] for row in rows], ddof=0)),
    }


def main() -> None:
    linker_models = load_linker_models()
    designed_layouts = load_layouts()
    best_designed = matched_designed_scores()
    rng = np.random.default_rng(CONTROL_SEED)

    control_layouts: dict[int, list[list[dict[str, Union[float, str]]]]] = {}
    for anchors in designed_layouts.values():
        count = len(anchors)
        if count not in control_layouts:
            control_layouts[count] = [
                random_control_layout(count, rng) for _ in range(CONTROL_REPLICATES)
            ]

    rows: list[dict[str, str]] = []
    for aptamer_count, replicate_layouts in sorted(control_layouts.items()):
        for diameter in EV_DIAMETERS_NM:
            for density_name, density_per_1000_nm2 in CD133_DENSITIES.items():
                replicate_metrics = [
                    score_layout(
                        control_anchors,
                        linker_models,
                        diameter,
                        density_per_1000_nm2,
                        rng,
                    )
                    for control_anchors in replicate_layouts
                ]
                metrics = summarize_replicates(replicate_metrics)
                designed = best_designed[(aptamer_count, diameter, density_name)]
                designed_score = float(designed["capture_score"])
                control_score = metrics["capture_score"]
                rows.append(
                    {
                        "control": f"random_matched_{aptamer_count}",
                        "aptamer_count": str(aptamer_count),
                        "replicates": str(CONTROL_REPLICATES),
                        "ev_diameter_nm": f"{diameter:.0f}",
                        "cd133_density": density_name,
                        "receptor_count": f"{metrics['receptor_count']:.0f}",
                        "mean_contacts": f"{metrics['mean_contacts']:.3f}",
                        "max_contacts": f"{metrics['max_contacts']:.0f}",
                        "p_at_least_3_contacts": f"{metrics['p_at_least_3_contacts']:.4f}",
                        "p_at_least_6_contacts": f"{metrics['p_at_least_6_contacts']:.4f}",
                        "capture_score": f"{metrics['capture_score']:.4f}",
                        "capture_score_sd": f"{metrics['capture_score_sd']:.4f}",
                        "best_matched_design": designed["layout"],
                        "best_matched_design_score": f"{designed_score:.4f}",
                        "score_delta_vs_best_matched_design": f"{designed_score - control_score:.4f}",
                    }
                )

    with open(OUT_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    top_control = max(rows, key=lambda row: float(row["capture_score"]))
    summary = {
        "model": "matched random aptamer-position controls on the same 90 x 60 nm tile",
        "control_seed": CONTROL_SEED,
        "replicates_per_aptamer_count": CONTROL_REPLICATES,
        "control_output": OUT_CSV.name,
        "top_control": top_control,
        "notes": [
            "Controls use the same polyT20/A15 reach model and finite-receptor EV scorer as the designed layouts.",
            "Controls randomize anchor positions but do not model scrambled/nonbinding aptamer chemistry.",
            "Positive score deltas mean the best designed layout for the same aptamer count, EV size, and density outscored the matched random control.",
        ],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")

    print(f"Wrote {OUT_CSV.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(
        f"Top random control: {top_control['control']} EV={top_control['ev_diameter_nm']} "
        f"density={top_control['cd133_density']} score={top_control['capture_score']}"
    )


if __name__ == "__main__":
    main()
