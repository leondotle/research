#!/usr/bin/env python3
"""Run next-stage EV capture sensitivity and control analysis.

This script extends the first-pass geometry screen without overwriting the
baseline outputs. It asks three follow-up questions:

1. Does the lead design survive linker-length sensitivity?
2. Does dense_24 still beat matched random 24-anchor layouts with more control
   replicates?
3. What score should a scrambled/nonbinding aptamer chemistry control receive
   under the same reporting schema?
"""

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
    load_layouts,
    load_linker_models,
    score_layout,
)

ROOT = Path(__file__).resolve().parent
OUT_SENSITIVITY_CSV = ROOT / "ev_capture_linker_sensitivity.csv"
OUT_RANDOM_CSV = ROOT / "ev_capture_dense24_random_controls_50rep.csv"
OUT_SCRAMBLED_CSV = ROOT / "ev_capture_scrambled_controls.csv"
OUT_JSON = ROOT / "ev_capture_next_stage_controls_summary.json"

LINKER_CONSTRUCTS = ("polyT15", "polyT20", "polyT30")
LEAD_LAYOUT = "dense_24"
RANDOM_CONTROL_REPLICATES = 50
TILE_WIDTH_NM = 90.0
TILE_HEIGHT_NM = 60.0
CONTROL_SEED = RNG_SEED + 202

Anchor = dict[str, Union[float, str]]


def with_linker(anchors: list[Anchor], construct: str) -> list[Anchor]:
    return [
        {
            **anchor,
            "linker_construct": construct,
            "linker_reach_nm": float(anchor.get("linker_reach_nm", 15.0)),
        }
        for anchor in anchors
    ]


def random_control_layout(
    aptamer_count: int, construct: str, rng: np.random.Generator
) -> list[Anchor]:
    xs = rng.uniform(-TILE_WIDTH_NM / 2.0, TILE_WIDTH_NM / 2.0, size=aptamer_count)
    ys = rng.uniform(-TILE_HEIGHT_NM / 2.0, TILE_HEIGHT_NM / 2.0, size=aptamer_count)
    return [
        {
            "x_nm": float(x),
            "y_nm": float(y),
            "linker_reach_nm": 15.0,
            "linker_construct": construct,
        }
        for x, y in zip(xs, ys)
    ]


def format_metric(value: float, places: int = 4) -> str:
    return f"{value:.{places}f}"


def score_to_row(
    label: str,
    aptamer_count: int,
    construct: str,
    diameter: float,
    density_name: str,
    metrics: dict[str, float],
) -> dict[str, str]:
    return {
        "layout": label,
        "aptamer_count": str(aptamer_count),
        "linker_construct": construct,
        "ev_diameter_nm": f"{diameter:.0f}",
        "cd133_density": density_name,
        "receptor_count": f"{metrics['receptor_count']:.0f}",
        "mean_contacts": format_metric(metrics["mean_contacts"], 3),
        "max_contacts": f"{metrics['max_contacts']:.0f}",
        "p_at_least_3_contacts": format_metric(metrics["p_at_least_3_contacts"]),
        "p_at_least_6_contacts": format_metric(metrics["p_at_least_6_contacts"]),
        "capture_score": format_metric(metrics["capture_score"]),
    }


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with open(path, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


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


def run_linker_sensitivity(
    layouts: dict[str, list[Anchor]],
    linker_models,
    rng: np.random.Generator,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for layout, anchors in sorted(layouts.items()):
        for construct in LINKER_CONSTRUCTS:
            print(f"Scoring sensitivity: {layout} {construct}", flush=True)
            variant = with_linker(anchors, construct)
            for diameter in EV_DIAMETERS_NM:
                for density_name, density_per_1000_nm2 in CD133_DENSITIES.items():
                    metrics = score_layout(
                        variant,
                        linker_models,
                        diameter,
                        density_per_1000_nm2,
                        rng,
                    )
                    rows.append(
                        score_to_row(
                            layout,
                            len(variant),
                            construct,
                            diameter,
                            density_name,
                            metrics,
                        )
                    )
    return rows


def run_dense24_random_controls(
    lead_anchors: list[Anchor],
    linker_models,
    rng: np.random.Generator,
    sensitivity_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    best_dense24_by_case = {
        (row["linker_construct"], row["ev_diameter_nm"], row["cd133_density"]): row
        for row in sensitivity_rows
        if row["layout"] == LEAD_LAYOUT
    }

    rows: list[dict[str, str]] = []
    for construct in LINKER_CONSTRUCTS:
        print(
            f"Building {RANDOM_CONTROL_REPLICATES} dense_24 random controls for {construct}",
            flush=True,
        )
        control_layouts = [
            random_control_layout(len(lead_anchors), construct, rng)
            for _ in range(RANDOM_CONTROL_REPLICATES)
        ]
        for diameter in EV_DIAMETERS_NM:
            for density_name, density_per_1000_nm2 in CD133_DENSITIES.items():
                print(
                    "Scoring random controls: "
                    f"{construct} EV={diameter:.0f} density={density_name}",
                    flush=True,
                )
                replicate_metrics = [
                    score_layout(
                        control_anchors,
                        linker_models,
                        diameter,
                        density_per_1000_nm2,
                        rng,
                    )
                    for control_anchors in control_layouts
                ]
                metrics = summarize_replicates(replicate_metrics)
                designed = best_dense24_by_case[
                    (construct, f"{diameter:.0f}", density_name)
                ]
                designed_score = float(designed["capture_score"])
                control_score = metrics["capture_score"]
                rows.append(
                    {
                        "control": f"random_matched_24_{construct}",
                        "aptamer_count": str(len(lead_anchors)),
                        "linker_construct": construct,
                        "replicates": str(RANDOM_CONTROL_REPLICATES),
                        "ev_diameter_nm": f"{diameter:.0f}",
                        "cd133_density": density_name,
                        "receptor_count": f"{metrics['receptor_count']:.0f}",
                        "mean_contacts": format_metric(metrics["mean_contacts"], 3),
                        "max_contacts": f"{metrics['max_contacts']:.0f}",
                        "p_at_least_3_contacts": format_metric(metrics["p_at_least_3_contacts"]),
                        "p_at_least_6_contacts": format_metric(metrics["p_at_least_6_contacts"]),
                        "capture_score": format_metric(control_score),
                        "capture_score_sd": format_metric(metrics["capture_score_sd"]),
                        "matched_design": LEAD_LAYOUT,
                        "matched_design_score": format_metric(designed_score),
                        "score_delta_vs_matched_design": format_metric(designed_score - control_score),
                    }
                )
    return rows


def run_scrambled_controls(sensitivity_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in sensitivity_rows:
        scrambled = dict(row)
        scrambled["control"] = f"scrambled_{row['layout']}"
        scrambled["binding_activity"] = "0.0"
        scrambled["mean_contacts"] = "0.000"
        scrambled["max_contacts"] = "0"
        scrambled["p_at_least_3_contacts"] = "0.0000"
        scrambled["p_at_least_6_contacts"] = "0.0000"
        scrambled["capture_score"] = "0.0000"
        rows.append(scrambled)
    return rows


def main() -> None:
    layouts = load_layouts()
    linker_models = load_linker_models()
    rng = np.random.default_rng(CONTROL_SEED)

    print("Running linker sensitivity screen", flush=True)
    sensitivity_rows = run_linker_sensitivity(layouts, linker_models, rng)
    write_csv(OUT_SENSITIVITY_CSV, sensitivity_rows)

    print("Running dense_24 matched random controls", flush=True)
    random_rows = run_dense24_random_controls(
        layouts[LEAD_LAYOUT], linker_models, rng, sensitivity_rows
    )
    write_csv(OUT_RANDOM_CSV, random_rows)

    print("Writing scrambled/nonbinding controls", flush=True)
    scrambled_rows = run_scrambled_controls(sensitivity_rows)
    write_csv(OUT_SCRAMBLED_CSV, scrambled_rows)

    best_sensitivity = max(sensitivity_rows, key=lambda row: float(row["capture_score"]))
    top_random = max(random_rows, key=lambda row: float(row["capture_score"]))
    best_dense24_medium = max(
        (
            row
            for row in sensitivity_rows
            if row["layout"] == LEAD_LAYOUT and row["cd133_density"] == "medium"
        ),
        key=lambda row: float(row["capture_score"]),
    )

    summary = {
        "model": "next-stage linker sensitivity plus dense_24 controls",
        "rng_seed": CONTROL_SEED,
        "linker_constructs": list(LINKER_CONSTRUCTS),
        "random_control_replicates": RANDOM_CONTROL_REPLICATES,
        "outputs": {
            "linker_sensitivity": OUT_SENSITIVITY_CSV.name,
            "dense24_random_controls": OUT_RANDOM_CSV.name,
            "scrambled_controls": OUT_SCRAMBLED_CSV.name,
        },
        "best_overall_sensitivity_case": best_sensitivity,
        "best_dense24_medium_density_case": best_dense24_medium,
        "top_random_control": top_random,
        "interpretation": [
            "Linker sensitivity keeps aptamer count and coordinates fixed while swapping all anchors to one linker construct.",
            "Matched random controls are focused on dense_24 because it is the current lead design and the fabrication go/no-go candidate.",
            "Scrambled/nonbinding controls are represented as zero receptor-specific binding activity; they are report-schema controls, not geometry controls.",
        ],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")

    print(f"Wrote {OUT_SENSITIVITY_CSV.name}")
    print(f"Wrote {OUT_RANDOM_CSV.name}")
    print(f"Wrote {OUT_SCRAMBLED_CSV.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(
        "Best sensitivity case: "
        f"{best_sensitivity['layout']} {best_sensitivity['linker_construct']} "
        f"EV={best_sensitivity['ev_diameter_nm']} density={best_sensitivity['cd133_density']} "
        f"score={best_sensitivity['capture_score']}"
    )
    print(
        "Top random control: "
        f"{top_random['control']} EV={top_random['ev_diameter_nm']} "
        f"density={top_random['cd133_density']} score={top_random['capture_score']}"
    )


if __name__ == "__main__":
    main()
