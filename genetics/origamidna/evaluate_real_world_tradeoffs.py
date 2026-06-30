#!/usr/bin/env python3
"""Estimate useful capture after practical real-world penalties.

Beginner picture:
Good capture means catching the target EVs. Useful capture means catching the
target EVs without also catching too much junk or damaging the sample.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
IN_FIELD = ROOT / "multitile_capture_field_summary.csv"
OUT_CSV = ROOT / "real_world_tradeoff_scores.csv"
OUT_JSON = ROOT / "real_world_tradeoff_summary.json"
OUT_PLOT = ROOT / "real_world_tradeoff_scores.png"


SCENARIOS = {
    "clean_buffer": {
        "nonspecific_per_encounter": 0.015,
        "damage_per_encounter": 0.005,
        "activity_loss_per_encounter": 0.010,
        "false_capture_weight": 0.35,
        "damage_weight": 0.25,
    },
    "moderate_background": {
        "nonspecific_per_encounter": 0.035,
        "damage_per_encounter": 0.010,
        "activity_loss_per_encounter": 0.020,
        "false_capture_weight": 0.45,
        "damage_weight": 0.30,
    },
    "dirty_biofluid_like": {
        "nonspecific_per_encounter": 0.070,
        "damage_per_encounter": 0.020,
        "activity_loss_per_encounter": 0.040,
        "false_capture_weight": 0.60,
        "damage_weight": 0.40,
    },
}

ENCOUNTERS = range(1, 11)


def cumulative_probability(single_probability: float, encounters: int) -> float:
    return 1.0 - (1.0 - single_probability) ** encounters


def activity_adjusted_target_capture(base_capture: float, encounters: int, activity_loss: float) -> float:
    # Each later pass has slightly less active aptamer surface.
    miss_probability = 1.0
    for i in range(encounters):
        active_fraction = max(0.0, (1.0 - activity_loss) ** i)
        miss_probability *= 1.0 - base_capture * active_fraction
    return 1.0 - miss_probability


def main() -> None:
    with open(IN_FIELD, newline="", encoding="ascii") as f:
        field_rows = list(csv.DictReader(f))

    rows: list[dict[str, str]] = []
    for field in field_rows:
        base_capture = float(field["capture_probability"])
        total_aptamers = int(field["total_aptamers"])
        aptamer_density_factor = min(2.0, total_aptamers / 216.0)
        for scenario_name, scenario in SCENARIOS.items():
            for encounters in ENCOUNTERS:
                target_capture = activity_adjusted_target_capture(
                    base_capture,
                    encounters,
                    scenario["activity_loss_per_encounter"],
                )
                false_capture = cumulative_probability(
                    scenario["nonspecific_per_encounter"] * aptamer_density_factor,
                    encounters,
                )
                damage_loss = cumulative_probability(scenario["damage_per_encounter"], encounters)
                useful_score = (
                    target_capture
                    * (1.0 - scenario["false_capture_weight"] * false_capture)
                    * (1.0 - scenario["damage_weight"] * damage_loss)
                )
                rows.append(
                    {
                        "scenario": scenario_name,
                        "field": field["field"],
                        "unit_layout": field["unit_layout"],
                        "tile_copies": field["tile_copies"],
                        "total_aptamers": field["total_aptamers"],
                        "spacing_nm": field["spacing_nm"],
                        "encounters": str(encounters),
                        "single_encounter_target_capture": f"{base_capture:.4f}",
                        "target_capture_after_encounters": f"{target_capture:.4f}",
                        "false_capture_risk": f"{false_capture:.4f}",
                        "damage_or_loss_risk": f"{damage_loss:.4f}",
                        "useful_capture_score": f"{useful_score:.4f}",
                    }
                )

    rows.sort(key=lambda row: float(row["useful_capture_score"]), reverse=True)
    with open(OUT_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    best_by_scenario = {}
    for scenario_name in SCENARIOS:
        best_by_scenario[scenario_name] = max(
            (row for row in rows if row["scenario"] == scenario_name),
            key=lambda row: float(row["useful_capture_score"]),
        )

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    for scenario_name in SCENARIOS:
        best_field = best_by_scenario[scenario_name]["field"]
        subset = [
            row
            for row in rows
            if row["scenario"] == scenario_name and row["field"] == best_field
        ]
        subset.sort(key=lambda row: int(row["encounters"]))
        ax.plot(
            [int(row["encounters"]) for row in subset],
            [float(row["useful_capture_score"]) for row in subset],
            marker="o",
            label=f"{scenario_name}: {best_field}",
        )
    ax.set_xlabel("near-surface encounters")
    ax.set_ylabel("useful capture score")
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    ax.set_title("Useful capture after background, damage, and activity-loss penalties")
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)

    summary = {
        "model": "real-world useful capture penalty model",
        "input_field_summary": IN_FIELD.name,
        "scenarios": SCENARIOS,
        "best_by_scenario": best_by_scenario,
        "outputs": {
            "tradeoff_scores": OUT_CSV.name,
            "tradeoff_plot": OUT_PLOT.name,
        },
        "interpretation": [
            "Target capture rises with repeated encounters, but false capture, sample damage, and aptamer activity loss also rise.",
            "Clean buffer can tolerate more repeated encounters than dirty biofluid-like conditions.",
            "These are scenario penalties, not measured constants; they should be replaced with experimental controls when available.",
        ],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
    print(f"Wrote {OUT_CSV.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_PLOT.name}")
    for scenario_name, row in best_by_scenario.items():
        print(
            f"{scenario_name}: best={row['field']} encounters={row['encounters']} "
            f"target={row['target_capture_after_encounters']} useful={row['useful_capture_score']} "
            f"false={row['false_capture_risk']}"
        )


if __name__ == "__main__":
    main()
