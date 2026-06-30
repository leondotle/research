#!/usr/bin/env python3
"""Search layouts, spacer/surface formulations, and staged capture zones.

This is a comparative design model, not a laboratory claim. The formulation
effects below are explicit assumptions that must later be replaced by measured
values. "Reliable 90%" means at least 90% useful capture in the 10th-percentile
stress-test scenario, not merely a 90% average target-contact probability.
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

from clinical_layouts import clinical_candidate_layouts, ellipse_ring, grid_layout
from ev_population_generator import receptor_points
from robustness_sensitivity_screen import (
    PATTERNS,
    PATTERN_PROBABILITIES,
    draw_scenario,
    sample_receptor_count,
    simulate_one_ev,
)
from score_ev_capture_geometry import load_linker_models

ROOT = Path(__file__).resolve().parent
OUT_RESULTS = ROOT / "advanced_capture_design_results.csv"
OUT_SUMMARY = ROOT / "advanced_capture_design_summary.csv"
OUT_LAYOUTS = ROOT / "advanced_capture_top_layouts.csv"
OUT_JSON = ROOT / "advanced_capture_design_summary.json"
OUT_PLOT = ROOT / "advanced_capture_reliability.png"
OUT_LAYOUT_PLOT = ROOT / "advanced_capture_top_layouts.png"

RNG_SEED = 20260628
N_SCENARIOS = 36
EVS_PER_SCENARIO = 18
MAX_ZONES = 16
TARGET = 0.90

# These are hypotheses for screening, not measured material constants.
FORMULATIONS = {
    "polyT30_standard": {
        "reach_multiplier": 1.00,
        "active_fraction_bonus": 0.00,
        "kon_multiplier": 1.00,
        "koff_multiplier": 1.00,
        "nonspecific_multiplier": 1.00,
        "activity_loss_multiplier": 1.00,
        "damage_multiplier": 1.00,
        "description": "Current polyT30 spacer and unmodified modeled surface.",
    },
    "polyT30_PEG_passivated": {
        "reach_multiplier": 1.00,
        "active_fraction_bonus": 0.08,
        "kon_multiplier": 1.00,
        "koff_multiplier": 1.00,
        "nonspecific_multiplier": 0.45,
        "activity_loss_multiplier": 0.70,
        "damage_multiplier": 0.80,
        "description": "polyT30 with an assumed PEG-like antifouling surface.",
    },
    "rigid_spacer_PEG": {
        "reach_multiplier": 1.25,
        "active_fraction_bonus": 0.08,
        "kon_multiplier": 0.95,
        "koff_multiplier": 0.90,
        "nonspecific_multiplier": 0.50,
        "activity_loss_multiplier": 0.65,
        "damage_multiplier": 0.80,
        "description": "Longer, more upright modeled spacer on a PEG-like surface.",
    },
    "mixed_spacer_antifouling": {
        "reach_multiplier": 1.15,
        "active_fraction_bonus": 0.12,
        "kon_multiplier": 1.05,
        "koff_multiplier": 0.85,
        "nonspecific_multiplier": 0.30,
        "activity_loss_multiplier": 0.50,
        "damage_multiplier": 0.65,
        "description": "Optimistic mixed-length spacer plus strong antifouling layer.",
    },
}


def sunflower_layout(count: int, radius_nm: float) -> list[dict[str, float | str]]:
    anchors = []
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(count):
        radius = radius_nm * math.sqrt((i + 0.5) / count)
        theta = i * golden_angle
        anchors.append(
            {
                "x_nm": radius * math.cos(theta),
                "y_nm": 0.72 * radius * math.sin(theta),
                "linker_reach_nm": 15.0,
                "linker_construct": "polyT30",
            }
        )
    return anchors


def multi_patch(count: int, radius_nm: float) -> list[dict[str, float | str]]:
    centers = [(-radius_nm, 0.0), (radius_nm, 0.0), (0.0, -0.65 * radius_nm), (0.0, 0.65 * radius_nm)]
    per_patch = count // len(centers)
    extras = count - per_patch * len(centers)
    anchors = []
    for patch_index, (cx, cy) in enumerate(centers):
        local_count = per_patch + (1 if patch_index < extras else 0)
        for i in range(local_count):
            theta = 2.0 * math.pi * i / local_count
            anchors.append(
                {
                    "x_nm": cx + 5.0 * math.cos(theta),
                    "y_nm": cy + 4.0 * math.sin(theta),
                    "linker_reach_nm": 15.0,
                    "linker_construct": "polyT30",
                }
            )
    return anchors


def candidate_layouts() -> dict[str, list[dict[str, float | str]]]:
    clinical = clinical_candidate_layouts("polyT30")
    return {
        "rescue_24": clinical["clinical_rescue_24"],
        "hybrid_24": clinical["clinical_hybrid_24"],
        "broad_grid_20": grid_layout(5, 4, 48.0, 28.0, "polyT30"),
        "broad_grid_24": grid_layout(6, 4, 54.0, 30.0, "polyT30"),
        "triple_ring_24": ellipse_ring(8, 10.0, 0.72, 0.0, "polyT30")
        + ellipse_ring(8, 20.0, 0.72, 0.5, "polyT30")
        + ellipse_ring(8, 30.0, 0.72, 0.25, "polyT30"),
        "sunflower_20": sunflower_layout(20, 29.0),
        "sunflower_24": sunflower_layout(24, 32.0),
        "four_patch_24": multi_patch(24, 19.0),
    }


def apply_formulation(base: dict[str, float], formulation: dict[str, float | str]) -> dict[str, float]:
    scenario = dict(base)
    scenario["linker_reach_multiplier"] *= float(formulation["reach_multiplier"])
    scenario["active_aptamer_fraction"] = min(
        0.99,
        scenario["active_aptamer_fraction"] + float(formulation["active_fraction_bonus"]),
    )
    scenario["k_on_per_step"] *= float(formulation["kon_multiplier"])
    scenario["k_off_per_s"] *= float(formulation["koff_multiplier"])
    scenario["nonspecific_per_encounter"] *= float(formulation["nonspecific_multiplier"])
    scenario["activity_loss_per_encounter"] *= float(formulation["activity_loss_multiplier"])
    scenario["damage_per_encounter"] *= float(formulation["damage_multiplier"])
    return scenario


def staged_score(
    single_capture: float,
    zones: int,
    aptamer_count: int,
    scenario: dict[str, float],
) -> tuple[float, float, float, float]:
    miss = 1.0
    for zone in range(zones):
        activity = (1.0 - scenario["activity_loss_per_encounter"]) ** zone
        miss *= 1.0 - single_capture * activity
    target_capture = 1.0 - miss
    density_factor = aptamer_count / 18.0
    false_per_zone = min(0.25, scenario["nonspecific_per_encounter"] * density_factor)
    false_risk = 1.0 - (1.0 - false_per_zone) ** zones
    # A staged flow surface gives several contact opportunities during one pass.
    # Damage grows sublinearly because this is not repeated collection/transfer.
    effective_damage_exposures = 1.0 + 0.25 * max(0, zones - 1)
    damage_risk = 1.0 - (1.0 - scenario["damage_per_encounter"]) ** effective_damage_exposures
    useful = target_capture * (1.0 - 0.50 * false_risk) * (1.0 - 0.35 * damage_risk)
    return target_capture, false_risk, damage_risk, useful


def main() -> None:
    rng = np.random.default_rng(RNG_SEED)
    linker_models = load_linker_models()
    layouts = candidate_layouts()
    all_results: list[dict[str, str]] = []

    scenarios = []
    for _ in range(N_SCENARIOS):
        base = draw_scenario(rng)
        evs = []
        for _ in range(EVS_PER_SCENARIO):
            diameter = float(np.clip(rng.normal(base["mean_diameter_nm"], 12.0), 45.0, 115.0))
            count = sample_receptor_count(base["mean_receptors"], rng)
            pattern = str(rng.choice(PATTERNS, p=PATTERN_PROBABILITIES))
            evs.append((diameter, receptor_points(pattern, count, diameter / 2.0, rng)))
        scenarios.append((base, evs))

    design_index = 0
    for formulation_name, formulation in FORMULATIONS.items():
        for layout_name, layout in layouts.items():
            design_index += 1
            for scenario_index, (base, evs) in enumerate(scenarios):
                scenario = apply_formulation(base, formulation)
                active_mask = rng.random(len(layout)) < scenario["active_aptamer_fraction"]
                if not np.any(active_mask):
                    active_mask[rng.integers(0, len(layout))] = True
                active_layout = [anchor for anchor, active in zip(layout, active_mask) if active]
                anchors = np.asarray(
                    [[float(anchor["x_nm"]), float(anchor["y_nm"]), 0.0] for anchor in active_layout]
                )
                captures = [
                    simulate_one_ev(anchors, receptors, diameter, scenario, linker_models, rng)
                    for diameter, receptors in evs
                ]
                single_capture = float(np.mean(captures))
                candidates = [
                    (zones,) + staged_score(single_capture, zones, len(layout), scenario)
                    for zones in range(1, MAX_ZONES + 1)
                ]
                best = max(candidates, key=lambda result: result[-1])
                all_results.append(
                    {
                        "scenario": str(scenario_index + 1),
                        "layout": layout_name,
                        "formulation": formulation_name,
                        "aptamer_count": str(len(layout)),
                        "active_aptamers": str(len(active_layout)),
                        "single_capture": f"{single_capture:.5f}",
                        "best_zones": str(best[0]),
                        "target_capture": f"{best[1]:.5f}",
                        "false_capture_risk": f"{best[2]:.5f}",
                        "damage_risk": f"{best[3]:.5f}",
                        "useful_capture": f"{best[4]:.5f}",
                    }
                )
            print(f"Completed design {design_index}/{len(layouts) * len(FORMULATIONS)}: {layout_name} + {formulation_name}", flush=True)

    with open(OUT_RESULTS, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_results[0].keys()))
        writer.writeheader()
        writer.writerows(all_results)

    summary_rows = []
    for formulation_name in FORMULATIONS:
        for layout_name, layout in layouts.items():
            subset = [
                row for row in all_results
                if row["layout"] == layout_name and row["formulation"] == formulation_name
            ]
            useful = np.asarray([float(row["useful_capture"]) for row in subset])
            target = np.asarray([float(row["target_capture"]) for row in subset])
            single = np.asarray([float(row["single_capture"]) for row in subset])
            zones = np.asarray([int(row["best_zones"]) for row in subset])
            summary_rows.append(
                {
                    "layout": layout_name,
                    "formulation": formulation_name,
                    "aptamer_count": str(len(layout)),
                    "median_single_capture": f"{np.median(single):.4f}",
                    "median_target_capture": f"{np.median(target):.4f}",
                    "p10_target_capture": f"{np.percentile(target, 10):.4f}",
                    "median_useful_capture": f"{np.median(useful):.4f}",
                    "p10_useful_capture": f"{np.percentile(useful, 10):.4f}",
                    "fraction_useful_ge_90": f"{np.mean(useful >= TARGET):.4f}",
                    "fraction_target_ge_90": f"{np.mean(target >= TARGET):.4f}",
                    "median_best_zones": f"{np.median(zones):.1f}",
                    "reliable_90_useful": "yes" if np.percentile(useful, 10) >= TARGET else "no",
                }
            )
    summary_rows.sort(key=lambda row: float(row["p10_useful_capture"]), reverse=True)
    with open(OUT_SUMMARY, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    top_layout_names = []
    for row in summary_rows:
        if row["layout"] not in top_layout_names:
            top_layout_names.append(row["layout"])
        if len(top_layout_names) == 6:
            break
    with open(OUT_LAYOUTS, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=["layout", "anchor_id", "x_nm", "y_nm"])
        writer.writeheader()
        for name in top_layout_names:
            for i, anchor in enumerate(layouts[name], 1):
                writer.writerow(
                    {"layout": name, "anchor_id": i, "x_nm": f"{float(anchor['x_nm']):.3f}", "y_nm": f"{float(anchor['y_nm']):.3f}"}
                )

    top = summary_rows[:12][::-1]
    fig, ax = plt.subplots(figsize=(11, 6.5), constrained_layout=True)
    labels = [f"{r['layout']} | {r['formulation']}" for r in top]
    values = [float(r["p10_useful_capture"]) for r in top]
    ax.barh(labels, values, color="#187a72")
    ax.axvline(TARGET, color="crimson", linestyle="--", label="reliable 90% target")
    ax.set_xlim(0, 1)
    ax.set_xlabel("10th-percentile useful capture")
    ax.set_title("Best designs under difficult simulated conditions")
    ax.grid(axis="x", alpha=0.25)
    ax.legend()
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(11, 7), constrained_layout=True)
    for ax, name in zip(axes.flat, top_layout_names):
        xy = np.asarray([[float(a["x_nm"]), float(a["y_nm"])] for a in layouts[name]])
        ax.add_patch(plt.Rectangle((-45, -30), 90, 60, fill=False, color="0.3"))
        ax.scatter(xy[:, 0], xy[:, 1], s=28, color="#187a72")
        ax.set_title(f"{name} ({len(layouts[name])} aptamers)")
        ax.set_xlim(-50, 50)
        ax.set_ylim(-35, 35)
        ax.set_aspect("equal")
        ax.grid(alpha=0.2)
    fig.savefig(OUT_LAYOUT_PLOT, dpi=220)
    plt.close(fig)

    reliable = [row for row in summary_rows if row["reliable_90_useful"] == "yes"]
    report = {
        "model": "joint layout, spacer/surface formulation, and staged-zone search",
        "reliable_90_definition": "10th-percentile useful capture is at least 0.90",
        "scenarios": N_SCENARIOS,
        "evs_per_scenario": EVS_PER_SCENARIO,
        "designs_tested": len(summary_rows),
        "formulation_assumptions": FORMULATIONS,
        "best_design": summary_rows[0],
        "reliable_90_design_count": len(reliable),
        "reliable_90_designs": reliable,
        "limitations": [
            "Formulation multipliers are hypotheses, not experimental measurements.",
            "The same A15 recognition sequence is modeled throughout; this does not compare validated aptamer sequences.",
            "Capture zones are treated as approximately independent opportunities during one surface pass.",
            "A 90% simulated score is a design target, not proof of 90% laboratory or clinical performance.",
        ],
        "outputs": {
            "all_results": OUT_RESULTS.name,
            "ranked_summary": OUT_SUMMARY.name,
            "top_layout_coordinates": OUT_LAYOUTS.name,
            "reliability_plot": OUT_PLOT.name,
            "layout_plot": OUT_LAYOUT_PLOT.name,
        },
    }
    OUT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="ascii")
    print(f"Wrote {OUT_RESULTS.name}")
    print(f"Wrote {OUT_SUMMARY.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Best: {summary_rows[0]['layout']} + {summary_rows[0]['formulation']} p10 useful={summary_rows[0]['p10_useful_capture']}")
    print(f"Reliable 90% designs: {len(reliable)}")


if __name__ == "__main__":
    main()
