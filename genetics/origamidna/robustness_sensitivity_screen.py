#!/usr/bin/env python3
"""Stress-test clinical aptamer layouts across uncertain real-world conditions.

Beginner picture:
One simulation is one weather forecast. This script creates many different
"weather days" by changing EV size, receptor count, aptamer activity, motion,
binding strength, background sticking, and sample loss. A robust layout works
well across many days, not only on one convenient day.
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

from clinical_layouts import clinical_candidate_layouts
from ev_population_generator import receptor_points
from score_ev_capture_geometry import load_linker_models, reach_probability

ROOT = Path(__file__).resolve().parent
OUT_CSV = ROOT / "robustness_sensitivity_results.csv"
OUT_SUMMARY_CSV = ROOT / "robustness_layout_summary.csv"
OUT_JSON = ROOT / "robustness_sensitivity_summary.json"
OUT_PLOT = ROOT / "robustness_layout_comparison.png"
OUT_SENSITIVITY_PLOT = ROOT / "robustness_parameter_sensitivity.png"

RNG_SEED = 20260627
N_SCENARIOS = 80
EVS_PER_SCENARIO = 24
N_STEPS = 240
DT_SECONDS = 0.05
CAPTURE_DWELL_SECONDS = 1.0
SURFACE_CLEARANCE_NM = 2.0
INITIAL_GAP_NM = 8.0
MAX_GAP_NM = 24.0
MAX_ENCOUNTERS = 10

LAYOUT_NAMES = (
    "clinical_grid_12",
    "clinical_grid_18",
    "clinical_dense_24",
    "clinical_hybrid_24",
    "clinical_rescue_24",
)
PATTERNS = ("random", "single_cluster", "two_cluster", "bottom_cap")
PATTERN_PROBABILITIES = (0.40, 0.25, 0.20, 0.15)


def draw_scenario(rng: np.random.Generator) -> dict[str, float]:
    """Draw one plausible set of uncertain physical and practical values."""
    return {
        "mean_diameter_nm": float(rng.uniform(65.0, 82.0)),
        "mean_receptors": float(rng.uniform(3.0, 10.0)),
        "active_aptamer_fraction": float(rng.uniform(0.45, 0.95)),
        "linker_reach_multiplier": float(rng.uniform(0.75, 1.15)),
        "diffusion_nm2_per_s": float(rng.uniform(20.0, 55.0)),
        "k_on_per_step": float(rng.uniform(0.08, 0.24)),
        "k_off_per_s": float(rng.uniform(0.06, 0.30)),
        "nonspecific_per_encounter": float(rng.uniform(0.01, 0.07)),
        "damage_per_encounter": float(rng.uniform(0.003, 0.02)),
        "activity_loss_per_encounter": float(rng.uniform(0.005, 0.04)),
    }


def sample_receptor_count(mean: float, rng: np.random.Generator) -> int:
    # A negative-binomial draw allows both sparse and unusually rich EVs.
    dispersion = 4.0
    probability = dispersion / (dispersion + mean)
    return int(np.clip(rng.negative_binomial(dispersion, probability), 1, 16))


def break_bonds(
    distances: np.ndarray,
    anchor_to_receptor: np.ndarray,
    receptor_to_anchor: np.ndarray,
    k_off_per_s: float,
    rng: np.random.Generator,
) -> None:
    p_off = 1.0 - math.exp(-k_off_per_s * DT_SECONDS)
    for anchor_index, receptor_index in enumerate(anchor_to_receptor):
        if receptor_index == -1:
            continue
        if distances[anchor_index, receptor_index] > 14.0 or rng.random() < p_off:
            anchor_to_receptor[anchor_index] = -1
            receptor_to_anchor[receptor_index] = -1


def try_form_bonds(
    distances: np.ndarray,
    linker_models,
    reach_multiplier: float,
    k_on_per_step: float,
    anchor_to_receptor: np.ndarray,
    receptor_to_anchor: np.ndarray,
    rng: np.random.Generator,
) -> None:
    free_anchors = np.flatnonzero(anchor_to_receptor == -1)
    free_receptors = np.flatnonzero(receptor_to_anchor == -1)
    if len(free_anchors) == 0 or len(free_receptors) == 0:
        return
    sub_distances = distances[np.ix_(free_anchors, free_receptors)] / reach_multiplier
    probabilities = np.asarray(
        [
            [k_on_per_step * reach_probability(linker_models, "polyT30", float(distance)) for distance in row]
            for row in sub_distances
        ]
    )
    candidate_pairs = np.argwhere(rng.random(probabilities.shape) < probabilities)
    rng.shuffle(candidate_pairs)
    for local_anchor, local_receptor in candidate_pairs:
        anchor_index = int(free_anchors[local_anchor])
        receptor_index = int(free_receptors[local_receptor])
        if anchor_to_receptor[anchor_index] == -1 and receptor_to_anchor[receptor_index] == -1:
            anchor_to_receptor[anchor_index] = receptor_index
            receptor_to_anchor[receptor_index] = anchor_index


def simulate_one_ev(
    anchors: np.ndarray,
    receptor_body: np.ndarray,
    diameter_nm: float,
    scenario: dict[str, float],
    linker_models,
    rng: np.random.Generator,
) -> bool:
    radius = diameter_nm / 2.0
    center = np.array(
        [rng.uniform(-22.0, 22.0), rng.uniform(-22.0, 22.0), radius + SURFACE_CLEARANCE_NM + INITIAL_GAP_NM]
    )
    min_z = radius + SURFACE_CLEARANCE_NM
    max_z = min_z + MAX_GAP_NM
    anchor_to_receptor = np.full(len(anchors), -1, dtype=int)
    receptor_to_anchor = np.full(len(receptor_body), -1, dtype=int)
    required_steps = int(round(CAPTURE_DWELL_SECONDS / DT_SECONDS))
    consecutive = 0

    for _ in range(N_STEPS):
        contacts = int(np.count_nonzero(anchor_to_receptor != -1))
        mobility = max(0.10, 1.0 / (1.0 + 0.75 * contacts))
        sigma = math.sqrt(2.0 * scenario["diffusion_nm2_per_s"] * mobility * DT_SECONDS)
        center[:2] += rng.normal(0.0, sigma, size=2)
        center[2] += rng.normal(0.0, sigma * 0.45) - 0.035 * contacts
        center[2] = float(np.clip(center[2], min_z, max_z))
        lateral_norm = float(np.linalg.norm(center[:2]))
        if lateral_norm > 70.0:
            center[:2] *= 70.0 / lateral_norm
        distances = np.linalg.norm(anchors[:, None, :] - (receptor_body + center)[None, :, :], axis=2)
        break_bonds(
            distances,
            anchor_to_receptor,
            receptor_to_anchor,
            scenario["k_off_per_s"],
            rng,
        )
        try_form_bonds(
            distances,
            linker_models,
            scenario["linker_reach_multiplier"],
            scenario["k_on_per_step"],
            anchor_to_receptor,
            receptor_to_anchor,
            rng,
        )
        contacts = int(np.count_nonzero(anchor_to_receptor != -1))
        consecutive = consecutive + 1 if contacts >= 1 else 0
        if consecutive >= required_steps:
            return True
    return False


def repeated_target_capture(single_capture: float, encounters: int, activity_loss: float) -> float:
    miss_probability = 1.0
    for encounter in range(encounters):
        activity = (1.0 - activity_loss) ** encounter
        miss_probability *= 1.0 - single_capture * activity
    return 1.0 - miss_probability


def useful_score(
    single_capture: float,
    encounters: int,
    aptamer_count: int,
    scenario: dict[str, float],
) -> tuple[float, float, float, float]:
    target = repeated_target_capture(single_capture, encounters, scenario["activity_loss_per_encounter"])
    density_factor = aptamer_count / 18.0
    false_risk = 1.0 - (1.0 - min(0.25, scenario["nonspecific_per_encounter"] * density_factor)) ** encounters
    damage_risk = 1.0 - (1.0 - scenario["damage_per_encounter"]) ** encounters
    useful = target * (1.0 - 0.50 * false_risk) * (1.0 - 0.35 * damage_risk)
    return target, false_risk, damage_risk, useful


def rank_correlation(values: np.ndarray, outcomes: np.ndarray) -> float:
    value_ranks = np.argsort(np.argsort(values)).astype(float)
    outcome_ranks = np.argsort(np.argsort(outcomes)).astype(float)
    return float(np.corrcoef(value_ranks, outcome_ranks)[0, 1])


def main() -> None:
    rng = np.random.default_rng(RNG_SEED)
    linker_models = load_linker_models()
    raw_layouts = clinical_candidate_layouts("polyT30")
    rows: list[dict[str, str]] = []

    for scenario_index in range(N_SCENARIOS):
        scenario = draw_scenario(rng)
        evs = []
        for _ in range(EVS_PER_SCENARIO):
            diameter = float(np.clip(rng.normal(scenario["mean_diameter_nm"], 12.0), 45.0, 115.0))
            count = sample_receptor_count(scenario["mean_receptors"], rng)
            pattern = str(rng.choice(PATTERNS, p=PATTERN_PROBABILITIES))
            evs.append((diameter, receptor_points(pattern, count, diameter / 2.0, rng)))

        for layout_name in LAYOUT_NAMES:
            layout = raw_layouts[layout_name]
            active_mask = rng.random(len(layout)) < scenario["active_aptamer_fraction"]
            if not np.any(active_mask):
                active_mask[rng.integers(0, len(layout))] = True
            active_layout = [anchor for anchor, active in zip(layout, active_mask) if active]
            anchors = np.asarray(
                [[float(anchor["x_nm"]), float(anchor["y_nm"]), 0.0] for anchor in active_layout],
                dtype=float,
            )
            captures = [
                simulate_one_ev(anchors, receptors, diameter, scenario, linker_models, rng)
                for diameter, receptors in evs
            ]
            single_capture = float(np.mean(captures))
            encounter_results = []
            for encounters in range(1, MAX_ENCOUNTERS + 1):
                encounter_results.append(
                    (encounters,) + useful_score(single_capture, encounters, len(layout), scenario)
                )
            best = max(encounter_results, key=lambda result: result[-1])
            row = {
                "scenario": str(scenario_index + 1),
                "layout": layout_name,
                "aptamer_count": str(len(layout)),
                "active_aptamers": str(len(active_layout)),
                **{key: f"{value:.5f}" for key, value in scenario.items()},
                "single_encounter_capture": f"{single_capture:.5f}",
                "best_encounters": str(best[0]),
                "target_capture": f"{best[1]:.5f}",
                "false_capture_risk": f"{best[2]:.5f}",
                "damage_risk": f"{best[3]:.5f}",
                "useful_capture_score": f"{best[4]:.5f}",
            }
            rows.append(row)
        print(f"Completed uncertain scenario {scenario_index + 1}/{N_SCENARIOS}", flush=True)

    with open(OUT_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary_rows = []
    for layout_name in LAYOUT_NAMES:
        subset = [row for row in rows if row["layout"] == layout_name]
        useful = np.asarray([float(row["useful_capture_score"]) for row in subset])
        target = np.asarray([float(row["target_capture"]) for row in subset])
        single = np.asarray([float(row["single_encounter_capture"]) for row in subset])
        encounters = np.asarray([int(row["best_encounters"]) for row in subset])
        summary_rows.append(
            {
                "layout": layout_name,
                "aptamer_count": subset[0]["aptamer_count"],
                "median_single_capture": f"{np.median(single):.4f}",
                "median_target_capture": f"{np.median(target):.4f}",
                "median_useful_score": f"{np.median(useful):.4f}",
                "p10_useful_score": f"{np.percentile(useful, 10):.4f}",
                "p90_useful_score": f"{np.percentile(useful, 90):.4f}",
                "probability_useful_ge_0_80": f"{np.mean(useful >= 0.80):.4f}",
                "failure_probability_below_0_50": f"{np.mean(useful < 0.50):.4f}",
                "median_best_encounters": f"{np.median(encounters):.1f}",
            }
        )
    summary_rows.sort(
        key=lambda row: (float(row["p10_useful_score"]), float(row["median_useful_score"])),
        reverse=True,
    )
    with open(OUT_SUMMARY_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    # Correlations show which uncertain inputs most strongly move the winner's score.
    winner = summary_rows[0]["layout"]
    winner_rows = [row for row in rows if row["layout"] == winner]
    outcomes = np.asarray([float(row["useful_capture_score"]) for row in winner_rows])
    parameter_names = list(draw_scenario(np.random.default_rng(0)).keys())
    sensitivity = {
        name: rank_correlation(
            np.asarray([float(row[name]) for row in winner_rows]),
            outcomes,
        )
        for name in parameter_names
    }

    fig, ax = plt.subplots(figsize=(10, 5.4), constrained_layout=True)
    data = [
        [float(row["useful_capture_score"]) for row in rows if row["layout"] == name]
        for name in LAYOUT_NAMES
    ]
    short_names = [name.replace("clinical_", "") for name in LAYOUT_NAMES]
    ax.boxplot(data, labels=short_names, showfliers=False)
    ax.axhline(0.80, color="tab:green", linestyle="--", label="0.80 useful-score target")
    ax.set_ylim(0, 1)
    ax.set_ylabel("useful capture score")
    ax.set_title("Layout reliability across 80 uncertain real-world scenarios")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)

    ordered = sorted(sensitivity.items(), key=lambda item: abs(item[1]), reverse=True)
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    labels = [name.replace("_", " ") for name, _ in ordered]
    values = [value for _, value in ordered]
    colors = ["tab:green" if value > 0 else "tab:red" for value in values]
    ax.barh(labels[::-1], values[::-1], color=colors[::-1])
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlim(-1, 1)
    ax.set_xlabel("rank correlation with useful capture")
    ax.set_title(f"What most changes {winner.replace('clinical_', '')} performance")
    fig.savefig(OUT_SENSITIVITY_PLOT, dpi=220)
    plt.close(fig)

    summary = {
        "model": "Monte Carlo robustness and sensitivity screen",
        "rng_seed": RNG_SEED,
        "scenarios": N_SCENARIOS,
        "evs_per_scenario": EVS_PER_SCENARIO,
        "layout_ranking_rule": "highest 10th-percentile useful score, then highest median",
        "best_robust_layout": summary_rows[0],
        "layout_summary": summary_rows,
        "winner_parameter_correlations": sensitivity,
        "important_limitations": [
            "The uncertainty ranges are plausible modeling assumptions, not measured experimental constants.",
            "Repeated encounters are projected from one Brownian encounter and assume encounters are independent.",
            "Aptamer deactivation is random; fabrication may create spatially patterned failures.",
            "The results compare designs and do not prove a laboratory capture percentage.",
        ],
        "outputs": {
            "all_results": OUT_CSV.name,
            "layout_summary": OUT_SUMMARY_CSV.name,
            "layout_plot": OUT_PLOT.name,
            "sensitivity_plot": OUT_SENSITIVITY_PLOT.name,
        },
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
    print(f"Wrote {OUT_CSV.name}")
    print(f"Wrote {OUT_SUMMARY_CSV.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_PLOT.name}")
    print(f"Wrote {OUT_SENSITIVITY_PLOT.name}")
    print(f"Most robust layout: {winner}")


if __name__ == "__main__":
    main()
