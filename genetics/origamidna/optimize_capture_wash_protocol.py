#!/usr/bin/env python3
"""Optimize capture, antifouling, and washing for useful EV recovery.

The model separates two kinds of attachments:

* target bonds: desired aptamer-CD133 attachments that should survive washing
* background bonds: weaker unwanted attachments that washing should remove

Wash strength is dimensionless because the project does not yet contain flow
rate or measured bond-lifetime data. A value of 1.0 means a moderate reference
wash; laboratory calibration is required before converting it to mL/min.
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

ROOT = Path(__file__).resolve().parent
IN_RESULTS = ROOT / "advanced_capture_design_results.csv"
OUT_ALL = ROOT / "capture_wash_protocol_results.csv"
OUT_BEST = ROOT / "capture_wash_protocol_best.csv"
OUT_JSON = ROOT / "capture_wash_protocol_summary.json"
OUT_PLOT = ROOT / "capture_wash_protocol_comparison.png"
OUT_CURVE = ROOT / "capture_wash_protocol_tradeoff.png"

RNG_SEED = 20260629
TARGET_USEFUL = 0.90
RELIABILITY_PERCENTILE = 10

PASSIVATION_LEVELS = np.arange(0.0, 0.91, 0.10)
WASH_STRENGTHS = np.arange(0.25, 2.01, 0.25)
WASH_DURATIONS_SECONDS = (5, 10, 20, 30, 45, 60, 90, 120)
WASH_CYCLES = (1, 2, 3, 4)
CONTACT_BOOSTS = (0.0, 0.10, 0.20, 0.30)

# Uncertain kinetic ranges used for scenario stress testing. These are model
# assumptions, not experimentally measured constants for this system.
TARGET_KOFF_RANGE_PER_S = (0.00015, 0.00150)
BACKGROUND_KOFF_RANGE_PER_S = (0.010, 0.080)
MAX_FLOW_DAMAGE_PER_CYCLE = 0.004


def p10(values: np.ndarray) -> float:
    return float(np.percentile(values, RELIABILITY_PERCENTILE))


def boost_target_capture(target: np.ndarray, boost: float) -> np.ndarray:
    return 1.0 - np.power(1.0 - target, 1.0 + boost)


def simulate_protocol(
    target: np.ndarray,
    false_risk: np.ndarray,
    existing_damage: np.ndarray,
    target_koff: np.ndarray,
    background_koff: np.ndarray,
    passivation: float,
    strength: float,
    duration_seconds: float,
    cycles: int,
    contact_boost: float,
) -> dict[str, float]:
    total_time = duration_seconds * cycles

    # Stronger flow accelerates bond loss. Unwanted bonds are modeled as more
    # force-sensitive than specific target bonds.
    target_survival = np.exp(-target_koff * total_time * strength**1.15)
    background_survival = np.exp(-background_koff * total_time * strength**1.55)

    target_before_wash = boost_target_capture(target, contact_boost)
    target_after_wash = target_before_wash * target_survival
    false_after_passivation = false_risk * (1.0 - passivation)
    false_after_wash = false_after_passivation * background_survival

    wash_damage = 1.0 - np.power(
        1.0 - MAX_FLOW_DAMAGE_PER_CYCLE * min(strength, 2.0) / 2.0,
        cycles,
    )
    total_damage = 1.0 - (1.0 - existing_damage) * (1.0 - wash_damage)
    useful = target_after_wash * (1.0 - 0.50 * false_after_wash) * (1.0 - 0.35 * total_damage)

    return {
        "p10_useful_capture": p10(useful),
        "median_useful_capture": float(np.median(useful)),
        "p10_target_retained": p10(target_after_wash),
        "median_target_retained": float(np.median(target_after_wash)),
        "p90_background_remaining": float(np.percentile(false_after_wash, 90)),
        "median_background_remaining": float(np.median(false_after_wash)),
        "median_target_wash_survival": float(np.median(target_survival)),
        "median_background_removal": float(1.0 - np.median(background_survival)),
        "reliable_90_fraction": float(np.mean(useful >= TARGET_USEFUL)),
    }


def protocol_cost(
    passivation: float,
    strength: float,
    duration: float,
    cycles: int,
    contact_boost: float,
) -> float:
    """Favor simpler protocols when reliability is otherwise similar."""
    return (
        0.9 * passivation
        + 0.35 * contact_boost
        + 0.04 * strength
        + 0.0005 * duration * cycles
        + 0.01 * max(0, cycles - 1)
    )


def main() -> None:
    with open(IN_RESULTS, newline="", encoding="ascii") as f:
        source = list(csv.DictReader(f))

    designs = sorted({(row["layout"], row["formulation"]) for row in source})
    all_rows: list[dict[str, str]] = []
    best_rows: list[dict[str, str]] = []

    for design_index, (layout, formulation) in enumerate(designs):
        subset = [
            row for row in source
            if row["layout"] == layout and row["formulation"] == formulation
        ]
        target = np.asarray([float(row["target_capture"]) for row in subset])
        false_risk = np.asarray([float(row["false_capture_risk"]) for row in subset])
        damage = np.asarray([float(row["damage_risk"]) for row in subset])

        # Every design uses the same scenario-level kinetic draws so layouts are
        # compared under equivalent good and bad wash conditions.
        scenario_ids = np.asarray([int(row["scenario"]) for row in subset])
        target_koff = np.empty(len(subset))
        background_koff = np.empty(len(subset))
        for i, scenario_id in enumerate(scenario_ids):
            scenario_rng = np.random.default_rng(RNG_SEED + scenario_id)
            target_koff[i] = math.exp(
                scenario_rng.uniform(
                    math.log(TARGET_KOFF_RANGE_PER_S[0]),
                    math.log(TARGET_KOFF_RANGE_PER_S[1]),
                )
            )
            background_koff[i] = math.exp(
                scenario_rng.uniform(
                    math.log(BACKGROUND_KOFF_RANGE_PER_S[0]),
                    math.log(BACKGROUND_KOFF_RANGE_PER_S[1]),
                )
            )

        design_rows = []
        for passivation in PASSIVATION_LEVELS:
            for strength in WASH_STRENGTHS:
                for duration in WASH_DURATIONS_SECONDS:
                    for cycles in WASH_CYCLES:
                        for boost in CONTACT_BOOSTS:
                            metrics = simulate_protocol(
                                target,
                                false_risk,
                                damage,
                                target_koff,
                                background_koff,
                                float(passivation),
                                float(strength),
                                float(duration),
                                cycles,
                                boost,
                            )
                            cost = protocol_cost(
                                float(passivation), float(strength), float(duration), cycles, boost
                            )
                            row = {
                                "layout": layout,
                                "formulation": formulation,
                                "passivation_effectiveness": f"{passivation:.2f}",
                                "wash_strength_relative": f"{strength:.2f}",
                                "wash_duration_seconds": str(duration),
                                "wash_cycles": str(cycles),
                                "contact_boost": f"{boost:.2f}",
                                **{key: f"{value:.4f}" for key, value in metrics.items()},
                                "protocol_cost_index": f"{cost:.4f}",
                                "reliable_90": "yes" if metrics["p10_useful_capture"] >= TARGET_USEFUL else "no",
                            }
                            design_rows.append(row)

        reliable = [row for row in design_rows if row["reliable_90"] == "yes"]
        if reliable:
            best = min(
                reliable,
                key=lambda row: (
                    float(row["protocol_cost_index"]),
                    -float(row["p10_useful_capture"]),
                ),
            )
        else:
            best = max(design_rows, key=lambda row: float(row["p10_useful_capture"]))
        best_rows.append(best)

        # Keep the most informative protocols rather than writing all 327,680 rows.
        design_rows.sort(key=lambda row: float(row["p10_useful_capture"]), reverse=True)
        top_performance = design_rows[:20]
        simplest_reliable = sorted(
            reliable,
            key=lambda row: float(row["protocol_cost_index"]),
        )[:20]
        unique = {
            (
                row["passivation_effectiveness"],
                row["wash_strength_relative"],
                row["wash_duration_seconds"],
                row["wash_cycles"],
                row["contact_boost"],
            ): row
            for row in top_performance + simplest_reliable
        }
        all_rows.extend(unique.values())
        print(
            f"Optimized {design_index + 1}/{len(designs)}: {layout} + {formulation}; "
            f"best p10={best['p10_useful_capture']}",
            flush=True,
        )

    best_rows.sort(
        key=lambda row: (
            row["reliable_90"] == "yes",
            -float(row["protocol_cost_index"]) if row["reliable_90"] == "yes" else float(row["p10_useful_capture"]),
        ),
        reverse=True,
    )

    with open(OUT_ALL, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)
    with open(OUT_BEST, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(best_rows[0].keys()))
        writer.writeheader()
        writer.writerows(best_rows)

    top = best_rows[:10][::-1]
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    labels = [f"{row['layout']} | {row['formulation']}" for row in top]
    values = [float(row["p10_useful_capture"]) for row in top]
    colors = ["#187a72" if row["reliable_90"] == "yes" else "#777777" for row in top]
    ax.barh(labels, values, color=colors)
    ax.axvline(TARGET_USEFUL, color="crimson", linestyle="--", label="reliable 90% target")
    ax.set_xlim(0, 1)
    ax.set_xlabel("10th-percentile useful capture after washing")
    ax.set_title("Optimized capture-and-wash protocols")
    ax.grid(axis="x", alpha=0.25)
    ax.legend()
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)

    winner = best_rows[0]
    winner_candidates = [
        row for row in all_rows
        if row["layout"] == winner["layout"]
        and row["formulation"] == winner["formulation"]
        and row["contact_boost"] == winner["contact_boost"]
        and row["passivation_effectiveness"] == winner["passivation_effectiveness"]
    ]
    winner_candidates.sort(key=lambda row: float(row["median_background_removal"]))
    fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    ax.scatter(
        [float(row["median_background_removal"]) for row in winner_candidates],
        [float(row["median_target_wash_survival"]) for row in winner_candidates],
        c=[float(row["p10_useful_capture"]) for row in winner_candidates],
        cmap="viridis",
        vmin=0.75,
        vmax=1.0,
        s=55,
    )
    ax.set_xlabel("median background removed by wash")
    ax.set_ylabel("median target bonds surviving wash")
    ax.set_title("Wash trade-off for the leading design")
    ax.grid(alpha=0.25)
    fig.savefig(OUT_CURVE, dpi=220)
    plt.close(fig)

    reliable_designs = [row for row in best_rows if row["reliable_90"] == "yes"]
    summary = {
        "model": "kinetic capture-and-wash protocol optimizer",
        "reliable_90_definition": "10th-percentile useful capture is at least 0.90",
        "wash_strength_definition": "dimensionless relative strength; 1.0 is a moderate reference wash",
        "kinetic_assumptions": {
            "target_koff_range_per_s": TARGET_KOFF_RANGE_PER_S,
            "background_koff_range_per_s": BACKGROUND_KOFF_RANGE_PER_S,
            "maximum_flow_damage_per_cycle": MAX_FLOW_DAMAGE_PER_CYCLE,
        },
        "protocols_evaluated_per_design": int(
            len(PASSIVATION_LEVELS)
            * len(WASH_STRENGTHS)
            * len(WASH_DURATIONS_SECONDS)
            * len(WASH_CYCLES)
            * len(CONTACT_BOOSTS)
        ),
        "designs_evaluated": len(designs),
        "best_protocol": winner,
        "reliable_design_count": len(reliable_designs),
        "limitations": [
            "Wash strength cannot be converted to a laboratory flow rate until bond lifetimes are measured.",
            "Target and background off-rate ranges are explicit assumptions, not measured values.",
            "The model assumes bound populations decay exponentially during each wash.",
            "A simulated 90% result is a protocol hypothesis, not experimental validation.",
        ],
        "outputs": {
            "selected_protocols": OUT_ALL.name,
            "best_by_design": OUT_BEST.name,
            "comparison_plot": OUT_PLOT.name,
            "tradeoff_plot": OUT_CURVE.name,
        },
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
    print(f"Wrote {OUT_ALL.name}")
    print(f"Wrote {OUT_BEST.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(
        f"Best protocol: {winner['layout']} + {winner['formulation']}; "
        f"p10 useful={winner['p10_useful_capture']} reliable={winner['reliable_90']}"
    )


if __name__ == "__main__":
    main()
