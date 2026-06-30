#!/usr/bin/env python3
"""Calculate what must improve to reach reliable 90% useful EV capture.

Beginner picture:
The forward model asks, "What score does this design get?" This inverse model
asks, "How much cleaner or stronger must this design become to get 90%?"

The calculations are design requirements, not predictions that a named
material will automatically achieve them.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
IN_RESULTS = ROOT / "advanced_capture_design_results.csv"
OUT_CSV = ROOT / "inverse_design_requirements.csv"
OUT_JSON = ROOT / "inverse_design_requirements_summary.json"
OUT_PLOT = ROOT / "inverse_design_requirements.png"
OUT_FRONTIER = ROOT / "inverse_design_tradeoff_frontier.csv"

TARGET = 0.90
RELIABILITY_PERCENTILE = 10
BACKGROUND_REDUCTION_GRID = np.linspace(0.0, 0.99, 100)
CONTACT_BOOST_GRID = np.linspace(0.0, 1.50, 61)
WASH_REMOVAL_GRID = np.linspace(0.0, 0.99, 100)
MAX_WASH_TARGET_LOSS = 0.02


def useful_capture(target: np.ndarray, false_risk: np.ndarray, damage: np.ndarray) -> np.ndarray:
    return target * (1.0 - 0.50 * false_risk) * (1.0 - 0.35 * damage)


def boosted_target(target: np.ndarray, boost: float) -> np.ndarray:
    """Add a fractional independent binding opportunity.

    boost=0 means no change. boost=1 approximates one equally strong,
    independent second binding route: the original miss probability is squared.
    """
    return 1.0 - np.power(1.0 - target, 1.0 + boost)


def p10(values: np.ndarray) -> float:
    return float(np.percentile(values, RELIABILITY_PERCENTILE))


def minimum_background_reduction(target: np.ndarray, false: np.ndarray, damage: np.ndarray) -> float | None:
    for reduction in BACKGROUND_REDUCTION_GRID:
        score = useful_capture(target, false * (1.0 - reduction), damage)
        if p10(score) >= TARGET:
            return float(reduction)
    return None


def minimum_contact_boost(target: np.ndarray, false: np.ndarray, damage: np.ndarray) -> float | None:
    for boost in CONTACT_BOOST_GRID:
        score = useful_capture(boosted_target(target, float(boost)), false, damage)
        if p10(score) >= TARGET:
            return float(boost)
    return None


def minimum_wash_removal(target: np.ndarray, false: np.ndarray, damage: np.ndarray) -> float | None:
    for removal in WASH_REMOVAL_GRID:
        # Stronger washing is assumed to lose up to 2% of specifically bound targets.
        retained_target = target * (1.0 - MAX_WASH_TARGET_LOSS * removal)
        score = useful_capture(retained_target, false * (1.0 - removal), damage)
        if p10(score) >= TARGET:
            return float(removal)
    return None


def combined_requirement(
    target: np.ndarray,
    false: np.ndarray,
    damage: np.ndarray,
) -> tuple[float, float, float, list[dict[str, float]]] | None:
    feasible = []
    for boost in CONTACT_BOOST_GRID:
        target_after = boosted_target(target, float(boost))
        for reduction in BACKGROUND_REDUCTION_GRID:
            score = useful_capture(target_after, false * (1.0 - reduction), damage)
            robust_score = p10(score)
            if robust_score >= TARGET:
                # Equal weighting gives an understandable "smallest total change" solution.
                burden = float(boost) + float(reduction)
                feasible.append(
                    {
                        "contact_boost": float(boost),
                        "background_reduction": float(reduction),
                        "p10_useful_capture": robust_score,
                        "burden": burden,
                    }
                )
                break
    if not feasible:
        return None
    feasible.sort(key=lambda row: (row["burden"], row["contact_boost"], row["background_reduction"]))
    best = feasible[0]
    return (
        best["contact_boost"],
        best["background_reduction"],
        best["p10_useful_capture"],
        feasible,
    )


def fmt_optional(value: float | None) -> str:
    return "not_achievable_in_range" if value is None else f"{value:.4f}"


def main() -> None:
    with open(IN_RESULTS, newline="", encoding="ascii") as f:
        source_rows = list(csv.DictReader(f))

    designs = sorted({(row["layout"], row["formulation"]) for row in source_rows})
    rows = []
    frontiers = []
    for layout, formulation in designs:
        subset = [
            row for row in source_rows
            if row["layout"] == layout and row["formulation"] == formulation
        ]
        target = np.asarray([float(row["target_capture"]) for row in subset])
        false = np.asarray([float(row["false_capture_risk"]) for row in subset])
        damage = np.asarray([float(row["damage_risk"]) for row in subset])
        baseline = useful_capture(target, false, damage)
        background_only = minimum_background_reduction(target, false, damage)
        contact_only = minimum_contact_boost(target, false, damage)
        wash_only = minimum_wash_removal(target, false, damage)
        combined = combined_requirement(target, false, damage)

        if combined is None:
            combined_boost = None
            combined_background = None
            combined_score = None
        else:
            combined_boost, combined_background, combined_score, feasible = combined
            for point in feasible:
                frontiers.append(
                    {
                        "layout": layout,
                        "formulation": formulation,
                        "contact_boost": f"{point['contact_boost']:.4f}",
                        "background_reduction": f"{point['background_reduction']:.4f}",
                        "p10_useful_capture": f"{point['p10_useful_capture']:.4f}",
                        "burden": f"{point['burden']:.4f}",
                    }
                )

        rows.append(
            {
                "layout": layout,
                "formulation": formulation,
                "baseline_p10_target": f"{p10(target):.4f}",
                "baseline_p10_false_risk": f"{p10(false):.4f}",
                "baseline_p10_useful": f"{p10(baseline):.4f}",
                "background_reduction_needed_alone": fmt_optional(background_only),
                "wash_removal_needed_alone": fmt_optional(wash_only),
                "contact_boost_needed_alone": fmt_optional(contact_only),
                "combined_contact_boost": fmt_optional(combined_boost),
                "combined_background_reduction": fmt_optional(combined_background),
                "combined_p10_useful": fmt_optional(combined_score),
                "combined_total_change": fmt_optional(
                    None if combined_boost is None or combined_background is None else combined_boost + combined_background
                ),
            }
        )

    def numeric_or_large(value: str) -> float:
        return float(value) if value != "not_achievable_in_range" else 999.0

    rows.sort(key=lambda row: numeric_or_large(row["combined_total_change"]))
    with open(OUT_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    frontiers.sort(key=lambda row: float(row["burden"]))
    with open(OUT_FRONTIER, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(frontiers[0].keys()))
        writer.writeheader()
        writer.writerows(frontiers)

    top = rows[:10][::-1]
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    labels = [f"{row['layout']} | {row['formulation']}" for row in top]
    backgrounds = [numeric_or_large(row["combined_background_reduction"]) for row in top]
    boosts = [numeric_or_large(row["combined_contact_boost"]) for row in top]
    y = np.arange(len(top))
    ax.barh(y, backgrounds, label="background reduction", color="#187a72")
    ax.barh(y, boosts, left=backgrounds, label="extra contact strength", color="#df8f44")
    ax.set_yticks(y, labels)
    ax.set_xlim(0, max(1.0, max(a + b for a, b in zip(backgrounds, boosts)) * 1.08))
    ax.set_xlabel("minimum combined modeled improvement")
    ax.set_title("Smallest changes calculated to reach reliable 90% useful capture")
    ax.grid(axis="x", alpha=0.25)
    ax.legend()
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)

    best = rows[0]
    summary = {
        "model": "inverse requirements for reliable useful capture",
        "target_useful_capture": TARGET,
        "reliability_definition": "10th-percentile useful capture across stress-test scenarios",
        "best_path": best,
        "interpretation": {
            "background_reduction": "Fraction of current nonspecific binding that must be prevented or removed.",
            "contact_boost": "Fractional extra independent binding opportunity; 1.0 approximates a second equally strong independent route.",
            "wash_removal": "Fraction of nonspecific material removed, assuming no more than 2% target loss at the strongest wash.",
        },
        "important_limitations": [
            "The calculator works backward from simulated outcomes, not laboratory measurements.",
            "Contact boost is a generic requirement and does not prove that a specific second aptamer will provide it.",
            "Background reduction may come from passivation, washing, blocking buffer, or a combination.",
            "The model does not include every possible interaction in patient biofluid.",
        ],
        "outputs": {
            "requirements": OUT_CSV.name,
            "tradeoff_frontier": OUT_FRONTIER.name,
            "requirements_plot": OUT_PLOT.name,
        },
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
    print(f"Wrote {OUT_CSV.name}")
    print(f"Wrote {OUT_FRONTIER.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_PLOT.name}")
    print(
        f"Best path: {best['layout']} + {best['formulation']}; "
        f"contact boost={best['combined_contact_boost']}, "
        f"background reduction={best['combined_background_reduction']}"
    )


if __name__ == "__main__":
    main()
