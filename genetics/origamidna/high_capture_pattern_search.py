#!/usr/bin/env python3
"""Search aggressive aptamer layouts for very high EV capture probability.

Beginner picture:
The earlier designs were like small sticky patches. This script tries bigger
"sticky nets" to see whether a single DNA tile can approach 90% capture in the
current model.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Union

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from clinical_layouts import ellipse_ring, grid_layout
from ev_population_generator import generate_population
from optimize_population_layouts import evaluate_layout
from score_ev_capture_geometry import load_linker_models
from validate_population_dynamics import simulate_one

ROOT = Path(__file__).resolve().parent
POPULATION_FILE = ROOT / "ev_population_optimization_run.npz"
OUT_CANDIDATES_CSV = ROOT / "high_capture_layout_scores.csv"
OUT_DYNAMIC_CSV = ROOT / "high_capture_dynamic_validation.csv"
OUT_LAYOUTS_CSV = ROOT / "high_capture_best_layouts.csv"
OUT_SUMMARY_JSON = ROOT / "high_capture_pattern_search_summary.json"
OUT_PLOT = ROOT / "high_capture_layout_scores.png"
OUT_LAYOUT_PNG = ROOT / "high_capture_best_layouts.png"

RNG_SEED = 20260618
LINKER = "polyT30"
N_POPULATION_EVS = 160
TOP_DYNAMIC_LAYOUTS = 8

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


def multi_patch_layout(
    patch_cols: int,
    patch_rows: int,
    patch_spacing_x: float,
    patch_spacing_y: float,
    local_cols: int,
    local_rows: int,
    local_span_x: float,
    local_span_y: float,
) -> list[Anchor]:
    anchors: list[Anchor] = []
    patch_xs = np.linspace(
        -patch_spacing_x * (patch_cols - 1) / 2.0,
        patch_spacing_x * (patch_cols - 1) / 2.0,
        patch_cols,
    )
    patch_ys = np.linspace(
        -patch_spacing_y * (patch_rows - 1) / 2.0,
        patch_spacing_y * (patch_rows - 1) / 2.0,
        patch_rows,
    )
    local = grid_layout(local_cols, local_rows, local_span_x, local_span_y, LINKER)
    for cx in patch_xs:
        for cy in patch_ys:
            for anchor in local:
                anchors.append(
                    {
                        "x_nm": float(np.clip(cx + float(anchor["x_nm"]), -45.0, 45.0)),
                        "y_nm": float(np.clip(cy + float(anchor["y_nm"]), -30.0, 30.0)),
                        "linker_reach_nm": 15.0,
                        "linker_construct": LINKER,
                    }
                )
    return anchors


def random_high_count_layout(name: str, count: int, rng: np.random.Generator) -> tuple[str, list[Anchor]]:
    # Mix broad coverage and center enrichment.
    n_center = count // 3
    n_broad = count - n_center
    center = rng.normal(0.0, [14.0, 9.0], size=(n_center, 2))
    broad = np.column_stack(
        (
            rng.uniform(-45.0, 45.0, size=n_broad),
            rng.uniform(-30.0, 30.0, size=n_broad),
        )
    )
    xy = np.vstack((center, broad))
    xy[:, 0] = np.clip(xy[:, 0], -45.0, 45.0)
    xy[:, 1] = np.clip(xy[:, 1], -30.0, 30.0)
    return name, [
        {"x_nm": float(x), "y_nm": float(y), "linker_reach_nm": 15.0, "linker_construct": LINKER}
        for x, y in xy
    ]


def candidate_layouts(rng: np.random.Generator) -> dict[str, list[Anchor]]:
    layouts: dict[str, list[Anchor]] = {
        "net_24_6x4": grid_layout(6, 4, 70.0, 38.0, LINKER),
        "net_36_9x4": grid_layout(9, 4, 82.0, 40.0, LINKER),
        "net_48_8x6": grid_layout(8, 6, 84.0, 48.0, LINKER),
        "net_60_10x6": grid_layout(10, 6, 88.0, 52.0, LINKER),
        "net_80_10x8": grid_layout(10, 8, 88.0, 56.0, LINKER),
        "center_net_48": grid_layout(8, 6, 58.0, 38.0, LINKER),
        "center_net_60": grid_layout(10, 6, 62.0, 40.0, LINKER),
        "multi_patch_48": multi_patch_layout(4, 3, 24.0, 18.0, 2, 2, 7.0, 5.0),
        "multi_patch_60": multi_patch_layout(5, 3, 20.0, 18.0, 2, 2, 7.0, 5.0),
        "ring_net_48": ellipse_ring(16, 12.0, 0.75, 0.0, LINKER)
        + ellipse_ring(16, 24.0, 0.75, 0.5, LINKER)
        + ellipse_ring(16, 36.0, 0.70, 0.25, LINKER),
    }
    for count in (36, 48, 60, 80):
        for i in range(16):
            name, anchors = random_high_count_layout(f"random_net_{count}_{i:02d}", count, rng)
            layouts[name] = anchors
    return {name: with_linker(anchors) for name, anchors in layouts.items()}


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


def write_layouts(layouts: dict[str, list[Anchor]], names: list[str]) -> None:
    with open(OUT_LAYOUTS_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["layout", "anchor_id", "x_nm", "y_nm", "linker_construct"],
        )
        writer.writeheader()
        for name in names:
            for i, anchor in enumerate(layouts[name], start=1):
                writer.writerow(
                    {
                        "layout": name,
                        "anchor_id": f"A{i:02d}",
                        "x_nm": f"{float(anchor['x_nm']):.3f}",
                        "y_nm": f"{float(anchor['y_nm']):.3f}",
                        "linker_construct": anchor["linker_construct"],
                    }
                )


def dynamic_validate(
    names: list[str],
    layouts: dict[str, list[Anchor]],
    population,
    linker_models,
) -> list[dict[str, str]]:
    rows = []
    for layout_index, name in enumerate(names):
        rng = np.random.default_rng(RNG_SEED + 500 + layout_index)
        per_ev = []
        for ev_index, diameter in enumerate(population["diameter_nm"]):
            receptor_count = int(population["receptor_count"][ev_index])
            receptor_body = population["receptor_points"][ev_index, :receptor_count, :]
            per_ev.append(simulate_one(layouts[name], linker_models, receptor_body, float(diameter), rng))
        rows.append(
            {
                "layout": name,
                "aptamer_count": str(len(layouts[name])),
                "population_evs": str(len(per_ev)),
                "dynamic_capture_probability": f"{np.mean([r['ever_captured'] for r in per_ev]):.4f}",
                "mean_contacts": f"{np.mean([r['mean_contacts'] for r in per_ev]):.4f}",
                "mean_max_contacts": f"{np.mean([r['mean_max_contacts'] for r in per_ev]):.4f}",
                "strong_fraction_p_ge_2": f"{np.mean([r['strong_fraction'] for r in per_ev]):.4f}",
            }
        )
        print(f"Dynamic validated {name}: {rows[-1]['dynamic_capture_probability']}", flush=True)
    return sorted(rows, key=lambda row: float(row["dynamic_capture_probability"]), reverse=True)


def plot_scores(rows: list[dict[str, str]]) -> None:
    top = rows[:12]
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    labels = [row["layout"] for row in top][::-1]
    values = [float(row["capture_probability"]) for row in top][::-1]
    ax.barh(labels, values, color="#1565c0")
    ax.axvline(0.90, color="crimson", ls="--", lw=1.2, label="90% target")
    ax.set_xlabel("snapshot capture probability")
    ax.set_xlim(0, 1)
    ax.set_title("Aggressive high-capture layout search")
    ax.grid(axis="x", alpha=0.25)
    ax.legend()
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)


def plot_layouts(layouts: dict[str, list[Anchor]], names: list[str]) -> None:
    n = min(6, len(names))
    fig, axes = plt.subplots(2, 3, figsize=(11, 7), constrained_layout=True)
    for ax, name in zip(axes.flat, names[:n]):
        xy = np.asarray([[float(a["x_nm"]), float(a["y_nm"])] for a in layouts[name]])
        ax.add_patch(plt.Rectangle((-45, -30), 90, 60, fill=False, color="0.25", lw=1.2))
        ax.scatter(xy[:, 0], xy[:, 1], s=28, color="#1565c0")
        ax.set_title(f"{name} ({len(layouts[name])})", fontsize=9)
        ax.set_xlim(-50, 50)
        ax.set_ylim(-35, 35)
        ax.set_aspect("equal")
        ax.grid(alpha=0.2)
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.savefig(OUT_LAYOUT_PNG, dpi=220)
    plt.close(fig)


def main() -> None:
    rng = np.random.default_rng(RNG_SEED)
    population = ensure_population()
    linker_models = load_linker_models()
    layouts = candidate_layouts(rng)
    score_rows = []
    for i, (name, anchors) in enumerate(layouts.items()):
        metrics = evaluate_layout(
            anchors,
            linker_models,
            population,
            np.random.default_rng(RNG_SEED + i),
        )
        score_rows.append(
            {
                "layout": name,
                "aptamer_count": str(len(anchors)),
                **{key: f"{value:.5f}" for key, value in metrics.items()},
            }
        )
    score_rows.sort(key=lambda row: float(row["capture_probability"]), reverse=True)
    with open(OUT_CANDIDATES_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(score_rows[0].keys()))
        writer.writeheader()
        writer.writerows(score_rows)
    top_names = [row["layout"] for row in score_rows[:TOP_DYNAMIC_LAYOUTS]]
    write_layouts(layouts, top_names)
    dynamic_rows = dynamic_validate(top_names, layouts, population, linker_models)
    with open(OUT_DYNAMIC_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(dynamic_rows[0].keys()))
        writer.writeheader()
        writer.writerows(dynamic_rows)
    plot_scores(score_rows)
    plot_layouts(layouts, top_names)
    summary = {
        "model": "aggressive high-capture aptamer net search",
        "target_capture_probability": 0.90,
        "rng_seed": RNG_SEED,
        "population_file": POPULATION_FILE.name,
        "layouts_tested": len(layouts),
        "best_snapshot_layout": score_rows[0],
        "best_dynamic_layout": dynamic_rows[0],
        "top_snapshot_rows": score_rows[:10],
        "dynamic_validation_rows": dynamic_rows,
        "outputs": {
            "snapshot_scores": OUT_CANDIDATES_CSV.name,
            "dynamic_validation": OUT_DYNAMIC_CSV.name,
            "best_layout_coordinates": OUT_LAYOUTS_CSV.name,
            "score_plot": OUT_PLOT.name,
            "layout_plot": OUT_LAYOUT_PNG.name,
        },
        "interpretation": [
            "These layouts are intentionally aggressive and may be harder to fabricate because they use many aptamers.",
            "Snapshot capture can approach or exceed 90%, but Brownian dynamic validation is the stricter test.",
            "If dynamic validation remains below 90%, a single tile is probably not enough under sparse clinical EV assumptions.",
        ],
    }
    OUT_SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
    print(f"Wrote {OUT_CANDIDATES_CSV.name}")
    print(f"Wrote {OUT_DYNAMIC_CSV.name}")
    print(f"Wrote {OUT_LAYOUTS_CSV.name}")
    print(f"Wrote {OUT_SUMMARY_JSON.name}")
    print(
        f"Best snapshot: {score_rows[0]['layout']} capture={score_rows[0]['capture_probability']}; "
        f"best dynamic: {dynamic_rows[0]['layout']} capture={dynamic_rows[0]['dynamic_capture_probability']}"
    )


if __name__ == "__main__":
    main()
