#!/usr/bin/env python3
"""Optimize aptamer layouts against a realistic sparse CD133+ EV population.

Beginner picture:
We make one shared test set of EVs. Then every layout tries to catch those same
EVs. The best layout is the one that catches the most EVs across the whole mix,
not just one perfect vesicle.
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

from clinical_layouts import clinical_candidate_layouts, grid_layout
from ev_population_generator import generate_population
from score_ev_capture_geometry import load_layouts, load_linker_models, reach_probability

ROOT = Path(__file__).resolve().parent
OUT_LAYOUTS_CSV = ROOT / "population_optimized_layouts.csv"
OUT_SCORES_CSV = ROOT / "population_layout_scores.csv"
OUT_SUMMARY_JSON = ROOT / "population_layout_optimization_summary.json"
OUT_PLOT = ROOT / "population_layout_top_scores.png"
OUT_LAYOUT_PNG = ROOT / "population_optimized_layouts.png"

RNG_SEED = 20260616
LINKER = "polyT30"
N_POPULATION_EVS = 160
N_RANDOM_LAYOUTS = 72
N_GENERATIONS = 5
KEEP_PER_GENERATION = 10
MUTANTS_PER_PARENT = 5
MAX_APTAMERS = 24
MIN_APTAMERS = 10
N_BINDING_TRIALS = 10
LATERAL_OFFSETS_NM = (
    (0.0, 0.0),
    (-14.0, 0.0),
    (14.0, 0.0),
    (0.0, -10.0),
    (0.0, 10.0),
)
SURFACE_CLEARANCE_NM = 2.0

Anchor = dict[str, Union[float, str]]


def with_linker(anchors: list[Anchor], construct: str = LINKER) -> list[Anchor]:
    return [
        {
            **anchor,
            "linker_construct": construct,
            "linker_reach_nm": float(anchor.get("linker_reach_nm", 15.0)),
        }
        for anchor in anchors
    ]


def anchor_array(anchors: list[Anchor]) -> np.ndarray:
    return np.asarray(
        [[float(anchor["x_nm"]), float(anchor["y_nm"]), 0.0] for anchor in anchors],
        dtype=float,
    )


def random_layout(name: str, aptamer_count: int, rng: np.random.Generator) -> tuple[str, list[Anchor]]:
    # Half broad coverage, half central-biased. This gives the search both
    # "catch wide" and "hold near center" ingredients.
    n_center = aptamer_count // 2
    n_broad = aptamer_count - n_center
    center = rng.normal(0.0, [12.0, 8.0], size=(n_center, 2))
    broad = np.column_stack(
        (
            rng.uniform(-42.0, 42.0, size=n_broad),
            rng.uniform(-28.0, 28.0, size=n_broad),
        )
    )
    xy = np.vstack((center, broad))
    xy[:, 0] = np.clip(xy[:, 0], -45.0, 45.0)
    xy[:, 1] = np.clip(xy[:, 1], -30.0, 30.0)
    anchors = [
        {
            "x_nm": float(x),
            "y_nm": float(y),
            "linker_reach_nm": 15.0,
            "linker_construct": LINKER,
        }
        for x, y in xy
    ]
    return name, anchors


def mutate_layout(name: str, parent: list[Anchor], rng: np.random.Generator) -> tuple[str, list[Anchor]]:
    xy = np.asarray([[float(a["x_nm"]), float(a["y_nm"])] for a in parent], dtype=float)
    xy += rng.normal(0.0, [5.0, 3.5], size=xy.shape)
    if rng.random() < 0.30 and len(xy) > MIN_APTAMERS:
        xy = np.delete(xy, int(rng.integers(0, len(xy))), axis=0)
    if rng.random() < 0.45 and len(xy) < MAX_APTAMERS:
        extra = rng.normal(0.0, [18.0, 12.0], size=(1, 2))
        xy = np.vstack((xy, extra))
    xy[:, 0] = np.clip(xy[:, 0], -45.0, 45.0)
    xy[:, 1] = np.clip(xy[:, 1], -30.0, 30.0)
    anchors = [
        {
            "x_nm": float(x),
            "y_nm": float(y),
            "linker_reach_nm": 15.0,
            "linker_construct": LINKER,
        }
        for x, y in xy
    ]
    return name, anchors


def seed_layouts(rng: np.random.Generator) -> dict[str, list[Anchor]]:
    layouts = {name: with_linker(anchors) for name, anchors in load_layouts().items()}
    layouts.update(clinical_candidate_layouts(LINKER))
    layouts["population_grid_24_wide"] = grid_layout(6, 4, 66.0, 34.0, LINKER)
    layouts["population_grid_18_wide"] = grid_layout(6, 3, 66.0, 30.0, LINKER)
    layouts["population_grid_15_mid"] = grid_layout(5, 3, 52.0, 24.0, LINKER)
    for i in range(N_RANDOM_LAYOUTS):
        count = int(rng.integers(MIN_APTAMERS, MAX_APTAMERS + 1))
        name, anchors = random_layout(f"evo_seed_{i:03d}", count, rng)
        layouts[name] = anchors
    return layouts


def max_bipartite_matches(edges: np.ndarray) -> int:
    if edges.size == 0:
        return 0
    n_anchors, n_receptors = edges.shape
    match_to_anchor = np.full(n_receptors, -1, dtype=int)
    adjacency = [np.flatnonzero(edges[i]) for i in range(n_anchors)]

    def assign(anchor_index: int, seen: np.ndarray) -> bool:
        for receptor_index in adjacency[anchor_index]:
            if seen[receptor_index]:
                continue
            seen[receptor_index] = True
            if match_to_anchor[receptor_index] == -1 or assign(match_to_anchor[receptor_index], seen):
                match_to_anchor[receptor_index] = anchor_index
                return True
        return False

    matches = 0
    for anchor_index in sorted(range(n_anchors), key=lambda i: len(adjacency[i])):
        if assign(anchor_index, np.zeros(n_receptors, dtype=bool)):
            matches += 1
    return matches


def layout_penalty(anchors: list[Anchor]) -> float:
    # A gentle crowding penalty. If many hooks sit almost on top of each other,
    # the score goes down slightly because real DNA/linkers occupy space.
    xy = np.asarray([[float(a["x_nm"]), float(a["y_nm"])] for a in anchors], dtype=float)
    if len(xy) < 2:
        return 0.0
    close_pairs = 0
    for i in range(len(xy)):
        d = np.linalg.norm(xy[i + 1 :] - xy[i], axis=1)
        close_pairs += int(np.sum(d < 5.0))
    return min(0.15, 0.004 * close_pairs + 0.002 * max(0, len(xy) - 18))


def evaluate_layout(
    anchors: list[Anchor],
    linker_models,
    population: np.lib.npyio.NpzFile,
    rng: np.random.Generator,
) -> dict[str, float]:
    anchor_xyz = anchor_array(anchors)
    diameters = population["diameter_nm"]
    receptor_counts = population["receptor_count"]
    receptor_points = population["receptor_points"]

    capture_values = []
    strong_values = []
    mean_contacts = []
    low_receptor_capture = []
    clustered_capture = []
    for ev_index, diameter in enumerate(diameters):
        receptor_count = int(receptor_counts[ev_index])
        if receptor_count <= 0:
            continue
        radius = float(diameter) / 2.0
        receptors_body = receptor_points[ev_index, :receptor_count, :]
        best_contacts = []
        best_p1 = 0.0
        best_p2 = 0.0
        for offset_x, offset_y in LATERAL_OFFSETS_NM:
            center = np.array([offset_x, offset_y, radius + SURFACE_CLEARANCE_NM])
            receptors = receptors_body + center
            distances = np.linalg.norm(anchor_xyz[:, None, :] - receptors[None, :, :], axis=2)
            probabilities = np.zeros_like(distances)
            for anchor_index, anchor in enumerate(anchors):
                construct = str(anchor["linker_construct"])
                probabilities[anchor_index, :] = [
                    reach_probability(linker_models, construct, float(distance))
                    for distance in distances[anchor_index, :]
                ]
            samples = [
                max_bipartite_matches(rng.random(probabilities.shape) < probabilities)
                for _ in range(N_BINDING_TRIALS)
            ]
            samples_arr = np.asarray(samples, dtype=float)
            p1 = float(np.mean(samples_arr >= 1))
            p2 = float(np.mean(samples_arr >= 2))
            if p1 + 0.25 * p2 > best_p1 + 0.25 * best_p2:
                best_p1 = p1
                best_p2 = p2
                best_contacts = samples
        mean_contact = float(np.mean(best_contacts)) if best_contacts else 0.0
        capture_values.append(best_p1)
        strong_values.append(best_p2)
        mean_contacts.append(mean_contact)
        if receptor_count <= 3:
            low_receptor_capture.append(best_p1)
        pattern = str(population["pattern"][ev_index])
        if "cluster" in pattern:
            clustered_capture.append(best_p1)

    capture = float(np.mean(capture_values))
    strong = float(np.mean(strong_values))
    mean_contact = float(np.mean(mean_contacts))
    low_capture = float(np.mean(low_receptor_capture)) if low_receptor_capture else 0.0
    cluster_capture = float(np.mean(clustered_capture)) if clustered_capture else 0.0
    penalty = layout_penalty(anchors)
    score = (
        0.45 * capture
        + 0.20 * strong
        + 0.15 * min(mean_contact / 2.0, 1.0)
        + 0.15 * low_capture
        + 0.05 * cluster_capture
        - penalty
    )
    return {
        "population_score": score,
        "capture_probability": capture,
        "p_at_least_2_contacts": strong,
        "mean_contacts": mean_contact,
        "low_receptor_capture": low_capture,
        "clustered_capture": cluster_capture,
        "crowding_penalty": penalty,
    }


def score_layouts(layouts: dict[str, list[Anchor]], linker_models, population, seed: int) -> list[dict[str, str]]:
    rows = []
    for index, (name, anchors) in enumerate(layouts.items()):
        rng = np.random.default_rng(seed + index)
        metrics = evaluate_layout(anchors, linker_models, population, rng)
        rows.append(
            {
                "layout": name,
                "aptamer_count": str(len(anchors)),
                **{key: f"{value:.5f}" for key, value in metrics.items()},
            }
        )
    return sorted(rows, key=lambda row: float(row["population_score"]), reverse=True)


def write_layouts(layouts: dict[str, list[Anchor]], top_names: list[str]) -> None:
    with open(OUT_LAYOUTS_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["layout", "anchor_id", "x_nm", "y_nm", "linker_construct"],
        )
        writer.writeheader()
        for name in top_names:
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


def plot_scores(rows: list[dict[str, str]]) -> None:
    top = rows[:12]
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    labels = [row["layout"] for row in top][::-1]
    scores = [float(row["population_score"]) for row in top][::-1]
    ax.barh(labels, scores, color="#1565c0")
    ax.set_xlabel("population score")
    ax.set_title("Best layouts against realistic sparse CD133+ EV population")
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)


def plot_layouts(layouts: dict[str, list[Anchor]], top_names: list[str]) -> None:
    n = min(6, len(top_names))
    fig, axes = plt.subplots(2, 3, figsize=(11, 7), constrained_layout=True)
    for ax, name in zip(axes.flat, top_names[:n]):
        anchors = layouts[name]
        xy = np.asarray([[float(a["x_nm"]), float(a["y_nm"])] for a in anchors], dtype=float)
        ax.add_patch(plt.Rectangle((-45, -30), 90, 60, fill=False, lw=1.2, color="0.25"))
        ax.scatter(xy[:, 0], xy[:, 1], s=42, color="#1565c0")
        ax.set_title(f"{name} ({len(anchors)})", fontsize=9)
        ax.set_xlim(-50, 50)
        ax.set_ylim(-35, 35)
        ax.set_aspect("equal")
        ax.grid(alpha=0.2)
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.savefig(OUT_LAYOUT_PNG, dpi=220)
    plt.close(fig)


def main() -> None:
    records, receptor_array = generate_population(N_POPULATION_EVS, RNG_SEED)
    np.savez_compressed(
        ROOT / "ev_population_optimization_run.npz",
        receptor_points=receptor_array,
        diameter_nm=np.asarray([r.diameter_nm for r in records], dtype=float),
        receptor_count=np.asarray([r.receptor_count for r in records], dtype=int),
        pattern=np.asarray([r.pattern for r in records]),
    )
    population = np.load(ROOT / "ev_population_optimization_run.npz")
    linker_models = load_linker_models()
    rng = np.random.default_rng(RNG_SEED)
    layouts = seed_layouts(rng)

    best_rows = score_layouts(layouts, linker_models, population, RNG_SEED + 1000)
    for generation in range(N_GENERATIONS):
        parents = [row["layout"] for row in best_rows[:KEEP_PER_GENERATION]]
        for parent_name in parents:
            parent = layouts[parent_name]
            for mutant_index in range(MUTANTS_PER_PARENT):
                child_name, child = mutate_layout(
                    f"evo_g{generation + 1}_{parent_name}_{mutant_index}",
                    parent,
                    rng,
                )
                layouts[child_name] = child
        best_rows = score_layouts(layouts, linker_models, population, RNG_SEED + 2000 + generation)
        print(
            f"Generation {generation + 1}: best={best_rows[0]['layout']} "
            f"score={best_rows[0]['population_score']}",
            flush=True,
        )

    with open(OUT_SCORES_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(best_rows[0].keys()))
        writer.writeheader()
        writer.writerows(best_rows)

    top_names = [row["layout"] for row in best_rows[:10]]
    write_layouts(layouts, top_names)
    plot_scores(best_rows)
    plot_layouts(layouts, top_names)

    summary = {
        "model": "EV population generator plus evolutionary layout competition",
        "rng_seed": RNG_SEED,
        "n_population_evs": N_POPULATION_EVS,
        "n_layouts_tested": len(layouts),
        "n_generations": N_GENERATIONS,
        "binding_trials_per_ev_offset": N_BINDING_TRIALS,
        "lateral_offsets_nm": LATERAL_OFFSETS_NM,
        "best_layout": best_rows[0],
        "top_10_layouts": best_rows[:10],
        "outputs": {
            "scores": OUT_SCORES_CSV.name,
            "top_layout_coordinates": OUT_LAYOUTS_CSV.name,
            "score_plot": OUT_PLOT.name,
            "layout_plot": OUT_LAYOUT_PNG.name,
        },
        "interpretation": [
            "This ranks layouts against a mixed EV population, not one ideal 73 nm vesicle.",
            "The score rewards at least one contact, two-contact events, low-receptor EV capture, and clustered-receptor capture.",
            "A small crowding penalty is included so overly packed layouts are not favored for free.",
        ],
    }
    OUT_SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")

    print(f"Wrote {OUT_SCORES_CSV.name}")
    print(f"Wrote {OUT_LAYOUTS_CSV.name}")
    print(f"Wrote {OUT_SUMMARY_JSON.name}")
    print(f"Wrote {OUT_PLOT.name}")
    print(f"Wrote {OUT_LAYOUT_PNG.name}")
    print(
        f"Best layout: {best_rows[0]['layout']} "
        f"score={best_rows[0]['population_score']} "
        f"capture={best_rows[0]['capture_probability']}"
    )


if __name__ == "__main__":
    main()
