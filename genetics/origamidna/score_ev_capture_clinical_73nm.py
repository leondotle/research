#!/usr/bin/env python3
"""Clinical-realism EV capture screen for sparse 73 nm CD133+ EVs.

This keeps the original optimistic density-based screen intact and adds a
patient-AH-inspired stress test:

* EV diameter fixed at 73 nm
* exact CD133 receptor counts of 2, 5, and 10 per EV
* sparse-capture metrics focused on at least 1, 2, and 3 contacts
* linker sensitivity across polyT10, polyT15, polyT20, and polyT30
* matched random controls for each designed aptamer count

The outputs are intended to answer whether dense_24/polyT30 still survives
when the vesicle is small, highly curved, and receptor-limited.
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

import score_ev_capture_geometry as scg
from clinical_layouts import clinical_candidate_layouts
from score_ev_capture_geometry import RNG_SEED, load_layouts, load_linker_models, score_layout

ROOT = Path(__file__).resolve().parent
OUT_SCORES_CSV = ROOT / "ev_capture_clinical_73nm_scores.csv"
OUT_RANDOM_CSV = ROOT / "ev_capture_clinical_73nm_random_controls.csv"
OUT_SUMMARY_JSON = ROOT / "ev_capture_clinical_73nm_summary.json"
OUT_PLOT = ROOT / "ev_capture_clinical_73nm_scores.png"

EV_DIAMETER_NM = 73.0
RECEPTOR_COUNTS = (2, 5, 10)
LINKER_CONSTRUCTS = ("polyT10", "polyT15", "polyT20", "polyT30")
RANDOM_CONTROL_REPLICATES = 12
SCORING_LATERAL_STEP_NM = 8.0
SCORING_RECEPTOR_REALIZATIONS = 4
SCORING_BINDING_TRIALS = 8
TILE_WIDTH_NM = 90.0
TILE_HEIGHT_NM = 60.0
CLINICAL_SEED = RNG_SEED + 404

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


def sparse_score(metrics: dict[str, float]) -> float:
    # In the 2-10 receptor regime, durable 1-3 contact formation is more
    # meaningful than the old >=6-contact fabrication screen.
    mean_contacts = metrics["mean_contacts"]
    p1 = metrics["p_at_least_1_contact"]
    p2 = metrics["p_at_least_2_contacts"]
    p3 = metrics["p_at_least_3_contacts"]
    normalized_contacts = min(mean_contacts / 3.0, 1.0)
    return 0.40 * p1 + 0.30 * p2 + 0.20 * p3 + 0.10 * normalized_contacts


def enrich_sparse_metrics(metrics: dict[str, float]) -> dict[str, float]:
    enriched = dict(metrics)
    enriched["clinical_sparse_score"] = float(sparse_score(enriched))
    return enriched


def metrics_to_row(
    label: str,
    aptamer_count: int,
    construct: str,
    receptor_count: int,
    metrics: dict[str, float],
) -> dict[str, str]:
    return {
        "layout": label,
        "aptamer_count": str(aptamer_count),
        "linker_construct": construct,
        "ev_diameter_nm": f"{EV_DIAMETER_NM:.0f}",
        "receptor_count": str(receptor_count),
        "mean_contacts": f"{metrics['mean_contacts']:.3f}",
        "max_contacts": f"{metrics['max_contacts']:.0f}",
        "p_at_least_1_contact": f"{metrics['p_at_least_1_contact']:.4f}",
        "p_at_least_2_contacts": f"{metrics['p_at_least_2_contacts']:.4f}",
        "p_at_least_3_contacts": f"{metrics['p_at_least_3_contacts']:.4f}",
        "p_at_least_6_contacts": f"{metrics['p_at_least_6_contacts']:.4f}",
        "legacy_capture_score": f"{metrics['capture_score']:.4f}",
        "clinical_sparse_score": f"{metrics['clinical_sparse_score']:.4f}",
    }


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with open(path, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def summarize_replicates(rows: list[dict[str, float]]) -> dict[str, float]:
    keys = [
        "receptor_count",
        "mean_contacts",
        "max_contacts",
        "p_at_least_1_contact",
        "p_at_least_2_contacts",
        "p_at_least_3_contacts",
        "p_at_least_6_contacts",
        "capture_score",
        "clinical_sparse_score",
    ]
    summary = {key: float(np.mean([row[key] for row in rows])) for key in keys}
    summary["max_contacts"] = float(np.max([row["max_contacts"] for row in rows]))
    summary["clinical_sparse_score_sd"] = float(
        np.std([row["clinical_sparse_score"] for row in rows], ddof=0)
    )
    return summary


def plot_scores(rows: list[dict[str, str]]) -> None:
    best_linker_rows = {}
    for row in rows:
        key = (row["layout"], row["receptor_count"])
        if key not in best_linker_rows or float(row["clinical_sparse_score"]) > float(
            best_linker_rows[key]["clinical_sparse_score"]
        ):
            best_linker_rows[key] = row

    layouts = sorted({row["layout"] for row in rows})
    x = np.arange(len(RECEPTOR_COUNTS))
    width = 0.18
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    for i, layout in enumerate(layouts):
        ys = [
            float(best_linker_rows[(layout, str(count))]["clinical_sparse_score"])
            for count in RECEPTOR_COUNTS
        ]
        ax.bar(x + (i - 1.5) * width, ys, width=width, label=layout)
    ax.set_xticks(x)
    ax.set_xticklabels([str(count) for count in RECEPTOR_COUNTS])
    ax.set_xlabel("CD133 receptors per 73 nm EV")
    ax.set_ylabel("clinical sparse score")
    ax.set_ylim(0, 1)
    ax.set_title("Best-linker layout score under clinical sparse-EV assumptions")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(title="layout")
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)


def main() -> None:
    scg.LATERAL_STEP_NM = SCORING_LATERAL_STEP_NM
    scg.N_RECEPTOR_REALIZATIONS = SCORING_RECEPTOR_REALIZATIONS
    scg.N_BINDING_TRIALS = SCORING_BINDING_TRIALS

    rng = np.random.default_rng(CLINICAL_SEED)
    layouts = load_layouts()
    layouts.update(clinical_candidate_layouts())
    linker_models = load_linker_models()

    rows: list[dict[str, str]] = []
    for layout, anchors in sorted(layouts.items()):
        for construct in LINKER_CONSTRUCTS:
            variant = with_linker(anchors, construct)
            for receptor_count in RECEPTOR_COUNTS:
                metrics = score_layout(
                    variant,
                    linker_models,
                    EV_DIAMETER_NM,
                    0.0,
                    rng,
                    fixed_receptor_count=receptor_count,
                )
                metrics = enrich_sparse_metrics(metrics)
                rows.append(
                    metrics_to_row(layout, len(variant), construct, receptor_count, metrics)
                )

    random_rows: list[dict[str, str]] = []
    for layout, anchors in sorted(layouts.items()):
        aptamer_count = len(anchors)
        for construct in LINKER_CONSTRUCTS:
            controls = [
                random_control_layout(aptamer_count, construct, rng)
                for _ in range(RANDOM_CONTROL_REPLICATES)
            ]
            for receptor_count in RECEPTOR_COUNTS:
                replicate_metrics = []
                for control in controls:
                    metrics = score_layout(
                        control,
                        linker_models,
                        EV_DIAMETER_NM,
                        0.0,
                        rng,
                        fixed_receptor_count=receptor_count,
                    )
                    replicate_metrics.append(enrich_sparse_metrics(metrics))
                metrics = summarize_replicates(replicate_metrics)
                matching_design = max(
                    (
                        row
                        for row in rows
                        if row["layout"] == layout
                        and row["linker_construct"] == construct
                        and row["receptor_count"] == str(receptor_count)
                    ),
                    key=lambda row: float(row["clinical_sparse_score"]),
                )
                random_row = metrics_to_row(
                    f"random_matched_{aptamer_count}",
                    aptamer_count,
                    construct,
                    receptor_count,
                    metrics,
                )
                random_row["replicates"] = str(RANDOM_CONTROL_REPLICATES)
                random_row["matched_design"] = layout
                random_row["matched_design_score"] = matching_design["clinical_sparse_score"]
                random_row["score_delta_vs_matched_design"] = (
                    f"{float(matching_design['clinical_sparse_score']) - metrics['clinical_sparse_score']:.4f}"
                )
                random_row["clinical_sparse_score_sd"] = f"{metrics['clinical_sparse_score_sd']:.4f}"
                random_rows.append(random_row)

    write_csv(OUT_SCORES_CSV, rows)
    write_csv(OUT_RANDOM_CSV, random_rows)
    plot_scores(rows)

    best_by_receptor_count = {
        str(count): max(
            (row for row in rows if row["receptor_count"] == str(count)),
            key=lambda row: float(row["clinical_sparse_score"]),
        )
        for count in RECEPTOR_COUNTS
    }
    dense24_polyT30 = [
        row
        for row in rows
        if row["layout"] == "dense_24" and row["linker_construct"] == "polyT30"
    ]
    summary = {
        "model": "clinical sparse 73 nm EV static reach screen",
        "ev_diameter_nm": EV_DIAMETER_NM,
        "fixed_receptor_counts": list(RECEPTOR_COUNTS),
        "linker_constructs": list(LINKER_CONSTRUCTS),
        "random_control_replicates": RANDOM_CONTROL_REPLICATES,
        "scoring_resolution": {
            "lateral_step_nm": SCORING_LATERAL_STEP_NM,
            "receptor_realizations": SCORING_RECEPTOR_REALIZATIONS,
            "binding_trials_per_realization": SCORING_BINDING_TRIALS,
        },
        "rng_seed": CLINICAL_SEED,
        "outputs": {
            "designed_scores": OUT_SCORES_CSV.name,
            "random_controls": OUT_RANDOM_CSV.name,
            "score_plot": OUT_PLOT.name,
        },
        "best_by_receptor_count": best_by_receptor_count,
        "dense24_polyT30_rows": dense24_polyT30,
        "interpretation": [
            "Exact receptor counts replace the previous density-derived receptor estimates.",
            "Clinical sparse score emphasizes P>=1, P>=2, and P>=3 contacts because P>=6 is usually unrealistic with 2-10 receptors.",
            "Sparse contact probabilities are sampled from the same finite one-to-one aptamer/CD133 trials as the original p>=3 and p>=6 metrics.",
        ],
    }
    OUT_SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")

    print(f"Wrote {OUT_SCORES_CSV.name}")
    print(f"Wrote {OUT_RANDOM_CSV.name}")
    print(f"Wrote {OUT_SUMMARY_JSON.name}")
    print(f"Wrote {OUT_PLOT.name}")
    for count, row in best_by_receptor_count.items():
        print(
            f"Best static 73 nm count={count}: {row['layout']} {row['linker_construct']} "
            f"score={row['clinical_sparse_score']} mean_contacts={row['mean_contacts']}"
        )


if __name__ == "__main__":
    main()
