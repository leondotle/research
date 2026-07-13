#!/usr/bin/env python3
"""Map aptamer layouts onto a simplified DNA-origami lattice.

Beginner picture:
The earlier layout is like drawing dots anywhere on paper.
This script changes that into placing dots on a real-ish DNA-origami grid.

Definitions:
* helix row: one parallel DNA double-helix track in the origami tile.
* base position: a location along that helix.
* lattice site: one allowed attachment point, described by helix row and base.
* snapping: moving a free dot to the nearest allowed lattice site.

This is still a simplified model. It is not a full caDNAno design and it does
not create staple sequences. Its job is to make the capture layout more honest:
aptamers can only sit on plausible DNA-origami attachment sites.
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

import score_ev_capture_geometry as scg
from advanced_capture_design_search import candidate_layouts
from score_ev_capture_clinical_73nm import enrich_sparse_metrics, sparse_score
from score_ev_capture_geometry import load_linker_models, score_layout
from score_origami_buildability import buildability_metrics

ROOT = Path(__file__).resolve().parent
OUT_SITES = ROOT / "origami_lattice_sites.csv"
OUT_MAPPED = ROOT / "origami_lattice_mapped_layouts.csv"
OUT_COMPARE = ROOT / "origami_lattice_mapping_comparison.csv"
OUT_JSON = ROOT / "origami_lattice_mapping_summary.json"
OUT_PLOT = ROOT / "origami_lattice_mapping.png"

TILE_WIDTH_NM = 90.0
TILE_HEIGHT_NM = 60.0

# Simplified origami grid assumptions.
# 2.5 nm is a common scale for spacing between neighboring DNA helices.
# 3.57 nm is about one full DNA helical turn, which is a convenient repeat for
# placing modifications that should point in a similar direction.
HELIX_SPACING_NM = 2.5
ATTACHMENT_REPEAT_NM = 3.57
EDGE_MARGIN_NM = 5.0

TARGET_LAYOUTS = ("broad_grid_24", "broad_grid_20", "triple_ring_24")
EV_DIAMETER_NM = 73.0
RECEPTOR_COUNTS = (2, 5, 10)
RNG_SEED = 20260711

Anchor = dict[str, float | str]


def make_lattice_sites() -> list[dict[str, float | int]]:
    """Create allowed attachment sites on a simple rectangular origami tile."""
    x_min = -TILE_WIDTH_NM / 2.0 + EDGE_MARGIN_NM
    x_max = TILE_WIDTH_NM / 2.0 - EDGE_MARGIN_NM
    y_min = -TILE_HEIGHT_NM / 2.0 + EDGE_MARGIN_NM
    y_max = TILE_HEIGHT_NM / 2.0 - EDGE_MARGIN_NM

    y_values = np.arange(y_min, y_max + 0.001, HELIX_SPACING_NM)
    x_values = np.arange(x_min, x_max + 0.001, ATTACHMENT_REPEAT_NM)

    sites: list[dict[str, float | int]] = []
    for helix_id, y in enumerate(y_values):
        for base_index, x in enumerate(x_values):
            sites.append(
                {
                    "site_id": len(sites) + 1,
                    "helix_id": int(helix_id),
                    "base_index": int(base_index),
                    "x_nm": float(x),
                    "y_nm": float(y),
                }
            )
    return sites


def nearest_available_site(
    x: float,
    y: float,
    sites: list[dict[str, float | int]],
    used_site_ids: set[int],
) -> dict[str, float | int]:
    available = [site for site in sites if int(site["site_id"]) not in used_site_ids]
    return min(
        available,
        key=lambda site: (float(site["x_nm"]) - x) ** 2 + (float(site["y_nm"]) - y) ** 2,
    )


def snap_layout_to_lattice(
    layout_name: str,
    anchors: list[Anchor],
    sites: list[dict[str, float | int]],
) -> tuple[list[Anchor], list[dict[str, str]]]:
    mapped: list[Anchor] = []
    rows: list[dict[str, str]] = []
    used_site_ids: set[int] = set()

    for anchor_id, anchor in enumerate(anchors, 1):
        original_x = float(anchor["x_nm"])
        original_y = float(anchor["y_nm"])
        site = nearest_available_site(original_x, original_y, sites, used_site_ids)
        used_site_ids.add(int(site["site_id"]))

        snapped_x = float(site["x_nm"])
        snapped_y = float(site["y_nm"])
        shift = math.hypot(snapped_x - original_x, snapped_y - original_y)

        mapped_anchor: Anchor = {
            "x_nm": snapped_x,
            "y_nm": snapped_y,
            "linker_reach_nm": float(anchor.get("linker_reach_nm", 15.0)),
            "linker_construct": str(anchor.get("linker_construct", "polyT30")),
            "helix_id": int(site["helix_id"]),
            "base_index": int(site["base_index"]),
            "site_id": int(site["site_id"]),
        }
        mapped.append(mapped_anchor)
        rows.append(
            {
                "layout": layout_name,
                "mapped_layout": f"{layout_name}_lattice",
                "anchor_id": str(anchor_id),
                "original_x_nm": f"{original_x:.3f}",
                "original_y_nm": f"{original_y:.3f}",
                "snapped_x_nm": f"{snapped_x:.3f}",
                "snapped_y_nm": f"{snapped_y:.3f}",
                "shift_nm": f"{shift:.3f}",
                "helix_id": str(int(site["helix_id"])),
                "base_index": str(int(site["base_index"])),
                "site_id": str(int(site["site_id"])),
                "linker_construct": str(mapped_anchor["linker_construct"]),
            }
        )
    return mapped, rows


def clinical_sparse_average(layout: list[Anchor]) -> dict[str, float]:
    """Fast 73 nm sparse-EV geometry check for the mapped layout."""
    scg.LATERAL_STEP_NM = 8.0
    scg.N_RECEPTOR_REALIZATIONS = 3
    scg.N_BINDING_TRIALS = 6

    rng = np.random.default_rng(RNG_SEED)
    linker_models = load_linker_models()
    scores = []
    mean_contacts = []
    p1_values = []
    p2_values = []
    p3_values = []
    for receptor_count in RECEPTOR_COUNTS:
        metrics = score_layout(
            layout,
            linker_models,
            EV_DIAMETER_NM,
            0.0,
            rng,
            fixed_receptor_count=receptor_count,
        )
        enriched = enrich_sparse_metrics(metrics)
        scores.append(enriched["clinical_sparse_score"])
        mean_contacts.append(enriched["mean_contacts"])
        p1_values.append(enriched["p_at_least_1_contact"])
        p2_values.append(enriched["p_at_least_2_contacts"])
        p3_values.append(enriched["p_at_least_3_contacts"])

    return {
        "mean_clinical_sparse_score": float(np.mean(scores)),
        "mean_contacts": float(np.mean(mean_contacts)),
        "mean_p_at_least_1_contact": float(np.mean(p1_values)),
        "mean_p_at_least_2_contacts": float(np.mean(p2_values)),
        "mean_p_at_least_3_contacts": float(np.mean(p3_values)),
    }


def fmt(value: float) -> str:
    if math.isclose(value, round(value)):
        return str(int(round(value)))
    return f"{value:.4f}"


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with open(path, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_mapping(
    original_layouts: dict[str, list[Anchor]],
    mapped_layouts: dict[str, list[Anchor]],
) -> None:
    fig, axes = plt.subplots(len(TARGET_LAYOUTS), 2, figsize=(9, 10), constrained_layout=True)
    for row_index, layout_name in enumerate(TARGET_LAYOUTS):
        original = np.asarray(
            [[float(anchor["x_nm"]), float(anchor["y_nm"])] for anchor in original_layouts[layout_name]]
        )
        mapped = np.asarray(
            [[float(anchor["x_nm"]), float(anchor["y_nm"])] for anchor in mapped_layouts[layout_name]]
        )
        for col_index, (title, xy) in enumerate(
            (("free coordinate layout", original), ("snapped lattice layout", mapped))
        ):
            ax = axes[row_index, col_index]
            ax.add_patch(plt.Rectangle((-45, -30), 90, 60, fill=False, color="0.3"))
            ax.scatter(xy[:, 0], xy[:, 1], s=38, color="#187a72")
            ax.set_title(f"{layout_name}\n{title}")
            ax.set_xlim(-50, 50)
            ax.set_ylim(-35, 35)
            ax.set_aspect("equal")
            ax.grid(alpha=0.2)
    fig.savefig(OUT_PLOT, dpi=220)
    plt.close(fig)


def main() -> None:
    sites = make_lattice_sites()
    all_layouts = candidate_layouts()
    selected = {name: all_layouts[name] for name in TARGET_LAYOUTS}

    site_rows = [
        {
            "site_id": str(int(site["site_id"])),
            "helix_id": str(int(site["helix_id"])),
            "base_index": str(int(site["base_index"])),
            "x_nm": f"{float(site['x_nm']):.3f}",
            "y_nm": f"{float(site['y_nm']):.3f}",
        }
        for site in sites
    ]

    mapped_rows: list[dict[str, str]] = []
    comparison_rows: list[dict[str, str]] = []
    mapped_layouts: dict[str, list[Anchor]] = {}

    for layout_name, original in selected.items():
        mapped, rows = snap_layout_to_lattice(layout_name, original, sites)
        mapped_layouts[layout_name] = mapped
        mapped_rows.extend(rows)

        original_build = buildability_metrics(original)
        mapped_build = buildability_metrics(mapped)
        original_sparse = clinical_sparse_average(original)
        mapped_sparse = clinical_sparse_average(mapped)
        shifts = np.asarray([float(row["shift_nm"]) for row in rows])

        comparison_rows.append(
            {
                "layout": layout_name,
                "mapped_layout": f"{layout_name}_lattice",
                "aptamer_count": str(len(original)),
                "mean_snap_shift_nm": f"{float(np.mean(shifts)):.3f}",
                "max_snap_shift_nm": f"{float(np.max(shifts)):.3f}",
                "original_origami_buildability": f"{original_build['origami_buildability_score']:.4f}",
                "mapped_origami_buildability": f"{mapped_build['origami_buildability_score']:.4f}",
                "original_min_spacing_nm": f"{original_build['min_spacing_nm']:.3f}",
                "mapped_min_spacing_nm": f"{mapped_build['min_spacing_nm']:.3f}",
                "original_sparse_score": f"{original_sparse['mean_clinical_sparse_score']:.4f}",
                "mapped_sparse_score": f"{mapped_sparse['mean_clinical_sparse_score']:.4f}",
                "sparse_score_change": f"{mapped_sparse['mean_clinical_sparse_score'] - original_sparse['mean_clinical_sparse_score']:.4f}",
                "mapped_mean_contacts": f"{mapped_sparse['mean_contacts']:.4f}",
                "mapped_p_at_least_1_contact": f"{mapped_sparse['mean_p_at_least_1_contact']:.4f}",
                "mapped_p_at_least_2_contacts": f"{mapped_sparse['mean_p_at_least_2_contacts']:.4f}",
                "mapped_p_at_least_3_contacts": f"{mapped_sparse['mean_p_at_least_3_contacts']:.4f}",
            }
        )

    comparison_rows.sort(
        key=lambda row: (
            float(row["mapped_origami_buildability"]),
            float(row["mapped_sparse_score"]),
        ),
        reverse=True,
    )

    write_csv(OUT_SITES, site_rows)
    write_csv(OUT_MAPPED, mapped_rows)
    write_csv(OUT_COMPARE, comparison_rows)
    plot_mapping(selected, mapped_layouts)

    summary = {
        "model": "simplified DNA-origami lattice mapper",
        "plain_language_summary": [
            "The old layout allowed aptamers to sit at any x/y coordinate.",
            "The new mapped layout moves each aptamer to the nearest allowed helix/base site.",
            "A small snap shift means the original pattern was already close to a realistic origami grid.",
        ],
        "lattice_assumptions": {
            "tile_width_nm": TILE_WIDTH_NM,
            "tile_height_nm": TILE_HEIGHT_NM,
            "helix_spacing_nm": HELIX_SPACING_NM,
            "attachment_repeat_nm": ATTACHMENT_REPEAT_NM,
            "edge_margin_nm": EDGE_MARGIN_NM,
            "site_count": len(sites),
        },
        "best_mapped_layout": comparison_rows[0],
        "ranked_mapped_layouts": comparison_rows,
        "limitations": [
            "This does not design staple sequences.",
            "This does not validate folding with full oxDNA or caDNAno routing.",
            "The lattice uses simplified helix/base spacing and should be replaced by a lab-specific origami design later.",
        ],
        "outputs": {
            "lattice_sites": OUT_SITES.name,
            "mapped_layout_coordinates": OUT_MAPPED.name,
            "mapping_comparison": OUT_COMPARE.name,
            "mapping_plot": OUT_PLOT.name,
        },
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")

    best = comparison_rows[0]
    print(f"Wrote {OUT_SITES.name}")
    print(f"Wrote {OUT_MAPPED.name}")
    print(f"Wrote {OUT_COMPARE.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_PLOT.name}")
    print(
        f"Best mapped layout: {best['mapped_layout']} "
        f"buildability={best['mapped_origami_buildability']} "
        f"sparse_score={best['mapped_sparse_score']} "
        f"mean_shift={best['mean_snap_shift_nm']} nm"
    )


if __name__ == "__main__":
    main()
