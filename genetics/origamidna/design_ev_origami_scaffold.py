#!/usr/bin/env python3
"""Generate first-pass DNA-origami aptamer layouts for CD133+ EV capture."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT_CSV = ROOT / "ev_origami_aptamer_layouts.csv"
OUT_PNG = ROOT / "ev_origami_aptamer_layouts.png"

TILE_WIDTH_NM = 90.0
TILE_HEIGHT_NM = 60.0
DEFAULT_LINKER_REACH_NM = 15.0
DEFAULT_LINKER_CONSTRUCT = "polyT20"


@dataclass(frozen=True)
class Anchor:
    layout: str
    anchor_id: str
    x_nm: float
    y_nm: float
    linker_reach_nm: float = DEFAULT_LINKER_REACH_NM
    linker_construct: str = DEFAULT_LINKER_CONSTRUCT
    receptor: str = "CD133 aptamer"


def ring_points(layout: str, n: int, radius_nm: float) -> list[Anchor]:
    anchors = []
    for i in range(n):
        theta = 2.0 * math.pi * i / n
        anchors.append(
            Anchor(
                layout=layout,
                anchor_id=f"A{i + 1:02d}",
                x_nm=radius_nm * math.cos(theta),
                y_nm=radius_nm * math.sin(theta),
            )
        )
    return anchors


def grid_points(layout: str, cols: int, rows: int, x_span: float, y_span: float) -> list[Anchor]:
    anchors = []
    idx = 1
    for iy in range(rows):
        y = -y_span / 2.0 if rows == 1 else -y_span / 2.0 + iy * y_span / (rows - 1)
        for ix in range(cols):
            x = -x_span / 2.0 if cols == 1 else -x_span / 2.0 + ix * x_span / (cols - 1)
            anchors.append(
                Anchor(
                    layout=layout,
                    anchor_id=f"A{idx:02d}",
                    x_nm=x,
                    y_nm=y,
                )
            )
            idx += 1
    return anchors


def dense_points() -> list[Anchor]:
    anchors: list[Anchor] = []
    idx = 1
    for radius, n in [(14.0, 8), (27.0, 16)]:
        for i in range(n):
            theta = 2.0 * math.pi * (i + 0.5 * (radius > 14.0)) / n
            anchors.append(
                Anchor(
                    layout="dense_24",
                    anchor_id=f"A{idx:02d}",
                    x_nm=radius * math.cos(theta),
                    y_nm=0.75 * radius * math.sin(theta),
                )
            )
            idx += 1
    return anchors


def build_layouts() -> list[Anchor]:
    return (
        ring_points("sparse_6", 6, 24.0)
        + ring_points("ring_12", 12, 25.0)
        + grid_points("grid_18", 6, 3, 70.0, 36.0)
        + dense_points()
    )


def write_csv(anchors: list[Anchor]) -> None:
    with open(OUT_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "layout",
                "anchor_id",
                "x_nm",
                "y_nm",
                "z_nm",
                "linker_reach_nm",
                "linker_construct",
                "receptor",
            ],
        )
        writer.writeheader()
        for a in anchors:
            writer.writerow(
                {
                    "layout": a.layout,
                    "anchor_id": a.anchor_id,
                    "x_nm": f"{a.x_nm:.3f}",
                    "y_nm": f"{a.y_nm:.3f}",
                    "z_nm": "0.000",
                    "linker_reach_nm": f"{a.linker_reach_nm:.3f}",
                    "linker_construct": a.linker_construct,
                    "receptor": a.receptor,
                }
            )


def plot_layouts(anchors: list[Anchor]) -> None:
    layouts = ["sparse_6", "ring_12", "grid_18", "dense_24"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for ax, layout in zip(axes.flat, layouts):
        subset = [a for a in anchors if a.layout == layout]
        ax.add_patch(
            plt.Rectangle(
                (-TILE_WIDTH_NM / 2.0, -TILE_HEIGHT_NM / 2.0),
                TILE_WIDTH_NM,
                TILE_HEIGHT_NM,
                fill=False,
                lw=1.5,
                color="0.25",
            )
        )
        ax.scatter([a.x_nm for a in subset], [a.y_nm for a in subset], s=55, color="#1565c0")
        for a in subset:
            ax.text(a.x_nm, a.y_nm + 2.5, a.anchor_id, ha="center", va="bottom", fontsize=7)
        ax.set_title(f"{layout}: {len(subset)} aptamers")
        ax.set_xlim(-50, 50)
        ax.set_ylim(-36, 36)
        ax.set_aspect("equal")
        ax.set_xlabel("x on origami tile (nm)")
        ax.set_ylabel("y on origami tile (nm)")
        ax.grid(alpha=0.2)
    fig.suptitle("Candidate DNA-origami CD133 aptamer layouts", fontsize=14)
    fig.savefig(OUT_PNG, dpi=220)
    plt.close(fig)


def main() -> None:
    anchors = build_layouts()
    write_csv(anchors)
    plot_layouts(anchors)
    print(f"Wrote {OUT_CSV.name}")
    print(f"Wrote {OUT_PNG.name}")
    for layout in ["sparse_6", "ring_12", "grid_18", "dense_24"]:
        print(f"  {layout}: {sum(a.layout == layout for a in anchors)} aptamers")


if __name__ == "__main__":
    main()
