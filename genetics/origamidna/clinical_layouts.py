#!/usr/bin/env python3
"""Clinical candidate aptamer layouts for 73 nm sparse CD133+ EVs.

These layouts keep the same 90 x 60 nm tile but move active aptamer anchors
toward the small central contact patch expected for a 73 nm vesicle.
"""

from __future__ import annotations

import math
from typing import Union

Anchor = dict[str, Union[float, str]]


def grid_layout(
    cols: int,
    rows: int,
    x_span_nm: float,
    y_span_nm: float,
    construct: str = "polyT30",
) -> list[Anchor]:
    anchors: list[Anchor] = []
    for iy in range(rows):
        y = 0.0 if rows == 1 else -y_span_nm / 2.0 + iy * y_span_nm / (rows - 1)
        for ix in range(cols):
            x = 0.0 if cols == 1 else -x_span_nm / 2.0 + ix * x_span_nm / (cols - 1)
            anchors.append(
                {
                    "x_nm": float(x),
                    "y_nm": float(y),
                    "linker_reach_nm": 15.0,
                    "linker_construct": construct,
                }
            )
    return anchors


def ellipse_ring(
    n: int,
    radius_nm: float,
    y_scale: float = 0.75,
    phase: float = 0.0,
    construct: str = "polyT30",
) -> list[Anchor]:
    anchors: list[Anchor] = []
    for i in range(n):
        theta = 2.0 * math.pi * (i + phase) / n
        anchors.append(
            {
                "x_nm": float(radius_nm * math.cos(theta)),
                "y_nm": float(y_scale * radius_nm * math.sin(theta)),
                "linker_reach_nm": 15.0,
                "linker_construct": construct,
            }
        )
    return anchors


def clinical_candidate_layouts(construct: str = "polyT30") -> dict[str, list[Anchor]]:
    """Return compact layouts tuned for tiny, receptor-sparse EVs."""
    hybrid_24 = (
        grid_layout(6, 3, 60.0, 30.0, construct)
        + grid_layout(3, 2, 20.0, 10.0, construct)
    )
    rescue_24 = (
        ellipse_ring(8, 13.0, 0.75, 0.0, construct)
        + ellipse_ring(8, 24.0, 0.75, 0.5, construct)
        + grid_layout(4, 2, 68.0, 28.0, construct)
    )
    return {
        # Same aptamer count as dense_24, but the anchors sit in a smaller patch.
        "clinical_dense_24": grid_layout(6, 4, 30.0, 18.0, construct),
        # Same count as dense_24, but combines broad "first contact" coverage
        # with a small central cluster for sparse clinical EVs.
        "clinical_hybrid_24": hybrid_24,
        "clinical_rescue_24": rescue_24,
        # Fewer active aptamers should reduce crowding while keeping coverage.
        "clinical_grid_18": grid_layout(6, 3, 30.0, 16.0, construct),
        "clinical_grid_12": grid_layout(4, 3, 24.0, 16.0, construct),
    }
