#!/usr/bin/env python3
"""Generate realistic sparse CD133+ EV populations for layout screening.

Beginner picture:
Each EV is a tiny ball. CD133 molecules are tiny handles on that ball. This
script makes many fake-but-realistic balls with different sizes, handle counts,
and handle patterns.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
OUT_CSV = ROOT / "ev_population_clinical.csv"
OUT_NPZ = ROOT / "ev_population_clinical.npz"
OUT_JSON = ROOT / "ev_population_clinical_summary.json"

RNG_SEED = 20260615
N_EVS = 160
MEAN_DIAMETER_NM = 73.5
DIAMETER_SD_NM = 12.0
MIN_DIAMETER_NM = 45.0
MAX_DIAMETER_NM = 115.0
MAX_RECEPTORS = 20


@dataclass(frozen=True)
class EVRecord:
    ev_id: int
    diameter_nm: float
    receptor_count: int
    pattern: str


def sample_diameter(rng: np.random.Generator) -> float:
    for _ in range(100):
        diameter = rng.normal(MEAN_DIAMETER_NM, DIAMETER_SD_NM)
        if MIN_DIAMETER_NM <= diameter <= MAX_DIAMETER_NM:
            return float(diameter)
    return float(np.clip(diameter, MIN_DIAMETER_NM, MAX_DIAMETER_NM))


def sample_receptor_count(rng: np.random.Generator) -> int:
    draw = rng.random()
    if draw < 0.30:
        return int(rng.integers(1, 4))
    if draw < 0.80:
        return int(rng.integers(4, 11))
    return int(rng.integers(11, 16))


def random_lower_hemisphere(n: int, radius_nm: float, rng: np.random.Generator) -> np.ndarray:
    theta = rng.uniform(0.0, 2.0 * math.pi, size=n)
    z_unit = -rng.uniform(0.0, 1.0, size=n)
    xy_unit = np.sqrt(np.clip(1.0 - z_unit * z_unit, 0.0, None))
    return np.column_stack(
        (
            radius_nm * xy_unit * np.cos(theta),
            radius_nm * xy_unit * np.sin(theta),
            radius_nm * z_unit,
        )
    )


def random_direction_lower(rng: np.random.Generator) -> np.ndarray:
    point = random_lower_hemisphere(1, 1.0, rng)[0]
    return point / np.linalg.norm(point)


def tangent_basis(direction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    helper = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(direction, helper))) > 0.90:
        helper = np.array([1.0, 0.0, 0.0])
    u = np.cross(direction, helper)
    u /= np.linalg.norm(u)
    v = np.cross(direction, u)
    v /= np.linalg.norm(v)
    return u, v


def cluster_points(
    n: int,
    radius_nm: float,
    rng: np.random.Generator,
    center: np.ndarray | None = None,
    spread_rad: float = 0.20,
) -> np.ndarray:
    center_dir = random_direction_lower(rng) if center is None else center / np.linalg.norm(center)
    u, v = tangent_basis(center_dir)
    points = []
    for _ in range(n):
        offset = rng.normal(0.0, spread_rad, size=2)
        direction = center_dir + offset[0] * u + offset[1] * v
        if direction[2] > 0:
            direction[2] *= -1
        direction /= np.linalg.norm(direction)
        points.append(radius_nm * direction)
    return np.asarray(points, dtype=float)


def receptor_points(pattern: str, receptor_count: int, radius_nm: float, rng: np.random.Generator) -> np.ndarray:
    if receptor_count == 0:
        return np.zeros((0, 3), dtype=float)
    if pattern == "random":
        return random_lower_hemisphere(receptor_count, radius_nm, rng)
    if pattern == "single_cluster":
        return cluster_points(receptor_count, radius_nm, rng, spread_rad=0.16)
    if pattern == "two_cluster":
        if receptor_count < 2:
            return cluster_points(receptor_count, radius_nm, rng, spread_rad=0.14)
        n1 = receptor_count // 2
        n2 = receptor_count - n1
        first = random_direction_lower(rng)
        second = random_direction_lower(rng)
        return np.vstack(
            (
                cluster_points(n1, radius_nm, rng, first, spread_rad=0.14),
                cluster_points(n2, radius_nm, rng, second, spread_rad=0.14),
            )
        )
    if pattern == "bottom_cap":
        return cluster_points(
            receptor_count,
            radius_nm,
            rng,
            center=np.array([0.0, 0.0, -1.0]),
            spread_rad=0.28,
        )
    raise ValueError(f"Unknown EV receptor pattern: {pattern}")


def generate_population(n_evs: int = N_EVS, seed: int = RNG_SEED) -> tuple[list[EVRecord], np.ndarray]:
    rng = np.random.default_rng(seed)
    patterns = ("random", "single_cluster", "two_cluster", "bottom_cap")
    pattern_prob = (0.40, 0.25, 0.20, 0.15)
    records: list[EVRecord] = []
    receptor_array = np.full((n_evs, MAX_RECEPTORS, 3), np.nan, dtype=float)
    for ev_id in range(n_evs):
        diameter = sample_diameter(rng)
        receptor_count = sample_receptor_count(rng)
        pattern = str(rng.choice(patterns, p=pattern_prob))
        points = receptor_points(pattern, receptor_count, diameter / 2.0, rng)
        records.append(EVRecord(ev_id, diameter, receptor_count, pattern))
        receptor_array[ev_id, :receptor_count, :] = points
    return records, receptor_array


def main() -> None:
    records, receptor_array = generate_population()
    with open(OUT_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["ev_id", "diameter_nm", "receptor_count", "pattern"],
        )
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "ev_id": record.ev_id,
                    "diameter_nm": f"{record.diameter_nm:.3f}",
                    "receptor_count": record.receptor_count,
                    "pattern": record.pattern,
                }
            )
    np.savez_compressed(
        OUT_NPZ,
        receptor_points=receptor_array,
        diameter_nm=np.asarray([r.diameter_nm for r in records], dtype=float),
        receptor_count=np.asarray([r.receptor_count for r in records], dtype=int),
        pattern=np.asarray([r.pattern for r in records]),
    )
    counts = np.asarray([r.receptor_count for r in records], dtype=float)
    diameters = np.asarray([r.diameter_nm for r in records], dtype=float)
    summary = {
        "model": "clinical CD133+ EV population generator",
        "n_evs": len(records),
        "rng_seed": RNG_SEED,
        "diameter_nm": {
            "mean": float(np.mean(diameters)),
            "p10": float(np.percentile(diameters, 10)),
            "p50": float(np.percentile(diameters, 50)),
            "p90": float(np.percentile(diameters, 90)),
        },
        "receptor_count": {
            "mean": float(np.mean(counts)),
            "p10": float(np.percentile(counts, 10)),
            "p50": float(np.percentile(counts, 50)),
            "p90": float(np.percentile(counts, 90)),
        },
        "pattern_counts": {
            pattern: sum(1 for record in records if record.pattern == pattern)
            for pattern in sorted({r.pattern for r in records})
        },
        "outputs": {
            "metadata_csv": OUT_CSV.name,
            "receptor_points_npz": OUT_NPZ.name,
        },
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
    print(f"Wrote {OUT_CSV.name}")
    print(f"Wrote {OUT_NPZ.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(
        f"Population: diameter mean={summary['diameter_nm']['mean']:.1f} nm, "
        f"receptor mean={summary['receptor_count']['mean']:.1f}"
    )


if __name__ == "__main__":
    main()
