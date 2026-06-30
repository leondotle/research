#!/usr/bin/env python3
""""Validate the intended restrained A15 stem conformation from an oxRNA trajectory."""""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
SIGMA_NM = 0.8518

APTAMER_NAME = "A15"
STEM_PAIRS = [(0, 14), (1, 13), (2, 12), (3, 11)]
BINDING_FACE_INDICES = (4, 5, 6, 7, 8, 9, 10)
CONJUGATION_INDEX = 0
STEM_CLOSED_DISTANCE_NM = 1.5

PREFERRED_TRAJECTORY = ROOT / "sim_aptamer_A15_restrained" / "trajectory.dat"
FALLBACK_TRAJECTORY = ROOT / "sim_aptamer_A15" / "trajectory.dat"
OUT_JSON = ROOT / "aptamer_a15_conformation_summary.json"
OUT_CSV = ROOT / "aptamer_a15_stem_distances.csv"


def parse_trajectory(path: Path):
    with open(path, encoding="ascii") as f:
        frame_positions = []
        box = None
        timestep = None
        for line in f:
            line = line.strip()
            if line.startswith("t = "):
                if frame_positions:
                    yield timestep, box, np.asarray(frame_positions, dtype=float)
                timestep = int(line.split("=")[1])
                frame_positions = []
            elif line.startswith("b = "):
                box = np.asarray([float(v) for v in line.split("=")[1].split()[:3]], dtype=float)
            elif line.startswith("E = "):
                continue
            elif line:
                vals = line.split()
                if len(vals) >= 3:
                    frame_positions.append([float(vals[0]), float(vals[1]), float(vals[2])])
        if frame_positions:
            yield timestep, box, np.asarray(frame_positions, dtype=float)


def minimum_image_vector(delta: np.ndarray, box: np.ndarray | None) -> np.ndarray:
    if box is None:
        return delta
    return delta - box * np.rint(delta / box)


def unwrap_backbone(positions: np.ndarray, box: np.ndarray | None) -> np.ndarray:
    if box is None or len(positions) == 0:
        return positions
    unwrapped = np.empty_like(positions)
    unwrapped[0] = positions[0]
    for i in range(1, len(positions)):
        unwrapped[i] = unwrapped[i - 1] + minimum_image_vector(positions[i] - positions[i - 1], box)
    return unwrapped


def trajectory_source() -> tuple[Path, str]:
    if PREFERRED_TRAJECTORY.exists():
        return PREFERRED_TRAJECTORY, "restrained"
    if FALLBACK_TRAJECTORY.exists():
        return FALLBACK_TRAJECTORY, "unrestrained"
    raise FileNotFoundError("No A15 trajectory found. Run the A15 oxRNA simulation first.")


def main() -> None:
    trajectory, source = trajectory_source()
    rows: list[dict[str, str]] = []
    all_pair_distances: list[list[float]] = []
    head_reaches: list[float] = []
    closed_frames = 0

    for frame_index, (timestep, box, positions) in enumerate(parse_trajectory(trajectory), start=1):
        if len(positions) < 15:
            continue
        unwrapped = unwrap_backbone(positions[:15], box)
        pair_distances = [
            float(np.linalg.norm(unwrapped[left] - unwrapped[right]) * SIGMA_NM)
            for left, right in STEM_PAIRS
        ]
        all_pair_distances.append(pair_distances)
        closed = all(distance <= STEM_CLOSED_DISTANCE_NM for distance in pair_distances)
        closed_frames += int(closed)
        anchor = unwrapped[CONJUGATION_INDEX]
        binding_face = unwrapped[list(BINDING_FACE_INDICES)].mean(axis=0)
        head_reach = float(np.linalg.norm(binding_face - anchor) * SIGMA_NM)
        head_reaches.append(head_reach)
        rows.append(
            {
                "frame": str(frame_index),
                "timestep": str(timestep),
                **{
                    f"pair_{left + 1}_{right + 1}_distance_nm": f"{distance:.4f}"
                    for (left, right), distance in zip(STEM_PAIRS, pair_distances)
                },
                "all_stem_pairs_closed": str(closed).lower(),
                "head_reach_nm": f"{head_reach:.4f}",
            }
        )

    if not all_pair_distances:
        raise RuntimeError(f"No usable frames found in {trajectory}")

    pair_array = np.asarray(all_pair_distances, dtype=float)
    head_array = np.asarray(head_reaches, dtype=float)
    summary = {
        "aptamer": APTAMER_NAME,
        "trajectory": str(trajectory.relative_to(ROOT)),
        "source": source,
        "stem_pairs_1based": [[left + 1, right + 1] for left, right in STEM_PAIRS],
        "stem_closed_distance_nm": STEM_CLOSED_DISTANCE_NM,
        "n_frames": int(len(pair_array)),
        "stem_closed_fraction": float(closed_frames / len(pair_array)),
        "pair_distance_mean_nm": {
            f"{left + 1}-{right + 1}": float(pair_array[:, i].mean())
            for i, (left, right) in enumerate(STEM_PAIRS)
        },
        "pair_distance_p90_nm": {
            f"{left + 1}-{right + 1}": float(np.percentile(pair_array[:, i], 90))
            for i, (left, right) in enumerate(STEM_PAIRS)
        },
        "head_reach_mean_nm": float(head_array.mean()),
        "head_reach_p10_nm": float(np.percentile(head_array, 10)),
        "head_reach_p50_nm": float(np.percentile(head_array, 50)),
        "head_reach_p90_nm": float(np.percentile(head_array, 90)),
    }

    with open(OUT_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")

    print(f"Wrote {OUT_CSV.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(
        f"A15 ({source}): closed_fraction={summary['stem_closed_fraction']:.3f}, "
        f"head_mean={summary['head_reach_mean_nm']:.3f} nm"
    )


if __name__ == "__main__":
    main()
