#!/usr/bin/env python3
"""Minimal Brownian dynamics model for DNA-origami EV capture.

This is a lightweight bridge between the static reach/overlap scorer and a
future HOOMD-blue or LAMMPS model. It keeps the DNA origami tile fixed, treats
the EV as a Brownian sphere near the tile, places finite CD133 receptors on the
lower EV hemisphere, and lets one-to-one aptamer/receptor bonds form and break
over time.

The model is phenomenological: the rates below are not fitted kinetic
constants. The purpose is to compare relative dynamic stability for the current
lead design and matched controls under the same assumptions.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Union

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from score_ev_capture_geometry import (
    CD133_DENSITIES,
    RNG_SEED,
    load_layouts,
    load_linker_models,
    random_lower_hemisphere,
    reach_probability,
    receptor_count,
)

ROOT = Path(__file__).resolve().parent
OUT_TRAJ_CSV = ROOT / "ev_capture_dynamics_trajectories.csv"
OUT_SUMMARY_CSV = ROOT / "ev_capture_dynamics_summary.csv"
OUT_JSON = ROOT / "ev_capture_dynamics_summary.json"
OUT_PNG = ROOT / "ev_capture_dynamics_contacts.png"

LEAD_LAYOUT = "dense_24"
LEAD_LINKER = "polyT30"
EV_DIAMETERS_NM = (150.0,)
DENSITY_NAMES = ("medium", "high")
RANDOM_CONTROL_REPLICATES = 12
TRAJECTORIES_PER_CASE = 80
N_STEPS = 600
DT_SECONDS = 0.05
CAPTURE_CONTACT_THRESHOLD = 3
STRONG_CONTACT_THRESHOLD = 6
CAPTURE_DWELL_SECONDS = 1.0
EV_SURFACE_CLEARANCE_NM = 2.0
INITIAL_GAP_NM = 8.0
MAX_GAP_NM = 24.0
LATERAL_START_NM = 22.0
LATERAL_ESCAPE_NM = 75.0
D_FREE_NM2_PER_S = 35.0
D_BOUND_FLOOR_FRACTION = 0.10
K_ON_PER_STEP = 0.18
K_OFF_PER_S = 0.12
RNG_DYNAMIC_SEED = RNG_SEED + 303

Anchor = dict[str, Union[float, str]]
LinkerModels = dict[str, tuple[np.ndarray, np.ndarray]]


def with_linker(anchors: list[Anchor], construct: str) -> list[Anchor]:
    return [
        {
            **anchor,
            "linker_construct": construct,
            "linker_reach_nm": float(anchor.get("linker_reach_nm", 15.0)),
        }
        for anchor in anchors
    ]


def random_control_layout(aptamer_count: int, rng: np.random.Generator) -> list[Anchor]:
    xs = rng.uniform(-45.0, 45.0, size=aptamer_count)
    ys = rng.uniform(-30.0, 30.0, size=aptamer_count)
    return [
        {
            "x_nm": float(x),
            "y_nm": float(y),
            "linker_reach_nm": 15.0,
            "linker_construct": LEAD_LINKER,
        }
        for x, y in zip(xs, ys)
    ]


def anchor_array(anchors: list[Anchor]) -> np.ndarray:
    return np.asarray(
        [[float(anchor["x_nm"]), float(anchor["y_nm"]), 0.0] for anchor in anchors],
        dtype=float,
    )


def try_form_bonds(
    distances: np.ndarray,
    anchors: list[Anchor],
    linker_models: LinkerModels,
    anchor_to_receptor: np.ndarray,
    receptor_to_anchor: np.ndarray,
    binding_activity: float,
    rng: np.random.Generator,
) -> None:
    if binding_activity <= 0.0:
        return

    free_anchor_indices = np.flatnonzero(anchor_to_receptor == -1)
    free_receptor_indices = np.flatnonzero(receptor_to_anchor == -1)
    if len(free_anchor_indices) == 0 or len(free_receptor_indices) == 0:
        return

    sub_distances = distances[np.ix_(free_anchor_indices, free_receptor_indices)]
    probabilities = np.zeros_like(sub_distances)
    for local_i, anchor_index in enumerate(free_anchor_indices):
        construct = str(anchors[int(anchor_index)]["linker_construct"])
        probabilities[local_i, :] = [
            binding_activity * K_ON_PER_STEP * reach_probability(linker_models, construct, float(d))
            for d in sub_distances[local_i, :]
        ]

    candidate_pairs = np.argwhere(rng.random(probabilities.shape) < probabilities)
    if len(candidate_pairs) == 0:
        return

    rng.shuffle(candidate_pairs)
    for local_anchor, local_receptor in candidate_pairs:
        anchor_index = int(free_anchor_indices[local_anchor])
        receptor_index = int(free_receptor_indices[local_receptor])
        if anchor_to_receptor[anchor_index] == -1 and receptor_to_anchor[receptor_index] == -1:
            anchor_to_receptor[anchor_index] = receptor_index
            receptor_to_anchor[receptor_index] = anchor_index


def break_bonds(
    distances: np.ndarray,
    anchor_to_receptor: np.ndarray,
    receptor_to_anchor: np.ndarray,
    rng: np.random.Generator,
) -> None:
    p_off = 1.0 - math.exp(-K_OFF_PER_S * DT_SECONDS)
    for anchor_index, receptor_index in enumerate(anchor_to_receptor):
        if receptor_index == -1:
            continue
        distance = distances[anchor_index, receptor_index]
        strained = distance > 14.0
        if strained or rng.random() < p_off:
            anchor_to_receptor[anchor_index] = -1
            receptor_to_anchor[receptor_index] = -1


def simulate_trajectory(
    case_label: str,
    anchors: list[Anchor],
    linker_models: LinkerModels,
    ev_diameter_nm: float,
    density_name: str,
    binding_activity: float,
    rng: np.random.Generator,
    fixed_receptor_count: int | None = None,
) -> dict[str, object]:
    radius = ev_diameter_nm / 2.0
    n_receptors = (
        max(1, int(fixed_receptor_count))
        if fixed_receptor_count is not None
        else receptor_count(radius, CD133_DENSITIES[density_name])
    )
    receptor_body = random_lower_hemisphere(n_receptors, radius, rng)
    anchor_xyz = anchor_array(anchors)

    center = np.array(
        [
            rng.uniform(-LATERAL_START_NM, LATERAL_START_NM),
            rng.uniform(-LATERAL_START_NM, LATERAL_START_NM),
            radius + EV_SURFACE_CLEARANCE_NM + INITIAL_GAP_NM,
        ],
        dtype=float,
    )
    min_center_z = radius + EV_SURFACE_CLEARANCE_NM
    max_center_z = radius + EV_SURFACE_CLEARANCE_NM + MAX_GAP_NM
    anchor_to_receptor = np.full(len(anchors), -1, dtype=int)
    receptor_to_anchor = np.full(n_receptors, -1, dtype=int)

    contacts_trace: list[int] = []
    captured_trace: list[int] = []
    consecutive_capture_steps = 0
    required_capture_steps = max(1, int(round(CAPTURE_DWELL_SECONDS / DT_SECONDS)))
    ever_captured = False

    for _ in range(N_STEPS):
        n_contacts = int(np.count_nonzero(anchor_to_receptor != -1))
        mobility_scale = max(D_BOUND_FLOOR_FRACTION, 1.0 / (1.0 + 0.75 * n_contacts))
        step_sigma = math.sqrt(2.0 * D_FREE_NM2_PER_S * mobility_scale * DT_SECONDS)
        center[:2] += rng.normal(0.0, step_sigma, size=2)
        center[2] += rng.normal(0.0, step_sigma * 0.45)
        center[2] -= 0.035 * n_contacts
        center[2] = float(np.clip(center[2], min_center_z, max_center_z))

        if np.linalg.norm(center[:2]) > LATERAL_ESCAPE_NM:
            # Keep escaped particles in the bookkeeping with no new capture.
            center[:2] *= LATERAL_ESCAPE_NM / np.linalg.norm(center[:2])

        receptors = receptor_body + center
        distances = np.linalg.norm(anchor_xyz[:, None, :] - receptors[None, :, :], axis=2)
        break_bonds(distances, anchor_to_receptor, receptor_to_anchor, rng)
        try_form_bonds(
            distances,
            anchors,
            linker_models,
            anchor_to_receptor,
            receptor_to_anchor,
            binding_activity,
            rng,
        )

        n_contacts = int(np.count_nonzero(anchor_to_receptor != -1))
        contacts_trace.append(n_contacts)
        if n_contacts >= CAPTURE_CONTACT_THRESHOLD:
            consecutive_capture_steps += 1
        else:
            consecutive_capture_steps = 0
        captured_now = consecutive_capture_steps >= required_capture_steps
        ever_captured = ever_captured or captured_now
        captured_trace.append(int(captured_now))

    contacts = np.asarray(contacts_trace, dtype=float)
    captured = np.asarray(captured_trace, dtype=float)
    return {
        "case": case_label,
        "ev_diameter_nm": ev_diameter_nm,
        "cd133_density": density_name,
        "binding_activity": binding_activity,
        "receptor_count": n_receptors,
        "ever_captured": float(ever_captured),
        "final_contacts": float(contacts[-1]),
        "max_contacts": float(np.max(contacts)),
        "mean_contacts": float(np.mean(contacts)),
        "mean_last_quarter_contacts": float(np.mean(contacts[int(0.75 * len(contacts)) :])),
        "strong_contact_fraction": float(np.mean(contacts >= STRONG_CONTACT_THRESHOLD)),
        "capture_dwell_fraction": float(np.mean(captured)),
        "contacts_trace": contacts_trace,
    }


def summarize(results: list[dict[str, object]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, float, str], list[dict[str, object]]] = defaultdict(list)
    for result in results:
        grouped[
            (
                str(result["case"]),
                float(result["ev_diameter_nm"]),
                str(result["cd133_density"]),
            )
        ].append(result)

    rows: list[dict[str, str]] = []
    for (case, diameter, density), group in sorted(grouped.items()):
        rows.append(
            {
                "case": case,
                "trajectories": str(len(group)),
                "ev_diameter_nm": f"{diameter:.0f}",
                "cd133_density": density,
                "receptor_count": f"{float(group[0]['receptor_count']):.0f}",
                "capture_probability": f"{np.mean([g['ever_captured'] for g in group]):.4f}",
                "mean_contacts": f"{np.mean([g['mean_contacts'] for g in group]):.3f}",
                "mean_last_quarter_contacts": f"{np.mean([g['mean_last_quarter_contacts'] for g in group]):.3f}",
                "mean_max_contacts": f"{np.mean([g['max_contacts'] for g in group]):.3f}",
                "strong_contact_fraction": f"{np.mean([g['strong_contact_fraction'] for g in group]):.4f}",
                "capture_dwell_fraction": f"{np.mean([g['capture_dwell_fraction'] for g in group]):.4f}",
            }
        )
    return rows


def write_trajectory_csv(results: list[dict[str, object]]) -> None:
    with open(OUT_TRAJ_CSV, "w", newline="", encoding="ascii") as f:
        fieldnames = [
            "case",
            "trajectory_id",
            "ev_diameter_nm",
            "cd133_density",
            "time_s",
            "contacts",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for trajectory_id, result in enumerate(results):
            for step, contacts in enumerate(result["contacts_trace"]):
                writer.writerow(
                    {
                        "case": result["case"],
                        "trajectory_id": trajectory_id,
                        "ev_diameter_nm": f"{float(result['ev_diameter_nm']):.0f}",
                        "cd133_density": result["cd133_density"],
                        "time_s": f"{step * DT_SECONDS:.2f}",
                        "contacts": contacts,
                    }
                )


def write_summary_csv(rows: list[dict[str, str]]) -> None:
    with open(OUT_SUMMARY_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_contact_traces(results: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, len(DENSITY_NAMES), figsize=(10, 4), sharey=True, constrained_layout=True)
    if len(DENSITY_NAMES) == 1:
        axes = [axes]

    for ax, density in zip(axes, DENSITY_NAMES):
        for case in ("dense_24_polyT30", "random_24_polyT30", "scrambled_dense_24_polyT30"):
            traces = [
                np.asarray(result["contacts_trace"], dtype=float)
                for result in results
                if result["case"] == case and result["cd133_density"] == density
            ]
            if not traces:
                continue
            mat = np.vstack(traces)
            mean_trace = np.mean(mat, axis=0)
            ax.plot(np.arange(N_STEPS) * DT_SECONDS, mean_trace, label=case)
        ax.axhline(CAPTURE_CONTACT_THRESHOLD, color="black", lw=1.0, ls="--", alpha=0.5)
        ax.axhline(STRONG_CONTACT_THRESHOLD, color="black", lw=1.0, ls=":", alpha=0.5)
        ax.set_title(f"150 nm EV, {density} CD133")
        ax.set_xlabel("time (s)")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("mean simultaneous contacts")
    axes[-1].legend(fontsize=8)
    fig.savefig(OUT_PNG, dpi=220)
    plt.close(fig)


def main() -> None:
    rng = np.random.default_rng(RNG_DYNAMIC_SEED)
    layouts = load_layouts()
    linker_models = load_linker_models()
    dense_anchors = with_linker(layouts[LEAD_LAYOUT], LEAD_LINKER)
    random_layouts = [
        random_control_layout(len(dense_anchors), rng)
        for _ in range(RANDOM_CONTROL_REPLICATES)
    ]

    cases = [
        ("dense_24_polyT30", dense_anchors, 1.0),
        ("scrambled_dense_24_polyT30", dense_anchors, 0.0),
    ]
    for i, control in enumerate(random_layouts, start=1):
        cases.append((f"random_24_polyT30_rep{i:02d}", control, 1.0))

    results: list[dict[str, object]] = []
    for diameter in EV_DIAMETERS_NM:
        for density in DENSITY_NAMES:
            print(f"Simulating EV={diameter:.0f} density={density}", flush=True)
            for trajectory_index in range(TRAJECTORIES_PER_CASE):
                result = simulate_trajectory(
                    "dense_24_polyT30",
                    dense_anchors,
                    linker_models,
                    diameter,
                    density,
                    1.0,
                    rng,
                )
                results.append(result)

                result = simulate_trajectory(
                    "scrambled_dense_24_polyT30",
                    dense_anchors,
                    linker_models,
                    diameter,
                    density,
                    0.0,
                    rng,
                )
                results.append(result)

                random_case, random_anchors, random_activity = cases[
                    2 + (trajectory_index % RANDOM_CONTROL_REPLICATES)
                ]
                result = simulate_trajectory(
                    random_case,
                    random_anchors,
                    linker_models,
                    diameter,
                    density,
                    random_activity,
                    rng,
                )
                result["case"] = "random_24_polyT30"
                results.append(result)

    summary_rows = summarize(results)
    write_trajectory_csv(results)
    write_summary_csv(summary_rows)
    plot_contact_traces(results)

    best = max(summary_rows, key=lambda row: float(row["capture_probability"]))
    summary = {
        "model": "minimal Brownian EV capture dynamics",
        "rng_seed": RNG_DYNAMIC_SEED,
        "lead_design": "dense_24/polyT30",
        "ev_diameters_nm": list(EV_DIAMETERS_NM),
        "cd133_densities": list(DENSITY_NAMES),
        "trajectories_per_case": TRAJECTORIES_PER_CASE,
        "steps": N_STEPS,
        "dt_seconds": DT_SECONDS,
        "capture_definition": {
            "contact_threshold": CAPTURE_CONTACT_THRESHOLD,
            "required_dwell_seconds": CAPTURE_DWELL_SECONDS,
            "strong_contact_threshold": STRONG_CONTACT_THRESHOLD,
        },
        "kinetic_assumptions": {
            "k_on_per_step_multiplier": K_ON_PER_STEP,
            "k_off_per_s": K_OFF_PER_S,
            "free_diffusion_nm2_per_s": D_FREE_NM2_PER_S,
            "binding_activity_scrambled": 0.0,
        },
        "outputs": {
            "trajectory_contacts": OUT_TRAJ_CSV.name,
            "summary_table": OUT_SUMMARY_CSV.name,
            "mean_contact_plot": OUT_PNG.name,
        },
        "best_capture_probability_case": best,
        "summary_rows": summary_rows,
        "notes": [
            "This model is for relative design/control comparison, not physical kinetic parameter inference.",
            "Random_24 results are aggregated over 12 matched random dense 24-anchor layouts.",
            "Scrambled controls preserve the same geometry but set receptor-specific binding activity to zero.",
        ],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")

    print(f"Wrote {OUT_TRAJ_CSV.name}")
    print(f"Wrote {OUT_SUMMARY_CSV.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_PNG.name}")
    print(
        "Best dynamic case: "
        f"{best['case']} EV={best['ev_diameter_nm']} density={best['cd133_density']} "
        f"capture_probability={best['capture_probability']}"
    )


if __name__ == "__main__":
    main()
