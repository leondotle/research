#!/usr/bin/env python3
"""Brownian dynamics stress test for 73 nm sparse CD133+ EV capture."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import simulate_ev_capture_dynamics as dyn
from clinical_layouts import clinical_candidate_layouts
from simulate_ev_capture_dynamics import (
    RNG_DYNAMIC_SEED,
    anchor_array,
    load_layouts,
    load_linker_models,
    random_control_layout,
    simulate_trajectory,
    with_linker,
)

ROOT = Path(__file__).resolve().parent
OUT_TRAJ_CSV = ROOT / "ev_capture_clinical_73nm_dynamics_trajectories.csv"
OUT_SUMMARY_CSV = ROOT / "ev_capture_clinical_73nm_dynamics_summary.csv"
OUT_JSON = ROOT / "ev_capture_clinical_73nm_dynamics_summary.json"
OUT_PNG = ROOT / "ev_capture_clinical_73nm_dynamics_contacts.png"

EV_DIAMETER_NM = 73.0
RECEPTOR_COUNTS = (2, 5, 10)
LINKER = "polyT30"
LAYOUTS_TO_TEST = (
    "sparse_6",
    "ring_12",
    "grid_18",
    "dense_24",
    "clinical_grid_12",
    "clinical_grid_18",
    "clinical_dense_24",
    "clinical_hybrid_24",
    "clinical_rescue_24",
)
RANDOM_CONTROL_REPLICATES = 12
TRAJECTORIES_PER_CASE = 120
CLINICAL_DYNAMIC_SEED = RNG_DYNAMIC_SEED + 404


def configure_sparse_capture() -> None:
    dyn.EV_DIAMETERS_NM = (EV_DIAMETER_NM,)
    dyn.CAPTURE_CONTACT_THRESHOLD = 1
    dyn.STRONG_CONTACT_THRESHOLD = 2
    dyn.CAPTURE_DWELL_SECONDS = 1.0
    dyn.N_STEPS = 600


def summarize(results: list[dict[str, object]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for result in results:
        grouped[(str(result["case"]), int(result["fixed_receptor_count"]))].append(result)

    rows: list[dict[str, str]] = []
    for (case, receptor_count), group in sorted(grouped.items()):
        rows.append(
            {
                "case": case,
                "trajectories": str(len(group)),
                "ev_diameter_nm": f"{EV_DIAMETER_NM:.0f}",
                "receptor_count": str(receptor_count),
                "capture_probability_p_ge_1_dwell": f"{np.mean([g['ever_captured'] for g in group]):.4f}",
                "mean_contacts": f"{np.mean([g['mean_contacts'] for g in group]):.3f}",
                "mean_last_quarter_contacts": f"{np.mean([g['mean_last_quarter_contacts'] for g in group]):.3f}",
                "mean_max_contacts": f"{np.mean([g['max_contacts'] for g in group]):.3f}",
                "strong_contact_fraction_p_ge_2": f"{np.mean([g['strong_contact_fraction'] for g in group]):.4f}",
                "capture_dwell_fraction": f"{np.mean([g['capture_dwell_fraction'] for g in group]):.4f}",
            }
        )
    return rows


def write_trajectory_csv(results: list[dict[str, object]]) -> None:
    with open(OUT_TRAJ_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "case",
                "trajectory_id",
                "ev_diameter_nm",
                "receptor_count",
                "time_s",
                "contacts",
            ],
        )
        writer.writeheader()
        for trajectory_id, result in enumerate(results):
            for step, contacts in enumerate(result["contacts_trace"]):
                writer.writerow(
                    {
                        "case": result["case"],
                        "trajectory_id": trajectory_id,
                        "ev_diameter_nm": f"{EV_DIAMETER_NM:.0f}",
                        "receptor_count": result["fixed_receptor_count"],
                        "time_s": f"{step * dyn.DT_SECONDS:.2f}",
                        "contacts": contacts,
                    }
                )


def write_summary_csv(rows: list[dict[str, str]]) -> None:
    with open(OUT_SUMMARY_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_contact_traces(results: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, len(RECEPTOR_COUNTS), figsize=(12, 4), sharey=True, constrained_layout=True)
    for ax, receptor_count in zip(axes, RECEPTOR_COUNTS):
        for case in (*LAYOUTS_TO_TEST, "random_24", "scrambled_dense_24"):
            traces = [
                np.asarray(result["contacts_trace"], dtype=float)
                for result in results
                if result["case"] == case and result["fixed_receptor_count"] == receptor_count
            ]
            if not traces:
                continue
            mean_trace = np.mean(np.vstack(traces), axis=0)
            ax.plot(np.arange(dyn.N_STEPS) * dyn.DT_SECONDS, mean_trace, label=case)
        ax.axhline(1, color="black", lw=1.0, ls="--", alpha=0.45)
        ax.axhline(2, color="black", lw=1.0, ls=":", alpha=0.45)
        ax.set_title(f"73 nm EV, {receptor_count} CD133")
        ax.set_xlabel("time (s)")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("mean simultaneous contacts")
    axes[-1].legend(fontsize=7)
    fig.savefig(OUT_PNG, dpi=220)
    plt.close(fig)


def main() -> None:
    configure_sparse_capture()
    rng = np.random.default_rng(CLINICAL_DYNAMIC_SEED)
    layouts = load_layouts()
    layouts.update(clinical_candidate_layouts(LINKER))
    linker_models = load_linker_models()

    designed_cases = [
        (layout, with_linker(layouts[layout], LINKER), 1.0) for layout in LAYOUTS_TO_TEST
    ]
    dense_anchors = with_linker(layouts["dense_24"], LINKER)
    random_layouts = [
        random_control_layout(len(dense_anchors), rng)
        for _ in range(RANDOM_CONTROL_REPLICATES)
    ]
    scrambled_case = ("scrambled_dense_24", dense_anchors, 0.0)

    results: list[dict[str, object]] = []
    for receptor_count in RECEPTOR_COUNTS:
        print(f"Simulating 73 nm EV with {receptor_count} CD133 receptors", flush=True)
        for trajectory_index in range(TRAJECTORIES_PER_CASE):
            for case_label, anchors, activity in designed_cases:
                result = simulate_trajectory(
                    case_label,
                    anchors,
                    linker_models,
                    EV_DIAMETER_NM,
                    f"count_{receptor_count}",
                    activity,
                    rng,
                    fixed_receptor_count=receptor_count,
                )
                result["fixed_receptor_count"] = receptor_count
                results.append(result)

            random_anchors = random_layouts[trajectory_index % RANDOM_CONTROL_REPLICATES]
            result = simulate_trajectory(
                "random_24",
                random_anchors,
                linker_models,
                EV_DIAMETER_NM,
                f"count_{receptor_count}",
                1.0,
                rng,
                fixed_receptor_count=receptor_count,
            )
            result["fixed_receptor_count"] = receptor_count
            results.append(result)

            result = simulate_trajectory(
                scrambled_case[0],
                scrambled_case[1],
                linker_models,
                EV_DIAMETER_NM,
                f"count_{receptor_count}",
                scrambled_case[2],
                rng,
                fixed_receptor_count=receptor_count,
            )
            result["fixed_receptor_count"] = receptor_count
            results.append(result)

    rows = summarize(results)
    write_trajectory_csv(results)
    write_summary_csv(rows)
    plot_contact_traces(results)

    best_by_receptor_count = {
        str(count): max(
            (row for row in rows if row["receptor_count"] == str(count)),
            key=lambda row: float(row["capture_probability_p_ge_1_dwell"]),
        )
        for count in RECEPTOR_COUNTS
    }
    summary = {
        "model": "clinical sparse 73 nm Brownian EV capture dynamics",
        "ev_diameter_nm": EV_DIAMETER_NM,
        "fixed_receptor_counts": list(RECEPTOR_COUNTS),
        "linker_construct": LINKER,
        "layouts_tested": list(LAYOUTS_TO_TEST),
        "random_control_replicates": RANDOM_CONTROL_REPLICATES,
        "trajectories_per_case": TRAJECTORIES_PER_CASE,
        "steps": dyn.N_STEPS,
        "dt_seconds": dyn.DT_SECONDS,
        "capture_definition": {
            "contact_threshold": dyn.CAPTURE_CONTACT_THRESHOLD,
            "required_dwell_seconds": dyn.CAPTURE_DWELL_SECONDS,
            "strong_contact_threshold": dyn.STRONG_CONTACT_THRESHOLD,
        },
        "rng_seed": CLINICAL_DYNAMIC_SEED,
        "outputs": {
            "trajectory_contacts": OUT_TRAJ_CSV.name,
            "summary_table": OUT_SUMMARY_CSV.name,
            "mean_contact_plot": OUT_PNG.name,
        },
        "best_by_receptor_count": best_by_receptor_count,
        "summary_rows": rows,
        "notes": [
            "Capture is redefined as at least one sustained contact for 1 s because the clinical receptor counts are only 2-10 per EV.",
            "Strong-contact fraction is the fraction of time with at least two simultaneous contacts.",
            "The model still lacks explicit steric/electrostatic repulsion, so dense layouts may remain optimistic.",
        ],
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")

    print(f"Wrote {OUT_TRAJ_CSV.name}")
    print(f"Wrote {OUT_SUMMARY_CSV.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_PNG.name}")
    for count, row in best_by_receptor_count.items():
        print(
            f"Best dynamic 73 nm count={count}: {row['case']} "
            f"capture={row['capture_probability_p_ge_1_dwell']} "
            f"mean_contacts={row['mean_contacts']}"
        )


if __name__ == "__main__":
    main()
