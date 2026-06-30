#!/usr/bin/env python3
"""Estimate aptamer linker reach distributions for EV capture scoring.

The reach of a tethered aptamer is modeled as the magnitude of the vector sum
of two independent isotropic contributions:

    composite reach = | r_linker + r_aptamer_head |

where r_linker is the polyT linker end-to-end vector (5' anchor on origami to
the 3' end conjugated to the aptamer) and r_aptamer_head is the vector from
the aptamer's 5' conjugation point to its binding face. Both magnitudes are
sampled from oxDNA / oxRNA trajectories when available.

If oxDNA trajectories for the linker exist (sim_linker_*/trajectory.dat) they
are used; otherwise a transparent worm-like-chain fallback for ssDNA is used.
If the restrained A15 aptamer trajectory exists
(sim_aptamer_A15_restrained/trajectory.dat) the head reach distribution is
taken from it. Otherwise the unrestrained A15 trajectory is used when present,
and the script falls back to a 3 nm scalar only if no aptamer trajectory exists.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT_CSV = ROOT / "linker_reach_models.csv"
OUT_JSON = ROOT / "linker_reach_summary.json"
OUT_PNG = ROOT / "linker_reach_survival.png"

SIGMA_NM = 0.8518
SSDNA_CONTOUR_PER_BASE_NM = 0.60
SSDNA_PERSISTENCE_LENGTH_NM = 1.50
PLACEHOLDER_APTAMER_HEAD_REACH_NM = 3.0

LINKERS = {
    "polyT10": 10,
    "polyT15": 15,
    "polyT20": 20,
    "polyT30": 30,
}

# A15 is a 4-bp stem (positions 1-4 paired with 12-15) and a 7-nt apical loop
# (positions 5-11). The 5' end (position 1) is the conjugation point to the
# linker; the binding face sits in the loop. Using the loop centroid as a
# representative binding-face position is more physically defensible than
# either 5'-to-3' end-to-end (which collapses for a closed hairpin) or the
# maximum atom distance (which is dominated by occasional fluctuations).
APTAMER_NAME = "A15"
APTAMER_LENGTH = 15
APTAMER_CONJUGATION_INDEX = 0  # 5' nucleotide
APTAMER_BINDING_FACE_INDICES = (4, 5, 6, 7, 8, 9, 10)  # loop region, 0-indexed
MONTE_CARLO_SAMPLES = 50000
RNG_SEED = 20260525

def parse_trajectory(path: Path):
    with open(path, encoding="ascii") as f:
        frame_positions = []
        box = None
        timestep = None
        for line in f:
            line = line.strip()
            if line.startswith("t = "):
                if frame_positions:
                    yield timestep, box, np.asarray(frame_positions)
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
            yield timestep, box, np.asarray(frame_positions)

def minimum_image_vector(delta: np.ndarray, box: np.ndarray | None) -> np.ndarray:
    if box is None:
        return delta
    return delta - box * np.rint(delta / box)

def unwrap_backbone(pos: np.ndarray, box: np.ndarray | None) -> np.ndarray:
    if box is None or len(pos) == 0:
        return pos
    unwrapped = np.empty_like(pos)
    unwrapped[0] = pos[0]
    for i in range(1, len(pos)):
        unwrapped[i] = unwrapped[i - 1] + minimum_image_vector(pos[i] - pos[i - 1], box)
    return unwrapped

def linker_distances_nm(name: str, bases: int) -> np.ndarray | None:
    traj = ROOT / f"sim_linker_{name}" / "trajectory.dat"
    if not traj.exists():
        return None
    distances = []
    for _, box, pos in parse_trajectory(traj):
        if len(pos) >= bases:
            unwrapped = unwrap_backbone(pos[:bases], box)
            distances.append(np.linalg.norm(unwrapped[bases - 1] - unwrapped[0]) * SIGMA_NM)
    if len(distances) < 10:
        return None
    burn_in = max(5, len(distances) // 10)
    return np.asarray(distances[burn_in:], dtype=float)

def aptamer_head_distances_nm() -> tuple[np.ndarray | None, str]:
    candidates = [
        (ROOT / f"sim_aptamer_{APTAMER_NAME}_restrained" / "trajectory.dat", "oxRNA_restrained_trajectory"),
        (ROOT / f"sim_aptamer_{APTAMER_NAME}" / "trajectory.dat", "oxRNA_unrestrained_trajectory"),
    ]
    traj = None
    source = "scalar_placeholder"
    for candidate, candidate_source in candidates:
        if candidate.exists():
            traj = candidate
            source = candidate_source
            break
    if traj is None:
        return None, "scalar_placeholder"
    distances = []
    for _, box, pos in parse_trajectory(traj):
        if len(pos) < APTAMER_LENGTH:
            continue
        unwrapped = unwrap_backbone(pos[:APTAMER_LENGTH], box)
        anchor = unwrapped[APTAMER_CONJUGATION_INDEX]
        binding_face = unwrapped[list(APTAMER_BINDING_FACE_INDICES)].mean(axis=0)
        distances.append(np.linalg.norm(binding_face - anchor) * SIGMA_NM)
    if len(distances) < 10:
        return None, "scalar_placeholder"
    burn_in = max(5, len(distances) // 10)
    return np.asarray(distances[burn_in:], dtype=float), source

def wlc_rms_nm(bases: int) -> float:
    contour = bases * SSDNA_CONTOUR_PER_BASE_NM
    lp = SSDNA_PERSISTENCE_LENGTH_NM
    mean_square = 2.0 * lp * contour * (1.0 - (lp / contour) * (1.0 - math.exp(-contour / lp)))
    return math.sqrt(max(mean_square, 0.0))

def wlc_sample_distances_nm(bases: int, n: int, rng: np.random.Generator) -> np.ndarray:
    """Draw end-to-end distance samples from a Maxwell-like WLC approximation."""
    rms = wlc_rms_nm(bases)
    if rms == 0.0:
        return np.zeros(n)
    sigma = rms / math.sqrt(3.0)
    components = rng.normal(loc=0.0, scale=sigma, size=(n, 3))
    return np.linalg.norm(components, axis=1)

def composite_survival(
    linker_samples: np.ndarray,
    head_samples: np.ndarray | None,
    head_scalar_nm: float,
    distance_nm: float,
    rng: np.random.Generator,
) -> float:
    n = MONTE_CARLO_SAMPLES
    l = rng.choice(linker_samples, size=n, replace=True)
    if head_samples is not None and len(head_samples) > 0:
        h = rng.choice(head_samples, size=n, replace=True)
    else:
        h = np.full(n, head_scalar_nm)
    cos_theta = rng.uniform(-1.0, 1.0, size=n)
    composite = np.sqrt(np.clip(l * l + h * h + 2.0 * l * h * cos_theta, 0.0, None))
    return float(np.mean(composite >= distance_nm))

def composite_samples(
    linker_samples: np.ndarray,
    head_samples: np.ndarray | None,
    head_scalar_nm: float,
    rng: np.random.Generator,
) -> np.ndarray:
    n = MONTE_CARLO_SAMPLES
    l = rng.choice(linker_samples, size=n, replace=True)
    if head_samples is not None and len(head_samples) > 0:
        h = rng.choice(head_samples, size=n, replace=True)
    else:
        h = np.full(n, head_scalar_nm)
    cos_theta = rng.uniform(-1.0, 1.0, size=n)
    return np.sqrt(np.clip(l * l + h * h + 2.0 * l * h * cos_theta, 0.0, None))

def build_models() -> tuple[list[dict[str, str]], dict[str, object]]:
    rng = np.random.default_rng(RNG_SEED)
    head_samples, head_source = aptamer_head_distances_nm()
    head_reference_nm = (
        float(np.mean(head_samples)) if head_samples is not None else PLACEHOLDER_APTAMER_HEAD_REACH_NM
    )

    rows: list[dict[str, str]] = []
    for name, bases in LINKERS.items():
        observed_linker = linker_distances_nm(name, bases)
        if observed_linker is not None:
            linker_source = "oxDNA_trajectory"
            linker_samples = observed_linker
        else:
            linker_source = "wlc_analytic_fallback"
            linker_samples = wlc_sample_distances_nm(bases, MONTE_CARLO_SAMPLES, rng)

        composite = composite_samples(
            linker_samples, head_samples, PLACEHOLDER_APTAMER_HEAD_REACH_NM, rng
        )
        rms = float(np.sqrt(np.mean(composite**2)))
        mean = float(np.mean(composite))
        p10, p50, p90 = np.percentile(composite, [10, 50, 90])

        survival_values = {}
        for d in range(2, 31, 2):
            p = float(np.mean(composite >= float(d)))
            survival_values[f"p_reach_{d:02d}_nm"] = f"{p:.6f}"

        rows.append(
            {
                "construct": name,
                "bases": str(bases),
                "linker_source": linker_source,
                "head_source": head_source,
                "head_aptamer": APTAMER_NAME,
                "contour_plus_head_nm": f"{bases * SSDNA_CONTOUR_PER_BASE_NM + head_reference_nm:.3f}",
                "rms_reach_nm": f"{rms:.3f}",
                "mean_reach_nm": f"{mean:.3f}",
                "p10_reach_nm": f"{p10:.3f}",
                "p50_reach_nm": f"{p50:.3f}",
                "p90_reach_nm": f"{p90:.3f}",
                **survival_values,
            }
        )
    head_summary = {
        "aptamer": APTAMER_NAME,
        "sequence": "CCCUCCUACAUAGGG",
        "source": head_source,
        "conjugation_index_0based": APTAMER_CONJUGATION_INDEX,
        "binding_face_indices_0based": list(APTAMER_BINDING_FACE_INDICES),
    }
    if head_samples is not None:
        head_summary.update(
            {
                "n_frames": int(len(head_samples)),
                "rms_nm": float(np.sqrt(np.mean(head_samples**2))),
                "mean_nm": float(np.mean(head_samples)),
                "p10_nm": float(np.percentile(head_samples, 10)),
                "p50_nm": float(np.percentile(head_samples, 50)),
                "p90_nm": float(np.percentile(head_samples, 90)),
            }
        )
    else:
        head_summary["fallback_scalar_nm"] = PLACEHOLDER_APTAMER_HEAD_REACH_NM

    return rows, head_summary

def write_csv(rows: list[dict[str, str]]) -> None:
    fieldnames = list(rows[0].keys())
    with open(OUT_CSV, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

def plot_survival(rows: list[dict[str, str]]) -> None:
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    xs = list(range(2, 31, 2))
    for row in rows:
        ys = [float(row[f"p_reach_{d:02d}_nm"]) for d in xs]
        label = f"{row['construct']} + {row['head_aptamer']} ({row['head_source']})"
        ax.plot(xs, ys, marker="o", label=label)
    ax.set_xlabel("required reach from origami anchor (nm)")
    ax.set_ylabel("P(linker + aptamer-head proxy can reach)")
    ax.set_ylim(-0.03, 1.03)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    ax.set_title("polyT linker + CD133 A15 aptamer tether reach")
    fig.savefig(OUT_PNG, dpi=220)
    plt.close(fig)

def main() -> None:
    rows, head_summary = build_models()
    write_csv(rows)
    plot_survival(rows)
    summary = {
        "aptamer_head": head_summary,
        "ssDNA_contour_per_base_nm": SSDNA_CONTOUR_PER_BASE_NM,
        "ssDNA_persistence_length_nm": SSDNA_PERSISTENCE_LENGTH_NM,
        "composition_model": "isotropic-orientation vector sum of linker and aptamer-head reach",
        "monte_carlo_samples": MONTE_CARLO_SAMPLES,
        "constructs": rows,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
    print(f"Wrote {OUT_CSV.name}")
    print(f"Wrote {OUT_JSON.name}")
    print(f"Wrote {OUT_PNG.name}")
    print()
    print(
        f"Aptamer head ({APTAMER_NAME}): source={head_summary['source']}, "
        f"mean={head_summary.get('mean_nm', 'NA')} nm"
    )
    for row in rows:
        print(
            f"  {row['construct']} + {APTAMER_NAME}: linker={row['linker_source']}, "
            f"rms={row['rms_reach_nm']} nm, mean={row['mean_reach_nm']} nm"
        )

if __name__ == "__main__":
    main()
