#!/usr/bin/env python3
"""Create oxDNA-ready poly-T linker and CD133 aptamer calibration simulations."""

from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATOR = ROOT / "oxDNA" / "utils" / "generate-sa.py"

LINKERS = {
    "polyT10": 10,
    "polyT15": 15,
    "polyT20": 20,
    "polyT30": 30,
}

# Shigdar A15 CD133 aptamer (Cancer Letters 2013, DOI 10.1016/j.canlet.2012.11.032).
# Sequence is RNA but written with T in seq.txt because oxDNA encodes T=U; the
# RNA force field is selected via interaction_type = RNA2 in the input file.
APTAMERS = {
    "A15": "CCCTCCTACATAGGG",
}

APTAMER_STEM_PAIRS = {
    "A15": [(0, 14), (1, 13), (2, 12), (3, 11)],
}


def _input_template(interaction_type: str, steps: str) -> str:
    return f"""########################################
####    PROGRAM PARAMETERS    ####
########################################
backend = CPU
seed = 20260525

########################################
####    SIM PARAMETERS        ####
########################################
interaction_type = {interaction_type}
salt_concentration = 1.0
sim_type = MD
ensemble = NVT
thermostat = john
newtonian_steps = 103
diff_coeff = 2.5
steps = {steps}
check_energy_every = 10000
check_energy_threshold = 1.e-4

T = 25C
dt = 0.005
verlet_skin = 0.20

########################################
####    INPUT / OUTPUT        ####
########################################
topology = generated.top
conf_file = generated.dat
trajectory_file = trajectory.dat
lastconf_file = last_conf.dat
refresh_vel = 1
no_stdout_energy = 1
restart_step_counter = 1
energy_file = energy.dat
print_conf_interval = 1000
print_energy_every = 1000
time_scale = linear
external_forces = 0
"""


LINKER_INPUT = _input_template("DNA2", "2e5")
APTAMER_INPUT = _input_template("RNA2", "5e5")
APTAMER_RESTRAINED_INPUT = (
    _input_template("RNA2", "5e5")
    .replace(
        "external_forces = 0",
        "external_forces = 1\nexternal_forces_file = a15_stem_restraints.dat",
    )
    .replace(
        "dt = 0.005",
        "dt = 0.001",
    )
)


def _mutual_trap_block(particle: int, ref_particle: int, stiff: float = 1.0, r0: float = 1.2) -> str:
    return f"""{{
type = mutual_trap
particle = {particle}
ref_particle = {ref_particle}
stiff = {stiff:g}
r0 = {r0:g}
PBC = 1
}}
"""


def _stem_restraints_text(pairs: list[tuple[int, int]]) -> str:
    blocks = []
    for left, right in pairs:
        blocks.append(_mutual_trap_block(left, right))
        blocks.append(_mutual_trap_block(right, left))
    return "\n".join(blocks)


def _generate(sim_dir: Path, sequence: str, input_text: str) -> None:
    sim_dir.mkdir(exist_ok=True)
    (sim_dir / "seq.txt").write_text(sequence + "\n", encoding="ascii")
    subprocess.run(
        ["python3", str(GENERATOR), "80", "seq.txt"],
        cwd=sim_dir,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    (sim_dir / "input").write_text(input_text, encoding="ascii")


def write_linker(name: str, bases: int) -> None:
    _generate(ROOT / f"sim_linker_{name}", "T" * bases, LINKER_INPUT)


def write_aptamer(name: str, sequence: str) -> None:
    _generate(ROOT / f"sim_aptamer_{name}", sequence, APTAMER_INPUT)
    restrained_dir = ROOT / f"sim_aptamer_{name}_restrained"
    _generate(restrained_dir, sequence, APTAMER_RESTRAINED_INPUT)
    (restrained_dir / "a15_stem_restraints.dat").write_text(
        _stem_restraints_text(APTAMER_STEM_PAIRS[name]), encoding="ascii"
    )


def main() -> None:
    for name, bases in LINKERS.items():
        write_linker(name, bases)
        print(f"Prepared sim_linker_{name}/ ({bases} thymidines, DNA2)")
    for name, sequence in APTAMERS.items():
        write_aptamer(name, sequence)
        print(f"Prepared sim_aptamer_{name}/ ({len(sequence)} nt, RNA2)")
        print(f"Prepared sim_aptamer_{name}_restrained/ ({len(sequence)} nt, RNA2 + stem restraints)")
    print("\nRun a calibration with, for example:")
    print("  cd sim_linker_polyT20")
    print("  ../../oxDNA/build/bin/oxDNA input")
    print("  cd ../sim_aptamer_A15")
    print("  ../../oxDNA/build/bin/oxDNA input")
    print("  cd ..")
    print("  python3 calibrate_linker_reach.py")


if __name__ == "__main__":
    main()
