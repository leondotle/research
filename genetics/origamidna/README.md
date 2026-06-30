# CD133+ EV DNA-Origami Capture Modeling

This project uses computational models to design a DNA-origami surface for
capturing small, receptor-sparse CD133+ extracellular vesicles (EVs).

An EV is a tiny membrane particle released by a cell. An aptamer is a short
nucleic-acid binder that acts like a molecular hook. The project asks where
those hooks should be placed, how many are useful, and how capture changes
when real-world problems such as sparse receptors, Brownian motion, unwanted
sticking, damaged samples, and washing are included.

## Clinical Design Conditions

- Mean EV diameter: approximately 73-74 nm
- Sparse CD133 counts: commonly tested at 2, 5, and 10 receptors per EV
- Mixed receptor patterns: random, clustered, two-cluster, and bottom-cap
- Main capture definition: at least one aptamer-CD133 contact lasting 1 second
- Summer validation scope: synthetic targets in standard buffers

## Current Leading Design

The current computational lead is `triple_ring_24`, which places 24 aptamers
in three rings. The strongest modeled deployment combines:

- mixed or more upright spacers,
- an antifouling surface,
- staged capture zones,
- one moderate 30-second wash, and
- an additional modeled 30% target-contact improvement.

Under the current assumptions, the difficult-condition useful-capture score is
approximately 90.05%. This is a simulation result, not proof of laboratory or
clinical performance. Material effects, wash kinetics, and penalty values must
be replaced with experimental measurements when available.

## Main Workflow

1. `calibrate_linker_reach.py` estimates how far each aptamer linker can reach.
2. `ev_population_generator.py` creates varied clinical-like EV populations.
3. `optimize_population_layouts.py` searches aptamer positions.
4. `validate_population_dynamics.py` tests layouts with Brownian motion.
5. `robustness_sensitivity_screen.py` changes uncertain real-world conditions.
6. `advanced_capture_design_search.py` compares layouts and surface/spacer ideas.
7. `inverse_design_requirements.py` calculates what must improve to exceed 90%.
8. `optimize_capture_wash_protocol.py` searches passivation and wash protocols.

## Important Results

- `robustness_layout_summary.csv`: reliability of the original clinical layouts
- `advanced_capture_design_summary.csv`: joint layout and material search
- `inverse_design_requirements.csv`: improvements required for reliable 90%
- `capture_wash_protocol_best.csv`: best modeled wash protocol for each design
- `Research_Report_2026-06-20_Annotated.md`: report with figure-placement notes

## Running

The scripts require Python 3, NumPy, and Matplotlib. Run them from this folder,
for example:

```bash
python3 robustness_sensitivity_screen.py
python3 advanced_capture_design_search.py
python3 inverse_design_requirements.py
python3 optimize_capture_wash_protocol.py
```

Several scripts use fixed random seeds so results can be reproduced.

## Scientific Limits

- The model compares design hypotheses; it does not validate a medical device.
- Antifouling, spacer, damage, and wash parameters include explicit assumptions.
- A simulated 90% score is not equivalent to 90% experimental recovery.
- Patient aqueous-humor samples are outside the current validation stage.
- Synthetic-target controls are needed before any clinical interpretation.
