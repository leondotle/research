# CD133+ EV DNA-Origami Capture Modeling

This project uses computer simulations to design a DNA-origami surface that
captures small, receptor-sparse CD133+ extracellular vesicles (EVs).

An **EV** is a tiny membrane particle released by a cell. An **aptamer** is a
short nucleic-acid binder that acts like a molecular hook. The simulations ask:

- Where should the aptamers be placed?
- How many aptamers should be used?
- What happens when an EV has only a few CD133 receptors?
- How do motion, unwanted sticking, sample loss, and washing affect capture?

## 1. Download the Project

Open Terminal and clone the GitHub repository:

```bash
git clone https://github.com/leondotle/research.git
```

Move into this project folder:

```bash
cd research/genetics/origamidna
```

Confirm that you are in the correct folder:

```bash
pwd
ls
```

You should see files such as `README.md`, `clinical_layouts.py`, and
`advanced_capture_design_search.py`.

## 2. Create a Python Environment

A virtual environment keeps this project's Python packages separate from other
projects on your computer.

Create it:

```bash
python3 -m venv .venv
```

Activate it on macOS or Linux:

```bash
source .venv/bin/activate
```

Install the two required packages:

```bash
python3 -m pip install --upgrade pip
python3 -m pip install numpy matplotlib
```

Check the installation:

```bash
python3 -c "import numpy, matplotlib; print('Setup works')"
```

## 3. Understand the Main Result

The current computational lead is `triple_ring_24`:

- `triple_ring` means the aptamers are arranged in three rings.
- `24` means one tile has 24 aptamer attachment points.
- The modeled surface uses mixed or more upright spacers.
- An antifouling treatment reduces unwanted sticking.
- The modeled protocol uses staged capture zones and one moderate wash.

The best stress-tested protocol produced approximately **90.05% useful
capture** under the current assumptions. This is a simulation result, not a
measured laboratory recovery rate.

## 4. Quickly Reproduce the Latest Analysis

The repository already contains the large upstream result tables. Therefore,
the quickest check only needs the final two analysis scripts.

Calculate what must improve to exceed 90% useful capture:

```bash
python3 inverse_design_requirements.py
```

Expected final message:

```text
Best path: triple_ring_24 + mixed_spacer_antifouling
```

Optimize passivation and washing:

```bash
python3 optimize_capture_wash_protocol.py
```

Expected final result:

```text
Best protocol: triple_ring_24 + mixed_spacer_antifouling
```

The exact result should show a difficult-condition useful-capture score close
to `0.9005`, which means 90.05% in the model.

## 5. Run the Main Portable Pipeline

Run these commands in order for a broader rebuild. Some searches take several
minutes because they simulate many EVs and uncertain conditions.

### Step 5.1: Create the Original Aptamer Layouts

```bash
python3 design_ev_origami_scaffold.py
```

This creates a table and image showing the original aptamer patterns.

### Step 5.2: Score the Original Geometry

```bash
python3 score_ev_capture_geometry.py
```

This performs a fast geometric test. It asks whether an aptamer can physically
reach a receptor in a still image of the EV and tile.

### Step 5.3: Test the 73 nm Clinical Conditions

```bash
python3 score_ev_capture_clinical_73nm.py
python3 simulate_ev_capture_clinical_73nm.py
```

These scripts use approximately 73 nm EVs with sparse receptor counts such as
2, 5, and 10. The second command adds Brownian motion, meaning random motion of
tiny particles in liquid.

The dynamics command creates a large trajectory CSV. Git ignores that file
because it can be regenerated.

### Step 5.4: Generate a Mixed EV Population

```bash
python3 ev_population_generator.py
```

This creates EVs with different sizes, receptor counts, and receptor patterns.
For example, one EV may have two scattered receptors while another has ten
receptors grouped into a cluster.

### Step 5.5: Search for Better Aptamer Positions

```bash
python3 optimize_population_layouts.py
```

This tests many layouts against the mixed EV population. It also creates the
temporary compressed population file required by the next command.

### Step 5.6: Validate the Best Layouts with Motion

```bash
python3 validate_population_dynamics.py
```

This checks whether the best layouts still work while the EV moves instead of
remaining frozen in one position.

### Step 5.7: Stress-Test Real-World Uncertainty

```bash
python3 robustness_sensitivity_screen.py
```

This deliberately changes uncertain conditions, including:

- EV diameter
- receptor count and clustering
- active aptamer percentage
- linker reach
- binding and unbinding rates
- unwanted sticking
- sample damage

A robust design works across many combinations, not only one ideal condition.

### Step 5.8: Search Layouts and Surface Formulations

```bash
python3 advanced_capture_design_search.py
```

This compares eight layouts and four modeled spacer/surface formulations. The
material multipliers are design assumptions, not measured material constants.

### Step 5.9: Calculate the Requirements for 90%

```bash
python3 inverse_design_requirements.py
```

This works backward from the results. Instead of asking, "What capture score
did we get?", it asks, "How much background reduction and binding improvement
are required to exceed 90%?"

### Step 5.10: Optimize Capture and Washing

```bash
python3 optimize_capture_wash_protocol.py
```

This searches wash strength, duration, cycle count, passivation, and target
contact improvement. Wash strength is a relative scale until laboratory flow
rates and bond lifetimes are measured.

### Step 5.11: Check DNA-Origami Buildability

```bash
python3 score_origami_buildability.py
```

This asks whether a layout is practical as a DNA-origami pattern, not just
whether it catches EVs in the model.

Simple example: two patterns may both catch EVs well, but a pattern with
aptamers spaced about 10 nm apart is easier to build and explain than a pattern
with several aptamers only 5-6 nm apart.

The current best practical lead is:

```text
Best practical layout: broad_grid_24
```

`broad_grid_24` keeps modeled useful capture above 90% while giving the
aptamer anchors more comfortable spacing on the DNA-origami tile.

### Step 5.12: Snap the Layout to an Origami Lattice

```bash
python3 map_layout_to_origami_lattice.py
```

This turns the best layouts from free x/y dots into simplified DNA-origami
attachment sites.

Simple example: instead of saying "put an aptamer at x = -27 nm, y = -15 nm,"
the mapped design says "put this aptamer on helix row 4, base position 4."

The current best mapped lead is:

```text
Best mapped layout: broad_grid_24_lattice
```

The average aptamer moved only about 1 nm during mapping. That means the
`broad_grid_24` pattern was already close to realistic origami attachment
positions.

### Step 5.13: Score Aptamer Orientation

```bash
python3 score_lattice_orientation.py
```

This checks whether each mapped aptamer site probably points upward toward the
EV or sideways/downward into a less useful direction.

Simple example: two aptamers can sit at almost the same x/y position, but the
one pointing upward is more useful than the one pointing sideways.

### Step 5.14: Optimize the Orientation Register

```bash
python3 optimize_lattice_orientation.py
```

This moves each aptamer to a nearby lattice site when that small move makes the
aptamer point upward.

The current best orientation-aware lead is:

```text
Best orientation-optimized layout: broad_grid_24_orientation_optimized
```

This design keeps the broad-grid shape, moves aptamers by about 2.5 nm on
average, and places all 24 aptamers in upward-facing sites under the simplified
orientation model.

### Step 5.15: Run the Simplified 3D Capture Screen

```bash
python3 simulate_3d_oriented_capture.py
```

This gives each aptamer a 3D base position and a 3D direction vector. It then
models each EV as a 73 nm sphere with sparse CD133 receptors on the lower
surface.

Simple example: an aptamer only gets a strong contact if the receptor is close
enough and inside the aptamer's 3D reach cone.

The current best 3D lead is:

```text
Best 3D layout: broad_grid_24_orientation_optimized
```

This is still a coarse 3D model, not full molecular dynamics. It is useful for
screening whether the layout survives a more realistic geometry check.

## 6. Optional Analyses

Test large, aggressive aptamer nets:

```bash
python3 high_capture_pattern_search.py
```

Test repeated tiles across a capture surface:

```bash
python3 simulate_multitile_capture_field.py
```

Apply background, activity-loss, and sample-damage penalties:

```bash
python3 evaluate_real_world_tradeoffs.py
```

These are supporting analyses. The aggressive 80-aptamer net is not the
recommended fabrication lead because crowding may make it physically
unrealistic.

## 7. Important Output Files

| File | Simple meaning |
|---|---|
| `robustness_layout_summary.csv` | Which original layout survives uncertain conditions? |
| `advanced_capture_design_summary.csv` | Which layout and surface idea work best together? |
| `inverse_design_requirements.csv` | What must improve to reach reliable 90% useful capture? |
| `capture_wash_protocol_best.csv` | Best modeled wash protocol for each design |
| `origami_buildability_scores.csv` | Which capture layout is also easiest to build on a DNA-origami tile? |
| `origami_lattice_mapping_comparison.csv` | How much each layout changes after snapping to origami sites |
| `origami_lattice_mapped_layouts.csv` | Helix/base attachment positions for the mapped layouts |
| `lattice_orientation_summary.csv` | How much aptamer direction changes capture after lattice mapping |
| `orientation_optimized_summary.csv` | Best nearby upward-facing lattice sites for each mapped design |
| `capture_3d_oriented_scores.csv` | Simplified 3D capture scores for orientation-optimized layouts |
| `advanced_capture_reliability.png` | Graph comparing advanced designs |
| `capture_wash_protocol_comparison.png` | Graph comparing optimized wash protocols |
| `origami_buildability_scores.png` | Graph comparing capture score against origami buildability |
| `origami_buildability_layouts.png` | Diagram of the top origami-friendly aptamer layouts |
| `origami_lattice_mapping.png` | Diagram comparing free layouts against lattice-snapped layouts |
| `lattice_orientation_scores.png` | Graph showing the orientation penalty |
| `orientation_optimized_layouts.png` | Diagram of the orientation-optimized aptamer sites |
| `capture_3d_oriented_scores.png` | Graph of the simplified 3D capture screen |
| `capture_3d_oriented_layouts.png` | Diagram of the 3D-oriented layout directions |
| `Research_Report_2026-06-20_Annotated.md` | Research report with figure-placement notes |

CSV means **comma-separated values**. It is a table that can be opened in
Excel, Google Sheets, Numbers, or a text editor.

JSON is a structured text format used to store detailed settings and summaries.

PNG is an image file containing a graph or layout diagram.

## 8. Linker Calibration Warning

The repository contains `linker_reach_models.csv`, which was generated using
the available oxDNA/oxRNA trajectories during development. The GitHub export
does not include those large trajectories or the full oxDNA source repository.

Do not rerun this command unless you have generated the required trajectory
folders:

```bash
python3 calibrate_linker_reach.py
```

Without those trajectories, the script uses a simpler mathematical fallback
and may overwrite the tracked calibrated table with different values.

## 9. Update Your Local Copy

To download future GitHub changes:

```bash
git pull origin main
```

To see files you changed locally:

```bash
git status
```

Generated `.npz` files, Python caches, and large trajectory tables are ignored
and should not be committed.

## 10. Scientific Limits

- The model compares design hypotheses; it does not validate a medical device.
- A simulated 90% score is not equivalent to 90% experimental recovery.
- Antifouling, spacer, damage, and wash values include explicit assumptions.
- Patient aqueous-humor samples are outside the current validation stage.
- Synthetic-target controls are required before clinical interpretation.
- Model assumptions should be replaced with measured values whenever possible.
