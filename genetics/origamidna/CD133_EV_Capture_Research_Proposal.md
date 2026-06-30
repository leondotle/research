# DNA-Origami-Scaffolded Aptamer Array Design for CD133+ Extracellular Vesicle Capture

## Project Summary

Extracellular vesicles (EVs) are nanoscale membrane-bound particles that carry molecular information from their cells of origin. CD133-positive EVs are of interest because CD133 is associated with stem-like and tumor-initiating cell states in several cancer contexts. However, EV capture remains technically challenging because target vesicles are heterogeneous, receptor copy number is finite, and multivalent binding depends strongly on nanoscale geometry.

This project proposes a computational design framework for DNA-origami-scaffolded CD133 aptamer arrays for EV capture. The recognition element is the Shigdar A15 CD133 aptamer, 5'-CCCUCCUACAUAGGG-3', represented in oxDNA/oxRNA with T used as U. The scaffold is modeled as a 90 nm x 60 nm addressable DNA-origami tile carrying defined aptamer anchor sites. Instead of treating aptamer placement as a qualitative design choice, the framework combines aptamer conformation modeling, poly-T linker reach calibration, and finite-receptor EV geometry scoring.

The current design stage evaluates four aptamer layouts: sparse_6, ring_12, grid_18, and dense_24. Each layout uses polyT20 linkers coupled to A15 aptamer heads. oxDNA trajectories are used to calibrate poly-T linker reach, restrained oxRNA simulations are used to estimate A15 head reach and stem stability, and a stochastic EV scoring model estimates multivalent capture probabilities for 50, 100, and 150 nm EVs across low, medium, and high CD133 surface densities.

At this stage, the project is not a validated EV isolation device or measured affinity assay. It is a computational and design-stage platform for identifying DNA-origami aptamer geometries that are worth building and testing experimentally.

[Figure 1 placeholder: insert overview schematic of the DNA-origami CD133 aptamer EV-capture concept.]

## Aim 1: Build an oxDNA/oxRNA-Guided A15 Aptamer and Linker Calibration Workflow

We will use oxDNA and oxRNA simulations to characterize the molecular reach of the CD133 aptamer capture construct. The A15 aptamer will be modeled with predicted stem restraints between base pairs 1-15, 2-14, 3-13, and 4-12. The simulated aptamer trajectory will be analyzed to estimate stem closure frequency, pairwise stem distances, and binding-face reach from the conjugation point.

In parallel, poly-T linkers of different lengths will be simulated to estimate linker end-to-end reach distributions. The current calibration includes polyT10, polyT15, polyT20, and polyT30. Linker reach and A15 head reach are combined by isotropic-orientation vector composition using Monte Carlo sampling, producing construct-level reach survival curves.

Preliminary work shows that the restrained A15 trajectory keeps all four predicted stem pairs closed in 63.4% of sampled frames under a 1.5 nm distance criterion. The A15 head reach has mean 2.17 nm, median 2.16 nm, and p90 2.72 nm. For the current polyT20/A15 construct, the combined reach distribution has mean 5.51 nm, median 5.55 nm, p90 7.65 nm, and approximately 6.28% probability of reaching 8 nm.

[Figure 2 placeholder: insert A15 restrained aptamer conformation or stem-distance summary.]

![Figure 3: calibrated linker-plus-aptamer reach survival curves](linker_reach_survival.png)

## Aim 1 Deliverables

- oxRNA-derived A15 stem-distance and head-reach distributions
- oxDNA-derived polyT10, polyT15, polyT20, and polyT30 reach models
- Combined linker-plus-aptamer reach survival curves
- Recommended linker lengths for scaffolded CD133 aptamer display

## Aim 2: Design DNA-Origami Aptamer Layouts for Multivalent CD133+ EV Capture

We will design addressable DNA-origami aptamer layouts that vary aptamer number, spacing, and spatial organization. The first-generation scaffold abstraction is a 90 nm x 60 nm tile with aptamer anchors placed at defined coordinates. The current layouts are sparse_6, ring_12, grid_18, and dense_24.

The scaffold is not intended to act only as a passive attachment surface. Its role is to control the nanoscale geometry of multivalent binding. By placing aptamers at known positions, the design can test whether dense arrays, grid-like arrays, or ring-like arrays better match the receptor-accessible surface of EVs with different diameters and CD133 densities.

Preliminary scoring favors the dense_24 design across the tested EV sizes and densities. Under the finite-receptor model, the best overall current case is dense_24 against a 150 nm high-CD133 EV, with mean contacts 4.334, P(at least 3 contacts) = 0.8121, P(at least 6 contacts) = 0.2681, and capture score 0.5676. At medium CD133 density, dense_24 gives P(at least 6 contacts) = 0.0001, 0.0035, and 0.0328 for 50, 100, and 150 nm EVs, respectively, showing that strong multivalent capture remains receptor-limited under the current assumptions.

![Figure 4: candidate DNA-origami aptamer layouts](ev_origami_aptamer_layouts.png)

![Figure 5: EV capture score heatmap across layout, EV size, and CD133 density](ev_capture_score_heatmap.png)

![Figure 6: strong multivalent-capture probability for medium CD133 density](ev_capture_multivalent_probability.png)

## Aim 2 Deliverables

- Scaffold layout table with aptamer coordinates and linker identities
- Schematic scaffold map for sparse_6, ring_12, grid_18, and dense_24
- EV capture score table across EV size and CD133 density
- Prioritized first-generation design for experimental fabrication

## Aim 3: Add Matched Controls and Extend the Model Toward Dynamic EV Capture Simulations

We will extend the current scoring workflow to include stronger negative and geometry controls. The first control class will randomize aptamer anchor positions on the same 90 nm x 60 nm tile while preserving aptamer count, polyT20/A15 reach, and finite-receptor EV scoring. These controls test whether designed layouts outperform matched random placement.

The second control class will represent scrambled or nonbinding aptamer chemistry. These controls are biologically distinct from random layout controls: random layouts test spatial organization, while scrambled aptamers test target-specific recognition. The scoring framework will be updated so nonbinding controls can be propagated through the same EV geometry model with reduced or absent receptor-specific binding probability.

Finally, the best designs will be exported toward coarse-grained dynamic simulation in LAMMPS or HOOMD-blue. This phase will move beyond static overlap/reach scoring by modeling EV approach, contact formation, receptor occupancy, and multivalent binding kinetics over time.

[Figure 7 placeholder: insert matched random-layout control comparison.]

[Figure 8 placeholder: insert dynamic EV capture simulation snapshot or trajectory summary.]

## Aim 3 Deliverables

- Matched random-layout control scores
- Scrambled/nonbinding aptamer control model
- Updated score comparisons for designed versus random and nonbinding controls
- Export-ready geometry files for LAMMPS or HOOMD-blue dynamic EV capture simulations
- Validation plan for no-EV, CD133-negative EV, CD133-positive EV, and mixed-EV conditions

## Significance

This research is valuable because it links molecular aptamer structure, linker mechanics, scaffold geometry, and EV-scale receptor statistics in one computational workflow. CD133+ EV capture depends on more than whether an aptamer can bind CD133 in isolation. It also depends on whether the aptamer can physically reach finite receptor sites on a curved vesicle surface, whether multiple receptors can be engaged simultaneously, and whether scaffold spacing helps or hurts multivalent contact formation.

DNA origami provides a route to control this geometry directly. A successful scaffolded aptamer array could support more selective EV capture, single-particle analysis, and multiplexed nanoscale presentation of capture ligands and controls. The framework can also reduce experimental trial-and-error by prioritizing geometries that are physically plausible before synthesis.

## Innovation

The project is innovative in three main ways:

- Aptamer-plus-linker reach calibration: poly-T linker trajectories and A15 aptamer head reach are combined into construct-level reach models rather than represented by a single arbitrary tether length.
- Finite-receptor EV scoring: CD133 sites are sampled as finite stochastic receptor positions on the lower EV hemisphere, and one receptor cannot contribute multiple simultaneous aptamer contacts.
- DNA-origami geometry optimization: aptamer count and anchor placement are treated as design variables that can be compared against matched random layouts and nonbinding chemistry controls.

## Preliminary Scaffold Design

The first scaffold design includes the following candidate layouts:

| Layout | Aptamer Count | Linker/Aptamer Construct | Geometry | Current Interpretation |
|---|---:|---|---|---|
| sparse_6 | 6 | polyT20/A15 | Six anchors on a broad ring | Weak under the current finite-receptor occupancy model |
| ring_12 | 12 | polyT20/A15 | Twelve-anchor circular layout | Improves coverage but remains weaker than grid_18 and dense_24 |
| grid_18 | 18 | polyT20/A15 | Rectangular grid across the tile | Current runner-up across most tested conditions |
| dense_24 | 24 | polyT20/A15 | Dense rectangular array | Current top design across tested EV sizes and CD133 densities |

The current best-scoring cases are:

| Layout | EV Diameter | CD133 Density | Mean Contacts | P>=3 Contacts | P>=6 Contacts | Capture Score |
|---|---:|---|---:|---:|---:|---:|
| dense_24 | 50 nm | high | 1.735 | 0.2385 | 0.0014 | 0.1512 |
| dense_24 | 100 nm | high | 3.198 | 0.6305 | 0.0968 | 0.3976 |
| dense_24 | 150 nm | high | 4.334 | 0.8121 | 0.2681 | 0.5676 |
| grid_18 | 150 nm | high | 2.584 | 0.5010 | 0.0213 | 0.2975 |

This layout set provides a compact first test of aptamer density, spatial organization, linker reach, and finite-receptor multivalent capture.

## Expected Outcomes

The expected outcome is a computationally grounded design package for DNA-origami-scaffolded CD133+ EV capture. The near-term deliverable is not a finished EV isolation product, but a reproducible design specification that identifies which aptamer layouts and linker constructs are most promising for experimental fabrication.

Longer term, the framework can support dynamic EV capture simulations, scrambled/nonbinding aptamer controls, mixed-EV selectivity tests, and broader DNA-origami capture platforms for EV subpopulation analysis.
