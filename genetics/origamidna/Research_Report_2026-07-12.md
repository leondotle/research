# Research Progress Report: Origami-Aware 3D Model for DNA-Origami CD133 Capture

Date: July 12, 2026

Project: DNA-origami-scaffolded A15 aptamer display for CD133+ extracellular vesicle capture

Current recommended design direction: `broad_grid_24_orientation_optimized`

## Executive Summary

This project is designing a DNA-origami aptamer surface to capture small, sparse CD133+ extracellular vesicles, or EVs.

The model originally focused mostly on whether aptamers could reach CD133 receptors in a simplified 2D/geometry screen. After feedback about real clinical EV size and receptor count, the project was rebuilt around more difficult clinical assumptions:

- EV diameter near 73 nm, not 150 nm.
- Sparse CD133 expression, often only 2, 5, or 10 receptors per EV.
- Mixed receptor patterns on the EV surface.
- Real-world penalties such as unwanted sticking, sample loss, and washing.
- DNA-origami buildability.
- Lattice-based attachment sites.
- Aptamer orientation.
- A simplified 3D reach-cone capture model.

The current main conclusion is:

`broad_grid_24_orientation_optimized` is the strongest current computational lead because it keeps a broad, evenly spaced aptamer pattern while placing all 24 aptamers in upward-facing lattice sites under the simplified orientation model.

In beginner terms: the best design is no longer just "many sticky hooks." It is now "24 hooks spread out evenly, placed on realistic origami sites, and adjusted so the hooks point upward toward the EV."

This is a stronger computational model than before, but it is still not a full experimental prediction. The current model is good for choosing which design to test next. It is not yet enough to guarantee a real laboratory capture percentage.

## Simple Definitions

EV: A tiny bubble-like particle released by cells. In this project, the EV is the target we want to catch.

CD133: A protein marker on some EVs. In the model, CD133 acts like a tiny handle on the EV surface.

Aptamer: A short DNA or RNA strand that binds a target. Here, the aptamer is the hook that tries to grab CD133.

DNA origami tile: A folded DNA platform. Think of it as a nanoscale pegboard where aptamer hooks can be placed.

Layout: The pattern of aptamer hooks on the DNA-origami tile.

Buildability: A score for whether the aptamer pattern is physically reasonable on a DNA-origami tile. For example, aptamers should not be too crowded.

Lattice site: A possible DNA-origami attachment point, described by a helix row and base position.

Orientation: The direction an aptamer points after attachment. An upward-facing aptamer is more useful than one pointing sideways or into the tile.

3D reach cone: The 3D space an aptamer can explore. A receptor must be close enough and in the right direction to count as a strong possible contact.

Useful capture: Target capture after subtracting practical problems such as unwanted binding, sample loss, or washing damage.

## Why The Project Changed

The earlier project direction favored dense or ring-like 24-aptamer layouts because those layouts performed well in simplified capture models.

However, the clinical target is very difficult:

- The EV is small, around 73-74 nm.
- CD133 receptors are sparse.
- A real DNA-origami tile cannot place aptamers at arbitrary perfect coordinates.
- Aptamers may point in different directions depending on the helix/base attachment site.

In simple language: it is not enough to ask, "Where should the hooks go?" We also need to ask, "Can the hooks actually be built there, and are they pointing toward the EV?"

The project now adds those origami-specific checks.

## Current Computational Workflow

### Step 1: Generate A Realistic Sparse EV Population

Output files:

- `ev_population_generator.py`
- `ev_population_clinical.csv`
- `ev_population_clinical_summary.json`

The model creates 160 synthetic EVs with realistic variation.

| Property | Value |
|---|---:|
| Number of EVs | 160 |
| Mean diameter | 73.3 nm |
| Median diameter | 73.4 nm |
| Diameter p10-p90 | 59.6-88.9 nm |
| Mean CD133 count | 6.4 |
| Median CD133 count | 6 |
| Receptor count p10-p90 | 2.0-11.1 |

The model also varies receptor patterns.

| CD133 Pattern | Count |
|---|---:|
| Random | 65 |
| Single cluster | 43 |
| Two cluster | 32 |
| Bottom cap | 20 |

Beginner example: instead of testing one perfect EV, the model creates many tiny balls with different sizes and different CD133 sticker patterns.

[put `ev_population_clinical_summary.json` table or population figure here]

Caption: Modeled clinical EV population centered near 73 nm with sparse CD133 counts and mixed receptor distributions.

### Step 2: Search Layouts And Surface Formulations

Output files:

- `advanced_capture_design_search.py`
- `advanced_capture_design_summary.csv`
- `advanced_capture_reliability.png`

This step compares layouts and surface/spacer assumptions. The goal is to test whether layout shape plus material treatment can improve useful capture.

Top advanced-screen results before wash optimization:

| Rank | Layout | Formulation | Aptamers | p10 Useful Capture |
|---:|---|---|---:|---:|
| 1 | `triple_ring_24` | `mixed_spacer_antifouling` | 24 | 0.7791 |
| 2 | `broad_grid_24` | `mixed_spacer_antifouling` | 24 | 0.7474 |
| 3 | `broad_grid_20` | `rigid_spacer_PEG` | 20 | 0.7415 |
| 4 | `sunflower_20` | `mixed_spacer_antifouling` | 20 | 0.7344 |

At this stage, `triple_ring_24` still looked strongest by capture alone.

Beginner explanation: this is like testing several hook patterns and surface coatings to see which one catches best under difficult simulated conditions.

[put `advanced_capture_reliability.png` here]

Caption: Advanced design screen comparing layout and surface/spacer combinations under difficult simulated conditions.

### Step 3: Optimize Capture And Washing

Output files:

- `optimize_capture_wash_protocol.py`
- `capture_wash_protocol_best.csv`
- `capture_wash_protocol_comparison.png`

This step asks how much passivation, contact improvement, and washing are needed to keep useful capture near or above 90%.

Several layouts can reach a modeled 90% useful-capture threshold after protocol optimization.

| Layout | Formulation | p10 Useful Capture | Reliable 90? |
|---|---|---:|---|
| `triple_ring_24` | `mixed_spacer_antifouling` | 0.9005 | yes |
| `broad_grid_24` | `rigid_spacer_PEG` | 0.9009 | yes |
| `hybrid_24` | `rigid_spacer_PEG` | 0.9006 | yes |
| `broad_grid_20` | `rigid_spacer_PEG` | 0.9027 | yes |
| `broad_grid_24` | `mixed_spacer_antifouling` | 0.9046 | yes |

This is important, but it is not the final answer because capture score alone does not prove the pattern is easy to build on DNA origami.

Beginner explanation: washing can remove unwanted weak sticking, but if the wash is too strong it may also remove target EVs. This step searches for a balance.

[put `capture_wash_protocol_comparison.png` here]

Caption: Optimized capture-and-wash protocols. Several designs cross the simulated 90% useful-capture line under modeled assumptions.

### Step 4: Add DNA-Origami Buildability

Output files:

- `score_origami_buildability.py`
- `origami_buildability_scores.csv`
- `origami_buildability_scores.png`
- `origami_buildability_layouts.png`

This step asks whether a layout is practical as a DNA-origami pattern.

The buildability score checks:

- aptamer spacing;
- crowding;
- edge margin;
- total aptamer count.

Best buildability-adjusted results:

| Rank | Layout | Capture Score Used | Origami Buildability | Practical Score |
|---:|---|---:|---:|---:|
| 1 | `broad_grid_24` | 0.9046 | 1.0000 | 0.9313 |
| 2 | `broad_grid_20` | 0.9027 | 1.0000 | 0.9299 |
| 3 | `hybrid_24` | 0.9006 | 0.8350 | 0.8822 |
| 6 | `triple_ring_24` | 0.9005 | 0.5478 | 0.8017 |

This changed the design direction.

`triple_ring_24` still catches well in the optimized model, but it has crowding issues. `broad_grid_24` also catches well and is much easier to defend as a buildable origami layout.

Beginner explanation: two hook patterns can both catch well, but the one with evenly spaced hooks is easier to build than one with hooks packed too close together.

[put `origami_buildability_scores.png` here]

Caption: Capture performance compared with DNA-origami buildability. `broad_grid_24` becomes the practical lead.

### Step 5: Snap The Layout To A DNA-Origami Lattice

Output files:

- `map_layout_to_origami_lattice.py`
- `origami_lattice_mapping_comparison.csv`
- `origami_lattice_mapped_layouts.csv`
- `origami_lattice_mapping.png`

This step changes a free x/y dot layout into simplified DNA-origami attachment sites.

Example:

Before:

`aptamer at x = -27 nm, y = -15 nm`

After:

`aptamer on helix row 4, base position 4`

Mapping results:

| Layout | Mapped Layout | Mean Snap Shift | Max Snap Shift | Buildability After Mapping | Sparse Score Change |
|---|---|---:|---:|---:|---:|
| `broad_grid_24` | `broad_grid_24_lattice` | 1.055 nm | 1.280 nm | 1.0000 | -0.0014 |
| `broad_grid_20` | `broad_grid_20_lattice` | 1.247 nm | 1.990 nm | 1.0000 | +0.0009 |
| `triple_ring_24` | `triple_ring_24_lattice` | 1.022 nm | 1.953 nm | 0.7441 | -0.0021 |

The average movement for `broad_grid_24` was only about 1 nm, meaning the original broad-grid pattern was already close to realistic origami lattice positions.

Beginner explanation: the dots did not have to move much to land on allowed pegboard holes.

[put `origami_lattice_mapping.png` here]

Caption: Free-coordinate layouts compared with lattice-snapped DNA-origami attachment sites.

### Step 6: Score Aptamer Orientation

Output files:

- `score_lattice_orientation.py`
- `lattice_orientation_summary.csv`
- `lattice_orientation_scores.png`
- `lattice_orientation_layouts.png`

This step asks whether each mapped aptamer points upward toward the EV.

Why this matters:

An aptamer can be present but still not useful if it points sideways or toward the tile.

Simple example:

- Upward-facing hook: likely useful.
- Sideways hook: partly useful.
- Downward or blocked hook: weak or mostly wasted.

The first nearest-site lattice mapping found a problem: `broad_grid_24_lattice` placed all 24 aptamers on poor-facing sites under the simplified orientation model.

This did not mean the broad grid was bad. It meant nearest-site mapping was too naive.

[put `lattice_orientation_scores.png` here]

Caption: Orientation penalty for nearest lattice-mapped layouts. The broad grid needed register optimization.

### Step 7: Optimize The Orientation Register

Output files:

- `optimize_lattice_orientation.py`
- `orientation_optimized_summary.csv`
- `orientation_optimized_mapped_layouts.csv`
- `orientation_optimized_layouts.png`

This step moves each aptamer to a nearby lattice site when that small move makes the aptamer point upward.

Best orientation-optimized results:

| Rank | Optimized Layout | Aptamers | Up-Facing | Mean Shift | Oriented Sparse Score |
|---:|---|---:|---:|---:|---:|
| 1 | `broad_grid_24_orientation_optimized` | 24 | 24/24 | 2.515 nm | 0.1750 |
| 2 | `triple_ring_24_orientation_optimized` | 24 | 16/24 | 2.409 nm | 0.1595 |
| 3 | `broad_grid_20_orientation_optimized` | 20 | 12/20 | 1.983 nm | 0.1544 |

The broad grid was rescued by small register shifts. It kept the broad-grid shape, but all 24 aptamers became upward-facing under the simplified model.

Beginner explanation: the design moved each hook from one pegboard hole to a nearby better-facing hole.

[put `orientation_optimized_layouts.png` here]

Caption: Orientation-optimized aptamer sites. `broad_grid_24_orientation_optimized` places all 24 aptamers in upward-facing positions.

### Step 8: Run A Simplified 3D Capture Screen

Output files:

- `simulate_3d_oriented_capture.py`
- `capture_3d_oriented_scores.csv`
- `capture_3d_oriented_summary.json`
- `capture_3d_oriented_scores.png`
- `capture_3d_oriented_layouts.png`

This is the newest realism upgrade.

The model now gives each aptamer:

- x position;
- y position;
- z height;
- 3D direction vector;
- 3D reach cone.

The EV is modeled as a 73 nm sphere with CD133 receptor points on the lower surface.

A receptor only counts strongly if it is close enough and inside the aptamer's 3D reach cone.

Best simplified 3D results:

| Layout | Up-Facing Aptamers | Mean 3D Sparse Score | Mean Chance of At Least 1 Contact |
|---|---:|---:|---:|
| `broad_grid_24_orientation_optimized` | 24/24 | 0.2707 | 0.3939 |
| `triple_ring_24_orientation_optimized` | 16/24 | 0.2426 | 0.3649 |
| `broad_grid_20_orientation_optimized` | 12/20 | 0.2170 | 0.3419 |

For the current lead, results by receptor count were:

| Receptors Per 73 nm EV | Mean Contacts | P(at least 1 contact) | P(at least 2 contacts) | P(at least 3 contacts) | 3D Sparse Score |
|---:|---:|---:|---:|---:|---:|
| 2 | 0.3042 | 0.2347 | 0.0695 | 0.0000 | 0.1249 |
| 5 | 0.7659 | 0.4003 | 0.1960 | 0.1059 | 0.2656 |
| 10 | 1.4777 | 0.5466 | 0.3644 | 0.2226 | 0.4217 |

The simplified 3D model still supports `broad_grid_24_orientation_optimized` as the current lead.

Beginner explanation: this step checks whether the hooks can reach the tiny handles in actual 3D space, not just on a flat drawing.

[put `capture_3d_oriented_scores.png` here]

Caption: Simplified 3D oriented capture screen. `broad_grid_24_orientation_optimized` performs best across 2, 5, and 10 CD133 receptor cases.

[put `capture_3d_oriented_layouts.png` here]

Caption: Orientation-optimized layouts with 3D direction cues.

## Current Best Design

The current computational lead is:

`broad_grid_24_orientation_optimized`

This means:

- `broad`: aptamers cover a wider contact area.
- `grid`: aptamers are arranged in rows and columns.
- `24`: one tile has 24 aptamer attachment points.
- `orientation_optimized`: aptamers were moved to nearby lattice sites so they point upward.

Why it is the best current lead:

- It reaches modeled 90% useful capture after protocol optimization.
- It has perfect buildability score under the current spacing/crowding rules.
- It maps cleanly onto simplified DNA-origami lattice sites.
- It can be adjusted so all 24 aptamers are upward-facing.
- It remains the best layout in the simplified 3D reach-cone screen.

## Accuracy And Limits

The current model is much more realistic than the original version, but it is still a screening model.

What the model now includes:

- realistic 73 nm EV size;
- sparse CD133 counts;
- mixed receptor distributions;
- repeated capture and washing logic;
- real-world penalty assumptions;
- DNA-origami spacing and crowding checks;
- helix/base-like lattice sites;
- aptamer orientation;
- simplified 3D reach cones.

What the model still does not include:

- real caDNAno scaffold routing;
- exact staple sequences;
- full oxDNA folding of the entire tile;
- explicit 3D folding of the A15 aptamer;
- measured CD133 binding rates for this exact construct;
- measured nonspecific binding in the real buffer;
- real microfluidic or bead-mixing flow geometry;
- experimental validation with synthetic targets.

Best accuracy statement:

The model is strong for ranking design hypotheses, but not yet strong enough to predict an exact wet-lab capture rate.

Simple version:

Good for choosing the next design to test. Not enough to guarantee the lab result.

## Recommended Next Step

The best next step is to convert `broad_grid_24_orientation_optimized` into a more explicit experimental design candidate.

Recommended next computational task:

Create a design handoff table that lists each aptamer attachment site as:

- anchor ID;
- x/y coordinate;
- helix row;
- base position;
- orientation class;
- spacer recommendation;
- aptamer modification note.

Recommended next experimental interpretation:

Use this as a synthetic-buffer validation candidate, not a clinical sample claim.

Beginner explanation: the model has chosen the best hook pattern. The next step is to turn that hook pattern into something a wet-lab person can actually review and eventually build.

## Final Takeaway

The project has progressed from a simple capture geometry model to an origami-aware, orientation-aware, simplified 3D screening workflow.

The current best design is:

`broad_grid_24_orientation_optimized`

The main scientific claim should be:

This design is the strongest current computational candidate because it combines high modeled capture, better DNA-origami buildability, realistic lattice placement, upward-facing aptamer orientation, and the best simplified 3D sparse-EV capture performance.

The claim should not be:

This design is proven to capture 90% of real clinical EVs.

The correct cautious wording is:

`broad_grid_24_orientation_optimized` is the best current design to advance toward synthetic-target validation.
