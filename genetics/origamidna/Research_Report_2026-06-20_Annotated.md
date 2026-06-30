# Research Report: Clinical Sparse-EV Model for DNA-Origami CD133 Capture

**Date:** June 20, 2026  
**Project:** DNA-origami-scaffolded A15 aptamer display for CD133+ extracellular vesicle capture  
**Current recommended design direction:** `clinical_grid_18` repeated-exposure capture under clean-buffer synthetic validation  

## Executive Summary

This project originally tested DNA-origami aptamer layouts against optimistic extracellular vesicle, or EV, assumptions. The first model used large EVs around 150 nm and high CD133 receptor availability. Under those assumptions, a dense 24-aptamer layout with a polyT30 linker looked strongest.

After advisor feedback, the model was rebuilt around more realistic patient aqueous humor constraints:

- CD133+ EV diameter near 73-74 nm instead of 150 nm.
- Sparse CD133 availability, often only 1-10 CD133 molecules per EV.
- Different possible CD133 patterns on the EV surface, including random and clustered placement.
- Repeated near-surface capture attempts, as would happen in bead mixing, incubation, or microfluidic recirculation.
- Practical penalties for nonspecific binding, sample damage/loss, and aptamer/origami activity decay.

The main conclusion changed.

The project should no longer present a single dense 24-aptamer tile as the final lead. A more realistic conclusion is:

**For sparse 73 nm CD133+ EVs, moderate-density `clinical_grid_18` tiles used with repeated clean-buffer encounters are more defensible than a single overpacked high-valency sticky net.**

In beginner terms: one very sticky trap can catch more, but it may also catch junk. A cleaner moderate trap used several times can be more useful.

## Simple Definitions

**EV:** A tiny bubble-like particle released by cells. In this project, the EV is the thing we want to catch.

**CD133:** A protein marker on some EVs. In the model, CD133 acts like a tiny handle on the EV.

**Aptamer:** A short DNA or RNA strand that binds a target. Here, the A15 aptamer is the hook that tries to bind CD133.

**DNA origami tile:** A tiny folded DNA platform. Think of it as a nanoscale board where we can place aptamer hooks.

**Layout:** The pattern of aptamer hooks on the DNA origami tile.

**Capture probability:** The chance that an EV gets caught by at least one aptamer/CD133 contact for long enough to count.

**Useful capture:** Target capture after subtracting practical problems, such as catching unwanted particles or damaging EVs.

## Why The Model Had To Change

The earlier model favored `dense_24/polyT30` because it assumed relatively large EVs and high CD133 availability. That made strong multivalent contact easier.

Advisor feedback changed the core assumptions:

- Real clinical CD133+ EVs in aqueous humor are much smaller, near 73-74 nm.
- Real CD133 counts are sparse, around 8.5-10.8 molecules per EV on average, and possibly 1-3 in low-expression settings.
- Patient samples are not available for near-term validation; summer testing should use synthetic targets in standard buffers.

In simple language: the old model was testing a big ball with many handles. The real target is a tiny ball with only a few handles.

[put ev_origami_aptamer_layouts.png here]

Caption idea: Original first-generation aptamer layouts on a 90 nm x 60 nm DNA origami tile. Each dot is one aptamer anchor.

## Updated Computational Workflow

### Step 1: Generate A Realistic EV Population

The new model creates a population of synthetic EVs instead of testing only one ideal EV.

Output files:

- `ev_population_generator.py`
- `ev_population_clinical.csv`
- `ev_population_clinical_summary.json`

The generated population had:

| Property | Value |
|---|---:|
| Number of EVs | 160 |
| Mean diameter | 73.3 nm |
| Median diameter | 73.4 nm |
| Diameter p10-p90 | 59.6-88.9 nm |
| Mean CD133 count | 6.4 |
| Median CD133 count | 6 |
| Receptor count p10-p90 | 2.0-11.1 |

The model also varied how CD133 molecules sit on each EV:

| CD133 Pattern | Count |
|---|---:|
| Random | 65 |
| Single cluster | 43 |
| Two cluster | 32 |
| Bottom cap | 20 |

Beginner example: instead of testing one perfect ball, we made 160 small balls with different sizes and different sticker patterns.

[optional: put ev_population_clinical_summary.json values into a small schematic or table here]

Caption idea: Modeled clinical EV population: small EVs centered near 73 nm with sparse CD133 counts and mixed receptor patterns.

### Step 2: Compare Many Aptamer Layouts

The layout optimizer tested 334 aptamer layouts against the same EV population.

Output files:

- `optimize_population_layouts.py`
- `population_layout_scores.csv`
- `population_layout_optimization_summary.json`
- `population_optimized_layouts.png`

The best snapshot-layout results were:

| Rank | Layout | Aptamers | Population Score | Snapshot Capture |
|---:|---|---:|---:|---:|
| 1 | clinical_dense_24 | 24 | 0.4751 | 0.5400 |
| 2 | clinical_grid_18 | 18 | 0.4656 | 0.5213 |
| 3 | clinical_grid_12 | 12 | 0.4543 | 0.5156 |
| 4 | evolved grid_18 variant | 18 | 0.4505 | 0.5181 |

Snapshot capture means the model asks: if the EV is near the tile, can the aptamers reach CD133?

This is useful, but not enough. EVs move.

[put population_layout_top_scores.png here]

Caption idea: Population-level layout scores across 160 modeled clinical EVs. Higher bars mean better average performance across the mixed EV population.

[put population_optimized_layouts.png here]

Caption idea: Top population-optimized aptamer layouts. Each dot is one aptamer anchor on the DNA origami tile.

### Step 3: Brownian Motion Validation

The top layouts were then checked with moving EVs.

Output files:

- `validate_population_dynamics.py`
- `population_layout_dynamics_validation.csv`
- `population_layout_dynamics_validation_summary.json`

Best dynamic results:

| Rank | Layout | Aptamers | Dynamic Capture |
|---:|---|---:|---:|
| 1 | clinical_dense_24 | 24 | 0.2750 |
| 2 | clinical_grid_12 | 12 | 0.2562 |
| 3 | clinical_grid_18 | 18 | 0.2500 |

Dynamic capture is stricter. It asks: can the EV get caught while it wiggles around?

The best single-tile dynamic result was only 27.5%. This means a single tile is not enough to reliably capture sparse 73 nm EVs.

[optional: put ev_capture_clinical_73nm_dynamics_contacts.png here]

Caption idea: Moving-EV contact traces under sparse 73 nm clinical conditions. This shows why snapshot reach is not enough; the EV must remain captured while it moves.

## Testing A Very Sticky Net

To see whether capture could be pushed above 90%, an aggressive high-capture search tested larger aptamer nets with up to 80 aptamers.

Output files:

- `high_capture_pattern_search.py`
- `high_capture_layout_scores.csv`
- `high_capture_dynamic_validation.csv`
- `high_capture_best_layouts.csv`
- `high_capture_pattern_search_summary.json`

Best aggressive single-tile result:

| Test | Best Layout | Aptamers | Capture |
|---|---|---:|---:|
| Snapshot | random_net_80_01 | 80 | 0.5656 |
| Dynamic | random_net_80_02 | 80 | 0.4062 |

Even an 80-aptamer net did not reach 90% in dynamic validation.

Beginner explanation: making one trap very sticky helps, but if the EV has only a few handles and moves away, one trap still misses many EVs.

The 80-aptamer net is also physically risky:

- too many aptamers may crowd each other;
- DNA charge may repel membranes;
- fabrication may be harder;
- nonspecific binding may increase.

So the high-valency net is useful as an upper-bound simulation, not as the best experimental design.

[put high_capture_layout_scores.png here]

Caption idea: Aggressive high-valency sticky-net search. The dashed 90% line shows the target; even the best single-tile net does not reach it dynamically.

[optional: put high_capture_best_layouts.png here]

Caption idea: Best aggressive sticky-net layouts. These show what high aptamer density looks like, but they may be physically crowded.

## Repeated Tile Capture Field

The next model tested repeated moderate-density tiles on a surface.

Output files:

- `simulate_multitile_capture_field.py`
- `multitile_capture_field_summary.csv`
- `multitile_capture_field_summary.json`
- `multitile_capture_field_layouts.csv`
- `multitile_capture_field_scores.png`

The best repeated-field result was:

| Field | Unit Tile | Tile Copies | Aptamers Per Tile | Total Aptamers | Spacing | Field Size | Single-Encounter Capture |
|---|---|---:|---:|---:|---:|---|---:|
| clinical_grid_18_4x3_70nm | clinical_grid_18 | 12 | 18 | 216 | 70 nm | 300 x 165 nm | 0.3063 |

This did not reach 90% in one near-surface encounter.

However, repeated encounters help. If the same EV gets several near-surface chances, projected cumulative capture improves:

| Near-Surface Encounters | Projected Capture |
|---:|---:|
| 1 | 0.3063 |
| 3 | 0.6662 |
| 5 | 0.8394 |
| 6 | 0.8886 |
| 7 | 0.9227 |
| 8 | 0.9464 |
| 10 | 0.9742 |

Beginner example: one pass over Velcro may miss. Several passes give more chances.

[put multitile_capture_field_scores.png here]

Caption idea: Repeated-tile capture fields. The best single near-surface encounter remains below 90%, but repeated encounters can cross 90%.

[put multitile_capture_field_layouts.png here]

Caption idea: Top repeated-tile field layouts. Each cluster of dots is a repeated DNA origami tile unit.

## Does Spacing Matter?

Yes. Spacing mattered a lot.

The model showed:

- closer tile spacing, around 70 nm, often helped;
- wider spacing, around 90-110 nm, sometimes created gaps;
- simply making the field bigger did not always help;
- too large a field can spread the same EV-starting probability over more area and reduce effective contact.

In simple language: a bigger capture area is not automatically better. If the sticky pads are too far apart, EVs can drift through empty space.

[refer back to multitile_capture_field_scores.png here]

Caption idea: Tile spacing comparison shows that 70 nm spacing often performs better than wider spacing because it reduces empty gaps.

## Real-World Penalty Model

A practical model was added to ask whether high capture remains useful after real-world problems.

Output files:

- `evaluate_real_world_tradeoffs.py`
- `real_world_tradeoff_scores.csv`
- `real_world_tradeoff_summary.json`
- `real_world_tradeoff_scores.png`

The penalty model included:

- nonspecific binding, meaning unwanted particles stick;
- EV/sample damage or loss;
- aptamer/origami activity loss over repeated encounters;
- higher false-capture risk for denser aptamer fields.

Three scenarios were tested:

| Scenario | Meaning |
|---|---|
| clean_buffer | Synthetic targets in controlled buffer |
| moderate_background | Some unwanted proteins or particles |
| dirty_biofluid_like | Harsher biological-fluid-like background |

After penalties, the best useful design was not the large 12-tile field. It was a smaller repeated-exposure design:

**`clinical_grid_18_1x1_90nm` with repeated near-surface encounters.**

Results:

| Scenario | Encounters | Target Capture | False Capture Risk | Damage/Loss Risk | Useful Score |
|---|---:|---:|---:|---:|---:|
| clean_buffer | 10 | 0.9296 | 0.0124 | 0.0489 | 0.9142 |
| moderate_background | 10 | 0.9198 | 0.0288 | 0.0956 | 0.8818 |
| dirty_biofluid_like | 10 | 0.8986 | 0.0568 | 0.1829 | 0.8044 |

This is the key practical conclusion.

**A moderate tile used repeatedly can be more useful than a giant sticky field, because giant sticky fields also increase background capture risk.**

[put real_world_tradeoff_scores.png here]

Caption idea: Useful capture after real-world penalties. The best practical design is not the largest sticky field; it is a moderate `clinical_grid_18` tile used repeatedly.

## Current Recommended Design

The current recommended experimental design is:

**Primary layout:** `clinical_grid_18/polyT30`  
**Deployment mode:** repeated near-surface exposure in clean buffer  
**Validation target:** synthetic CD133+ EV-like particles or synthetic CD133-displaying vesicles  
**Exposure series:** 1x, 3x, 5x, 8x, 10x equivalent encounters  

The model suggests that `clinical_grid_18` is a good compromise:

- fewer aptamers than dense_24;
- lower crowding risk;
- good population-level performance;
- better useful-capture score after penalties;
- easier to explain and defend experimentally.

In simple terms: it is not the stickiest trap, but it may be the cleanest useful trap.

## Proposed Wet-Lab Validation Controls

The model should be connected to synthetic validation with controls:

| Condition | Purpose |
|---|---|
| clinical_grid_18 + CD133+ target | Main capture test |
| scrambled aptamer + CD133+ target | Tests nonspecific aptamer effects |
| clinical_grid_18 + CD133-negative EV-like target | Tests marker specificity |
| blank origami + CD133+ target | Tests origami/background sticking |
| no-origami surface | Tests surface-only binding |
| exposure time series | Tests whether repeated encounters improve capture |

Beginner explanation: controls tell us whether we caught the target for the right reason.

## Limitations

This model is still not a full physical or clinical predictor.

Important limitations:

- Penalty values are scenario assumptions, not measured constants.
- EVs are treated as simplified spheres.
- Receptor mobility on the EV membrane is still simplified.
- DNA origami deformation is not fully modeled.
- Electrostatic repulsion is represented only indirectly.
- Nonspecific binding is modeled as a penalty, not with molecular detail.
- Patient aqueous humor samples are not modeled directly and should remain future-stage.

This model is best described as:

**a clean-buffer synthetic-target design prioritization tool.**

It should not be described as:

**a validated patient-sample performance predictor.**

## Recommended Next Steps

### 1. Use `clinical_grid_18/polyT30` As The Main Synthetic Test Layout

This layout is the best current balance between capture and realism.

### 2. Run A Synthetic Exposure Series

Test something like:

- 1x exposure;
- 3x exposure;
- 5x exposure;
- 8x exposure;
- 10x exposure.

The goal is to see whether repeated encounters improve target capture without unacceptable background.

### 3. Measure Nonspecific Binding Experimentally

The current model assumes nonspecific binding rates. These should be replaced with measured values from controls.

### 4. Keep `clinical_dense_24` As A Backup

`clinical_dense_24` performed well in snapshot and single-tile dynamic tests, but may carry higher crowding/background risk.

### 5. Avoid Promoting The 80-Aptamer Net As The Main Experimental Design

It improves capture in the model but is likely less realistic and may increase nonspecific binding.

## Final Conclusion

The project has moved from an optimistic single-EV geometry screen to a more realistic clinical sparse-EV modeling framework.

The best current story is:

**Sparse 73 nm CD133+ EVs are too difficult for one DNA-origami tile to capture reliably. Very dense sticky nets improve capture but are physically and experimentally risky. A more defensible strategy is repeated clean-buffer exposure to moderate-density `clinical_grid_18/polyT30` origami tiles, with proper scrambled, CD133-negative, blank-origami, and surface-only controls.**

This is a solid Aim 1 result and a much stronger foundation for advisor discussion than the original dense_24/150 nm model.
