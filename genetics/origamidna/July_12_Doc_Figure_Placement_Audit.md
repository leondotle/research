# July 12 Research Update: Figure And Caption Audit

Source document reviewed:

`/Users/leonle/Downloads/July 12 Research Update.docx`

Main issue:

The Word document is still mostly the older June 20 report. It still says the recommended design is `clinical_grid_18`, while the current project lead is now `broad_grid_24_orientation_optimized`.

So several old placeholders/captions should be removed or replaced. Do not paste all old figures into the document unchanged.

## Highest-Priority Text Fixes Before Figures

Change these before adding images:

1. Title/date/current design block

   Current text says:

   `Research Report: Clinical Sparse-EV Model for DNA-Origami CD133 Capture`

   `Date: June 20, 2026`

   `Current recommended design direction: clinical_grid_18 repeated-exposure capture...`

   Change to:

   `Research Progress Report: Origami-Aware 3D Model for DNA-Origami CD133 Capture`

   `Date: July 12, 2026`

   `Current recommended design direction: broad_grid_24_orientation_optimized`

2. Main conclusion section

   Current text says the main conclusion is `clinical_grid_18`.

   Change the conclusion to:

   `The current lead is broad_grid_24_orientation_optimized because it combines strong modeled capture, better DNA-origami buildability, lattice-compatible attachment sites, upward-facing aptamer orientation, and the best simplified 3D sparse-EV performance.`

## Placeholder-By-Placeholder Actions

### Paragraphs 31-32

Current placeholder:

`[put candidate_origami_layouts_diagram here]`

Current caption:

`Caption: Original first-generation aptamer layouts on a 90 nm x 60 nm DNA origami tile. Each dot is one aptamer anchor.`

Recommended action:

Remove this old placeholder/caption, or replace it with the current buildability layout figure.

Use this image if you keep the figure:

`/Users/leonle/dna-origami/report_figures_2026-07-12/10_origami_buildability_layouts_optional.png`

New caption:

`Caption: Current candidate aptamer layouts after origami buildability scoring. The broad-grid layout became the practical lead because it preserves strong modeled capture while reducing crowding.`

Why:

The old caption is about first-generation layouts. The current report should focus on the updated origami-aware design.

### Paragraph 43

Current caption:

`Caption: Modeled clinical EV population: small EVs centered near 73 nm with sparse CD133 counts and mixed receptor patterns.`

Recommended action:

Keep this caption, but insert the missing figure directly before it.

Use this image:

`/Users/leonle/dna-origami/report_figures_2026-07-12/01_ev_population_summary.png`

Updated caption:

`Caption: Modeled clinical EV population centered near 73 nm with sparse CD133 counts and mixed receptor distributions.`

Why:

This figure is still relevant. It explains the clinical sparse-EV model.

### Paragraphs 53-54

Current placeholder:

`[put population_layout_performance_chart here]`

Current caption:

`Caption: Population-level layout scores across 160 modeled clinical EVs. Higher bars mean better average performance across the mixed EV population.`

Recommended action:

Replace this old population-optimizer figure with the current origami buildability figure.

Use this image:

`/Users/leonle/dna-origami/report_figures_2026-07-12/04_origami_buildability_scores.png`

New caption:

`Caption: Capture performance compared with DNA-origami buildability. broad_grid_24 becomes the practical lead because it combines high modeled capture with better spacing and lower crowding.`

Why:

The current decision point is no longer just population layout score. The newer result is capture plus origami buildability.

### Paragraphs 55-56

Current placeholder:

`[put top_population_optimized_layouts_diagram here]`

Current caption:

`Caption: Top population-optimized aptamer layouts. Each dot is one aptamer anchor on the DNA origami tile.`

Recommended action:

Replace with the orientation-optimized layout figure.

Use this image:

`/Users/leonle/dna-origami/report_figures_2026-07-12/07_orientation_optimized_layouts.png`

New caption:

`Caption: Orientation-optimized aptamer sites. broad_grid_24_orientation_optimized keeps the broad-grid shape while placing all 24 aptamers in upward-facing lattice positions.`

Why:

This is much more current than the old population-optimized layout diagram.

### Paragraphs 65-66

Current placeholder:

`[put moving_ev_contact_traces_chart here]`

Current caption:

`Caption: Moving-EV contact traces under sparse 73 nm clinical conditions. This shows why snapshot reach is not enough; the EV must remain captured while it moves.`

Recommended action:

Remove this placeholder and caption if the report is being updated to the July 12 version.

Optional replacement:

Use the 3D oriented capture score figure instead:

`/Users/leonle/dna-origami/report_figures_2026-07-12/08_3d_oriented_capture_scores.png`

New caption if replaced:

`Caption: Simplified 3D oriented capture screen. broad_grid_24_orientation_optimized performs best across 2, 5, and 10 CD133 receptor cases.`

Why:

The moving-EV contact trace is older. The newer 3D reach-cone model is more relevant to the current project state.

### Paragraphs 83-84

Current placeholder:

`[put aggressive_high_valency_search_chart here]`

Current caption:

`Caption: Aggressive high-valency sticky-net search. The dashed 90% line shows the target; even the best single-tile net does not reach it dynamically.`

Recommended action:

Remove this placeholder and caption unless you keep a short historical section explaining why the 80-aptamer net was rejected.

If you keep the historical section, use:

`/Users/leonle/dna-origami/high_capture_layout_scores.png`

But for the July 12 current report, removal is better.

Why:

The high-valency net is no longer central. The current lead is not a sticky net; it is an origami-aware broad grid.

### Paragraphs 103-104

Current placeholder:

`[put tile_spacing_comparison_diagram here]`

Current caption:

`Caption: Tile spacing comparison shows that 70 nm spacing often performs better than wider spacing because it reduces empty gaps.`

Recommended action:

Remove this placeholder and caption unless you keep the older repeated-tile field section.

If you keep the older repeated-field section, use:

`/Users/leonle/dna-origami/multitile_capture_field_scores.png`

But for the July 12 current report, removal is better.

Why:

The July 12 update focuses on single-tile origami realism, lattice placement, orientation, and 3D capture. Tile-spacing fields are older supporting work.

## Recommended Current Figure Order

If you update the document to the current July 12 story, use this order:

1. EV population model

   Image:

   `/Users/leonle/dna-origami/report_figures_2026-07-12/01_ev_population_summary.png`

   Caption:

   `Caption: Modeled clinical EV population centered near 73 nm with sparse CD133 counts and mixed receptor distributions.`

2. Advanced design screen

   Image:

   `/Users/leonle/dna-origami/report_figures_2026-07-12/02_advanced_design_screen.png`

   Caption:

   `Caption: Advanced design screen comparing layout and surface/spacer combinations under difficult simulated conditions.`

3. Wash protocol optimization

   Image:

   `/Users/leonle/dna-origami/report_figures_2026-07-12/03_wash_protocol_optimization.png`

   Caption:

   `Caption: Optimized capture-and-wash protocols. Several designs cross the simulated 90% useful-capture line under modeled assumptions.`

4. Origami buildability

   Image:

   `/Users/leonle/dna-origami/report_figures_2026-07-12/04_origami_buildability_scores.png`

   Caption:

   `Caption: Capture performance compared with DNA-origami buildability. broad_grid_24 becomes the practical lead because it combines high modeled capture with better spacing and lower crowding.`

5. Origami lattice mapping

   Image:

   `/Users/leonle/dna-origami/report_figures_2026-07-12/05_origami_lattice_mapping.png`

   Caption:

   `Caption: Free-coordinate layouts compared with lattice-snapped DNA-origami attachment sites. broad_grid_24 maps cleanly with only small coordinate shifts.`

6. Orientation penalty

   Image:

   `/Users/leonle/dna-origami/report_figures_2026-07-12/06_lattice_orientation_penalty.png`

   Caption:

   `Caption: Orientation penalty for nearest lattice-mapped layouts. The first broad-grid lattice map needed register optimization because nearest sites did not necessarily point upward.`

7. Orientation-optimized layout

   Image:

   `/Users/leonle/dna-origami/report_figures_2026-07-12/07_orientation_optimized_layouts.png`

   Caption:

   `Caption: Orientation-optimized aptamer sites. broad_grid_24_orientation_optimized places all 24 aptamers in upward-facing positions while preserving the broad-grid shape.`

8. Simplified 3D capture scores

   Image:

   `/Users/leonle/dna-origami/report_figures_2026-07-12/08_3d_oriented_capture_scores.png`

   Caption:

   `Caption: Simplified 3D oriented capture screen. broad_grid_24_orientation_optimized performs best across 2, 5, and 10 CD133 receptor cases.`

9. 3D-oriented layout directions

   Image:

   `/Users/leonle/dna-origami/report_figures_2026-07-12/09_3d_oriented_layouts.png`

   Caption:

   `Caption: Orientation-optimized layouts with 3D direction cues. This figure shows how the simplified 3D model treats aptamer display direction.`

## Short Version: What To Remove

Remove or replace these old captions:

- Original first-generation aptamer layouts.
- Population-level layout scores across 160 modeled EVs.
- Top population-optimized aptamer layouts.
- Moving-EV contact traces.
- Aggressive high-valency sticky-net search.
- Tile spacing comparison.

Keep/update this one:

- Modeled clinical EV population.

Add these current captions:

- Origami buildability.
- Origami lattice mapping.
- Orientation penalty.
- Orientation-optimized layout.
- Simplified 3D oriented capture.
- 3D-oriented layout directions.
