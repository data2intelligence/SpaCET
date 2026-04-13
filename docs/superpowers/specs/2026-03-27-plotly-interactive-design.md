# Design: Standalone Plotly Interactive Mode for SpaCET

## Summary

Add `interactive = "plotly"` option to `SpaCET.visualize.spatialFeature()` that returns a standalone Plotly htmlwidget. Coexists with the existing Shiny app (`interactive = TRUE`) and static ggplot (`interactive = FALSE`).

## Motivation

The spatial-gpu Python package provides lightweight inline Plotly figures via `interactive=True`. SpaCET's current interactive mode requires a full Shiny server. A standalone Plotly widget is more portable: works in RStudio viewer, R Markdown, Jupyter, and can be saved as self-contained HTML.

## Design

### Parameter change

`interactive` changes from logical-only to accepting `FALSE`, `TRUE`, or `"plotly"`:

- `FALSE` -- static ggplot (unchanged)
- `TRUE` -- Shiny app (unchanged)
- `"plotly"` -- standalone Plotly htmlwidget (new)

Backward-compatible: all existing calls work unchanged.

### Supported spatial types (10 of 11)

All types except CellTypeComposition (pie charts). Specifically:
QualityControl, GeneExpression, CellFraction, MostAbundantCellType,
LRNetworkScore, Interface, GeneSetScore, SecretedProteinActivity,
SignalingPattern, metaData.

CellTypeComposition raises a clear error directing users to other modes.

### Implementation

1. **Data preparation**: Reuse existing code (lines 60-349) that constructs `visiualVector` from `spatialType`. This logic is mode-independent.

2. **New helper `visualSpatialPlotly()`**: Renders a single spatial feature as a Plotly trace.
   - `plotly::plot_ly()` with `type = "scattergl"` (WebGL for large datasets)
   - Continuous color: `colorscale` from the `colors` vector
   - Discrete color: one trace per category via `plotly::add_trace()`
   - Image background: `layout(images = ...)` with base64-encoded raster
   - Hover tooltip: spot ID, feature name, value
   - Aspect ratio locked via `xaxis.scaleanchor = "y"`
   - Coordinate transforms match existing `visualSpatial()` logic

3. **Multi-feature panels**: `plotly::subplot()` with `nrows` parameter, matching the `nrow` parameter behavior of patchwork.

4. **Return value**: A plotly htmlwidget object. Users can display inline or save via `htmlwidgets::saveWidget()`.

### Dependencies

`plotly` already in Suggests. No new dependencies needed.
