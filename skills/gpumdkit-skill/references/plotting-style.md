# GPUMDkit Plotting Style

Read this reference when creating a plotting script, changing the appearance of an existing plot, or reviewing a plot contribution. Also read `visualization.md`, `contributing.md`, and the reference governing the scientific quantity. This guide distills recurring choices in `Scripts/plt_scripts`; it does not replace script-specific behavior.

## Design Character

GPUMDkit figures are compact, publication-oriented scientific plots. They favor a white background, restrained Matplotlib colors, direct physical labels with units, and enough annotation to interpret a result without decorative UI elements. The data and its uncertainty are visually dominant.

There are two established families, and both are valid:

- General analysis and NEP plots use Matplotlib's object-oriented API, a sans-serif font stack, compact multi-panel layouts, and `tight_layout()`.
- Thermal-transport and phonon plots retain the older `pylab` style, inward ticks on all four sides, minor ticks, thicker final/average curves, and frameless legends.

Do not force an existing script from one family into the other. Match neighboring scripts in the same scientific workflow.

## Visual Vocabulary

### Typography

- Prefer `sans-serif` with `Arial`, `DejaVu Sans`, and `Liberation Sans` as fallbacks.
- Use math text and explicit units in axis labels, for example `Å$^2$/ps`, `eV/atom`, and `W/mK`.
- Keep label and tick sizes internally consistent within a figure. Existing figures commonly use roughly 10-13 pt, but this is an observation, not a mandatory global value.
- Use concise legends and annotations. Put fit statistics or derived values close to the relevant panel or curve.

### Color and line hierarchy

- Use Matplotlib cycle colors (`C0`, `C1`, `C2`, ...) for related components. Directional data conventionally maps `x`, `y`, and `z` to `C0`, `C1`, and `C2` when the target script already follows that mapping.
- Use grey or black for references, raw replicas, and identity/zero lines. A grey dashed line is the common parity reference.
- Show noisy replicas or distributions with low alpha; emphasize the mean, cumulative result, fit, or final trajectory with a thicker opaque line.
- Use translucent bands for uncertainty and open markers when fitted points must remain distinct from fit lines.
- Preserve a script's existing color semantics. Never change colors in a way that changes the meaning of components between related figures.

### Axes and layout

- Choose the layout from the scientific comparison: compact single panels, aligned 1xN comparisons, or 2x2 training summaries are all established patterns.
- Use `tight_layout()` or a deliberate `subplots_adjust(...)` when tight layout cannot preserve the intended panel geometry.
- Keep axes uncluttered. Grid lines are usually absent; legends are often frameless in transport plots.
- Parity panels use equal data limits and a diagonal reference. Log axes are appropriate for losses or quantities spanning orders of magnitude.
- Use inset axes only for a distinct scale or diagnostic that materially helps interpretation, such as a residual distribution or early-time behavior.

### Saving

- Preserve the existing flexible `save` argument position of each script; do not impose a global argument contract.
- Preserve per-script display DPI and save DPI. Saved figures commonly use PNG at 300 DPI, sometimes with `bbox_inches='tight'`, but DPI unification is explicitly rejected.
- Keep the established output filename and capitalization for an existing command. For a new command, use a stable, descriptive filename tied to the plot type.
- Interactive display remains the default where neighboring scripts behave that way; saving is opt-in through the routed command.

## Scientific Plotting Rules

- Never invent or silently change fit windows, smoothing windows, cutoffs, exclusions, component mappings, or unit conversions.
- Distinguish raw data, transformed data, fits, averages, and uncertainty visually and in labels.
- Put units on every dimensional axis and report derived quantities with their units.
- For parity plots, show the identity line and report the error metrics actually computed by the script. Do not imply validation quality from visual agreement alone.
- For convergence and transport plots, keep individual runs visible when they provide uncertainty context, then emphasize the aggregate result.
- Validate array shape, finite values, and expected columns before plotting. Stop with a user-facing error on missing or malformed inputs rather than emitting a misleading empty figure.

## Workflow for a New Plot

1. Identify the closest script in the same scientific family; reuse its figure grammar and CLI/save behavior.
2. Confirm inputs, columns, units, component mapping, fit/smoothing choices, and output name with project evidence or the user.
3. Follow the header, help, error-message, dependency, and import-order rules in `contributing.md`.
4. Separate data loading and validation from plotting helpers where that improves readability; add docstrings to helper functions.
5. Use the smallest layout that communicates the result and preserve raw-data context.
6. Test both interactive and save paths with representative data. Check labels, units, legends, clipping, empty data, NaN/Inf, and headless saving.
7. Run the validation checklist in `contributing.md` and update both English and Chinese user documentation for a user-visible feature.

## Review Checklist

- The plot matches a neighboring GPUMDkit workflow rather than introducing an unrelated house style.
- Color, marker, line, and alpha choices encode a documented distinction.
- Labels include quantities and units; legends are concise and do not obscure data.
- Identity, zero, fit, and uncertainty elements are visually distinguishable from observations.
- Scientific choices come from user input, documented defaults, or existing project settings.
- `-h` and missing-argument behavior work without optional/heavy dependencies for new scripts.
- Display and save paths both work; filename and DPI follow the target script rather than a global rewrite.
- No plot logic was changed under the guise of cosmetic cleanup.

## Dependency Reality

The plotting directory does not have one universal dependency set. Matplotlib and NumPy are common, while individual scripts may also require pandas, SciPy, seaborn, ASE, scikit-learn, UMAP, calorine, or `ferrodispcalc`. Inspect the target script before execution and keep optional imports after help handling in new scripts.
