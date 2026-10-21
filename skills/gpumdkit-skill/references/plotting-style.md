# GPUMDkit Plotting Style

Read this reference when creating a plotting script, changing the appearance of an existing plot, or reviewing a plot contribution. Also read `visualization.md`, `contributing.md`, and the reference governing the scientific quantity. This guide distills recurring choices in `Scripts/plt_scripts`; it does not replace script-specific behavior.

## Design Character

GPUMDkit figures are compact, publication-oriented scientific plots. They favor a white background, restrained Matplotlib colors, direct physical labels with units, and enough annotation to interpret a result without decorative UI elements. The data and its uncertainty are visually dominant.

Most plot scripts share one house style. It is written inline at the top of every script (a short `plt.rcParams.update({...})` block), so each script stays copyable on its own; do not move it into a shared module.

A few scripts keep their earlier look on purpose, for example `plt_xrd`, `plt_parity_density`, `plt_net_force`, `plt_vac`, `plt_pdos`, and the heat-transport scripts (`plt_emd`, `plt_emd2`, `plt_hnemd`, `plt_nemd`, `plt_viscosity`). Do not restyle them unless the maintainer asks, and do not copy their older conventions (inward ticks on all four sides, `pylab` imports, `tight_layout()`) into new scripts. When editing one of them, match the script itself.

## Visual Vocabulary

### Typography

- Use `sans-serif` with `Arial`, `DejaVu Sans`, and `Liberation Sans` as fallbacks, and `mathtext.fontset` set to `dejavusans`.
- House sizes: 9 pt tick labels, legends, and annotations; 10 pt axis labels; 12 pt panel labels such as `(a)` placed outside the upper-left corner of the panel.
- Capitalize the first word of axis labels and give explicit units, for example `DFT Energy (eV/atom)` and `Diffusion Coefficient (Å$^2$/ps)`.
- Use concise legends and annotations. Put fit statistics or derived values close to the relevant panel or curve.

### Color and line hierarchy

- Use saturated Matplotlib cycle colors (`C0`, `C1`, `C2`, ...) for related components. Directional data conventionally maps `x`, `y`, and `z` to `C0`, `C1`, and `C2` when the target script already follows that mapping. Avoid pale or yellowish choices unless the script already uses them (for example the directional Arrhenius plots keep their own palette).
- Do not use alpha on scatter points, legend handles, or ordinary data lines; open markers keep points distinct from fit lines. Use alpha only for filled uncertainty bands, KDE fills, and faint replicas drawn behind a thicker opaque mean, cumulative result, fit, or final trajectory.
- Use grey for references, raw replicas, and identity/zero lines. A grey dashed line (`0.55`, 1.2 pt) is the common parity reference.
- Preserve a script's existing color semantics. Never change colors in a way that changes the meaning of components between related figures.

### Axes and layout

- Choose the layout from the scientific comparison: compact single panels (about 4.2 x 3.5 in), aligned 1xN comparisons, or 2x2 training summaries (about 7.2 x 6.1 in) are all established patterns.
- Set margins and spacing with one explicit `fig.subplots_adjust(...)` line, commented so that users can tune it, instead of `tight_layout()`. Keep the left and bottom margins close to the right and top ones, and remember that matplotlib's `wspace` is the gap between columns and `hspace` the gap between rows.
- Ticks point outward on all axes; the frame is a thin 0.9 pt line. Keep axes uncluttered: no grid, and frameless legends placed where they hide no data (lower left is the usual choice).
- Parity panels use equal data limits and a diagonal reference. Log axes are appropriate for losses or quantities spanning orders of magnitude.
- Use inset axes only for a distinct scale or diagnostic that materially helps interpretation, such as a residual distribution or early-time behavior.

### Saving

- Preserve the existing flexible `save` argument position of each script; do not impose a global argument contract.
- Call `savefig(..., dpi=300)` without `bbox_inches="tight"`, so the geometry set by `subplots_adjust(...)` is kept in the saved file.
- Preserve existing scripts' display/save DPI. New plots use 150 DPI for display and at least 300 DPI for saved figures; Gallery examples use a separate 100 DPI export.
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

1. Identify the closest script in the same scientific family; reuse its figure grammar and CLI/save behavior, and copy the inline style block from a restyled script (for example `plt_train.py` or `plt_arrhenius_d.py`).
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
