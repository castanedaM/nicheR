# Plot a nicheR ellipsoid in environmental space

Plots the 2D boundary of a `nicheR_ellipsoid` object in environmental
space for a chosen pair of dimensions. Optionally overlays background
points or a prediction surface colored by a continuous variable (e.g.,
suitability). Use
[`add_data`](https://castanedam.github.io/nicheR/reference/add_data.md)
and
[`add_ellipsoid`](https://castanedam.github.io/nicheR/reference/add_ellipsoid.md)
to layer additional data onto the plot after calling this function.

## Usage

``` r
plot_ellipsoid(
  object,
  background = NULL,
  prediction = NULL,
  dim = c(1, 2),
  col_layer = NULL,
  pal = hcl.colors(100, palette = "Viridis"),
  rev_pal = FALSE,
  bg_sample = NULL,
  lty = 1,
  lwd = 1,
  col_ell = "#000000",
  col_bg = "#8A8A8A",
  pch = 1,
  alpha_bg = 1,
  alpha_ell = 1,
  cex_ell = 1,
  cex_bg = 1,
  fixed_lims = NULL,
  ...
)
```

## Arguments

- object:

  A `nicheR_ellipsoid` object containing at least `centroid`,
  `cov_matrix`, `chi2_cutoff`, and `var_names`.

- background:

  Optional data frame or matrix of background points to plot behind the
  ellipsoid. Rows are observations, columns are environmental variables.
  If provided, `prediction` is ignored.

- prediction:

  Optional data frame or matrix of prediction values to plot. Used when
  `background` is `NULL`. Can be colored by a continuous variable using
  `col_layer`.

- dim:

  Integer vector of length 2. Indices of the two dimensions to plot.
  Default is `c(1, 2)`.

- col_layer:

  Character or `NULL`. Name of a column in `prediction` to use for
  coloring points by a continuous variable. If `NULL` (default), all
  prediction points are drawn with `col_bg`.

- pal:

  A color palette function or character vector used when `col_layer` is
  provided. Default is `hcl.colors(100, palette = "Viridis")`.

- rev_pal:

  Logical. If `TRUE`, reverses the color palette. Default is `FALSE`.

- bg_sample:

  Integer or `NULL`. If provided and the number of background or
  prediction rows exceeds this value, a random subsample of this size is
  drawn before plotting. Useful for large data frames. Default is `NULL`
  (plot all points).

- lty:

  Integer. Line type for the ellipsoid boundary. Default is `1` (solid).

- lwd:

  Numeric. Line width for the ellipsoid boundary. Default is `1`.

- col_ell:

  Character. Color of the ellipsoid boundary line. Default is
  `"#000000"` (black).

- col_bg:

  Character. Color of background or prediction points when `col_layer`
  is `NULL`. Default is `"#8A8A8A"` (grey).

- pch:

  Integer or character. Point symbol for background or prediction
  points. Default is `1`.

- alpha_bg:

  Numeric in `[0, 1]`. Transparency of background or prediction points.
  Default is `1` (fully opaque).

- alpha_ell:

  Numeric in `[0, 1]`. Transparency of the ellipsoid boundary line.
  Default is `1` (fully opaque).

- cex_ell:

  Numeric. Size scaling for the ellipsoid boundary line. Default is `1`.

- cex_bg:

  Numeric. Size scaling for background or prediction points. Default is
  `1`.

- fixed_lims:

  A named list with elements `xlim` and `ylim`, each a numeric vector of
  length 2. When provided, overrides the limits computed by
  [`safe_lims()`](https://castanedam.github.io/nicheR/reference/safe_lims.md).
  Intended for use by
  [`plot_ellipsoid_pairs`](https://castanedam.github.io/nicheR/reference/plot_ellipsoid_pairs.md)
  to enforce consistent axis scales across panels, but can also be set
  manually by the user. Default is `NULL` (limits computed from data).

- ...:

  Additional graphical parameters passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Called for its side effect of creating a plot. Returns `NULL` invisibly.

## Details

The function has three display modes depending on what is provided:

1.  **Background only** (`background` is not `NULL`): plots background
    points in `col_bg` with the ellipsoid boundary overlaid.

2.  **Prediction surface** (`background` is `NULL`, `prediction` is not
    `NULL`): plots prediction points, optionally colored by `col_layer`
    using values mapped onto `pal`. When `col_layer` is provided, points
    outside the ellipsoid (zero or `NA` in `col_layer`, as produced by
    truncated prediction types) are drawn in `col_bg` behind the colored
    interior points. Axis limits are computed from the full `prediction`
    extent so the view is never collapsed to the ellipsoid interior.

3.  **Ellipsoid only** (both `NULL`): plots the ellipsoid boundary alone
    with no background.

## See also

[`add_data`](https://castanedam.github.io/nicheR/reference/add_data.md)
to overlay occurrence points,
[`add_ellipsoid`](https://castanedam.github.io/nicheR/reference/add_ellipsoid.md)
to overlay additional ellipsoid boundaries,
[`plot_ellipsoid_pairs`](https://castanedam.github.io/nicheR/reference/plot_ellipsoid_pairs.md)
for pairwise plots of all dimensions.
