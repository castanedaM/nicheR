# Plot a nicheR ellipsoid in 3D environmental space

Creates an interactive 3D plot of a `nicheR_ellipsoid` object with
support for background points or suitability surfaces.

## Usage

``` r
plot_ellipsoid_3d(
  object,
  dim = c(1, 2, 3),
  wire = FALSE,
  aspect = TRUE,
  background = NULL,
  prediction = NULL,
  col_layer = NULL,
  pal = hcl.colors(100, palette = "Viridis"),
  rev_pal = FALSE,
  bg_sample = NULL,
  col_ell = "#8b0000",
  alpha_ell = 1,
  alpha_bg = 1,
  col_bg = "#8A8A8A",
  fixed_lims = NULL,
  xlab = NULL,
  ylab = NULL,
  zlab = NULL,
  ...
)
```

## Arguments

- object:

  A `nicheR_ellipsoid` object constructed with at least 3 dimensions.

- dim:

  Integer vector of length 3. Indices of dimensions to plot.

- wire:

  Logical. If `TRUE`, plots wireframe. The default, `FALSE`, plots a
  shaded volume.

- aspect:

  Logical. If `TRUE`, maintains aspect ratio (1:1:1).

- background:

  Optional data frame/matrix of background points. This argument has
  priority over `prediction` if both are provided.

- prediction:

  Optional data frame/matrix for prediction surfaces.

- col_layer:

  Character or `NULL`. Column in `prediction` to use for coloring.

- pal:

  Color palette function or character vector.

- rev_pal:

  Logical. If `TRUE`, reverses the palette.

- bg_sample:

  Integer or `NULL`. Subsample size for large data.

- col_ell:

  Color of the ellipsoid. Default is `"#000000"`.

- alpha_ell:

  Transparency of the ellipsoid. Default is `1`. Not applied to
  wireframe mode.

- alpha_bg:

  Transparency of background points. Default is `1`. Also applied to
  prediction points.

- col_bg:

  Color for background points.

- fixed_lims:

  Named list with `xlim`, `ylim`, and `zlim`.

- xlab, :

  ylab, zlab Axis labels. The default, `NULL`, uses the variable names.

- ...:

  Additional graphical parameters.
