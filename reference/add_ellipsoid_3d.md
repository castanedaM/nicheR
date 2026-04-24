# Add an ellipsoid to an existing 3D E-space plot

Add an ellipsoid to an existing 3D E-space plot

## Usage

``` r
add_ellipsoid_3d(
  object,
  dim = c(1, 2, 3),
  wire = FALSE,
  col_ell = "#800000",
  alpha_ell = 1,
  ...
)
```

## Arguments

- object:

  A `nicheR_ellipsoid` object.

- dim:

  Integer vector of length 3.

- wire:

  Logical. If `TRUE`, plots wireframe, otherwise plots a shaded volume.
  Default is `FALSE`.

- col_ell:

  Color of the ellipsoid. Default is `"#000000"`.

- alpha_ell:

  Transparency of the ellipsoid. Default is `1`. Not applied if
  `wire = TRUE`.

- ...:

  Additional arguments passed to
  [`rgl::wire3d`](https://dmurdoch.github.io/rgl/dev/reference/shade3d.html)
  or
  [`rgl::shade3d`](https://dmurdoch.github.io/rgl/dev/reference/shade3d.html).
