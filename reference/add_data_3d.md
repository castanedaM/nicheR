# Add data to an existing 3D E-space plot

Add data to an existing 3D E-space plot

## Usage

``` r
add_data_3d(data, dim = c(1, 2, 3), col_layer = NULL, alpha = 1, ...)
```

## Arguments

- data:

  A data frame or matrix containing the points.

- dim:

  Integer vector of length 3. Indices of dimensions to plot.

- col_layer:

  Character or `NULL`. Column for coloring.

- alpha:

  Transparency for the points. Default is `1`.

- ...:

  Additional arguments passed to
  [`rgl::points3d`](https://dmurdoch.github.io/rgl/dev/reference/primitives.html).
