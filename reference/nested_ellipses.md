# Generate nested ellipses based on a reference ellipse

Creates a sequence of nested ellipses by scaling the covariance matrix
of a reference ellipse. The distribution of the nested ellipses can be
controlled using a bias exponent to cluster them toward the border or
the centroid.

## Usage

``` r
nested_ellipses(object, n = 10, smallest_proportion = 0.1, bias = 1)
```

## Arguments

- object:

  An object of class "nicheR_ellipsoid" describing an initial ellipse.
  Must contain `centroid`, `cov_matrix`, and `cl`.

- n:

  Integer. Number of nested ellipses to generate. Default is 10.

- smallest_proportion:

  Numeric scalar in (0, 1). The scale of the smallest ellipse relative
  to the original. Default is `0.1`.

- bias:

  Numeric. An exponent controlling the spacing of the nested ellipses.

  - `bias = 1`: Linear spacing (default).

  - `0 < bias < 1`: Clusters ellipses toward the border (outer original
    ellipse).

  - `bias > 1`: Clusters ellipses toward the centroid (inner smallest
    ellipse).

## Value

An object of class `nicheR_community` containing the generated ellipses,
the reference object, and generation metadata.

## Details

The largest ellipse corresponds to the original reference, and the
smallest is that scaled by `smallest_proportion`.

The function generates a sequence of scale factors \\k\\ using the
formula: \\k_i = \text{smallest\\proportion} + (1 -
\text{smallest\\proportion}) \times t_i^{\text{bias}}\\, where \\t_i\\
is a linear sequence from 1 down to 0.
