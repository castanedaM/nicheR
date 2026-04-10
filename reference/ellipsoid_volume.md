# Compute ellipsoid hypervolume

Computes the geometric volume (area in 2D, volume in 3D, hypervolume in
higher dimensions) of a \\p\\-dimensional ellipsoid defined by its
semi-axis lengths.

## Usage

``` r
ellipsoid_volume(n_dimensions, semi_axes_lengths)
```

## Arguments

- n_dimensions:

  Integer. Number of dimensions \\p\\.

- semi_axes_lengths:

  Numeric vector of length `n_dimensions` containing the ellipsoid
  semi-axis lengths.

## Value

Numeric. Geometric volume (or hypervolume) of the ellipsoid.

## Details

For semi-axes \\a_1, \dots, a_p\\, the volume is: \$\$ V_p =
\frac{\pi^{p/2}}{\Gamma(p/2 + 1)} \prod\_{i=1}^{p} a_i \$\$ where
\\\pi^{p/2} / \Gamma(p/2 + 1)\\ is the volume of the unit
\\p\\-dimensional ball.

In probabilistic niche models, semi-axis lengths are typically derived
from covariance eigenvalues and a chi-square cutoff.

## See also

[`build_ellipsoid`](https://castanedam.github.io/nicheR/reference/build_ellipsoid.md),
[`ellipsoid_calculator`](https://castanedam.github.io/nicheR/reference/ellipsoid_calculator.md)
