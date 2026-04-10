# Reference ellipse for virtual community examples

A pre-calculated `nicheR_ellipsoid` object representing a hypothetical
species niche based on Annual Mean Temperature (bio_1) and Annual
Precipitation (bio_12).

## Usage

``` r
ref_ellipse
```

## Format

An object of class `nicheR_ellipsoid` (which is a `list`) with 13
elements:

- dimensions:

  Integer. Number of dimensions (2).

- var_names:

  Character vector. Names of variables (`"bio_1"`, `"bio_12"`).

- centroid:

  Named numeric vector. The center of the niche (\\\mu\\).

- cov_matrix:

  Matrix. The \\2 \times 2\\ covariance matrix (\\\Sigma\\).

- Sigma_inv:

  Matrix. The precision matrix (inverse covariance).

- chol_Sigma:

  Matrix. Cholesky decomposition of the covariance.

- eigen:

  List. Eigenvectors and eigenvalues of the covariance.

- cl:

  Numeric. Confidence level used (e.g., 0.99).

- chi2_cutoff:

  Numeric. The chi-square quantile for the given `cl`.

- semi_axes_lengths:

  Numeric vector. Radii of the ellipsoid axes.

- axes_coordinates:

  List. Vertices (endpoints) for each ellipsoid axis.

- volume:

  Numeric. The hyper-volume of the ellipsoid.

- cov_limits:

  List. Axis-aligned minimum and maximum limits.

## Details

This object serves as a template for testing community simulation
functions like
[`conserved_ellipses`](https://castanedam.github.io/nicheR/reference/conserved_ellipses.md).
It was generated using
[`build_ellipsoid`](https://castanedam.github.io/nicheR/reference/build_ellipsoid.md)
with a centroid at (23.75, 1750) and specific covariance structures to
reflect a typical temperature-precipitation relationship.

## Examples

``` r
data(ref_ellipse)
print(ref_ellipse)
#> nicheR Ellipsoid Object
#> -----------------------
#> Dimensions:        2D
#> Chi-square cutoff: 9.21
#> Centroid (mu):     23.5, 1750
#> 
#> Covariance matrix:
#>           bio_1 bio_12
#> bio_1     1.361   -100
#> bio_12 -100.000  62500
#> 
#> Ellipsoid semi-axis lengths:
#>   758.715, 3.326
#> 
#> Ellipsoid axis endpoints:
#>  Axis 1:
#>           bio_1   bio_12
#> vertex_a 24.714  991.286
#> vertex_b 22.286 2508.714
#> 
#>  Axis 2:
#>           bio_1   bio_12
#> vertex_a 26.826 1750.005
#> vertex_b 20.174 1749.995
#> 
#> Ellipsoid volume:  7927.882
#> 

# Access the volume
ref_ellipse$volume
#> [1] 7927.882
```
