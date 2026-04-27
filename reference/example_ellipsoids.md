# Example niche ellipsoid objects for virtual communities

Pre-calculated `nicheR_ellipsoid` objects representing hypothetical
species niches based on bioclimatic variables. These objects are
primarily used in the package vignettes and examples to demonstrate
ellipsoid creation, covariance adjustments, visual comparisons, and
multidimensional niches.

## Usage

``` r
example_sp_1

example_sp_2

example_sp_3

example_sp_4
```

## Format

Objects of class `nicheR_ellipsoid` (which are `list`s) with 13
elements:

- dimensions:

  Integer. Number of dimensions (2 or 3).

- var_names:

  Character vector. Names of variables (e.g., `"bio_1"`, `"bio_12"`).

- centroid:

  Named numeric vector. The center of the niche (\\\mu\\).

- cov_matrix:

  Matrix. The covariance matrix (\\\Sigma\\).

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

  List. Axis-aligned minimum and maximum covariance limits.

An object of class `nicheR_ellipsoid` of length 13.

An object of class `nicheR_ellipsoid` of length 13.

An object of class `nicheR_ellipsoid` of length 14.

## Details

These objects serve as templates for testing community simulation and
projection functions. They were generated using
[`build_ellipsoid`](https://castanedam.github.io/nicheR/reference/build_ellipsoid.md)
and modified with
[`update_ellipsoid_covariance`](https://castanedam.github.io/nicheR/reference/update_ellipsoid_covariance.md)
to reflect distinct ecological strategies:

- **`example_sp_1`**: A 2D niche (`bio_1`, `bio_12`) with a broad,
  warm-climate preference and wide precipitation tolerance. It includes
  a positive covariance (750) between temperature and precipitation.

- **`example_sp_2`**: A 2D niche (`bio_1`, `bio_12`) shifted toward
  cooler and wetter environments compared to `example_sp_1`, with a
  positive covariance (500).

- **`example_sp_3`**: A 2D niche (`bio_1`, `bio_12`) representing a
  warm-adapted specialist restricted to dry environments. It has a much
  smaller niche volume and a slight positive covariance (120).

- **`example_sp_4`**: A 3D niche (`bio_1`, `bio_12`, `bio_15`) adapted
  to seasonally pulsed precipitation environments. It features a strong
  negative covariance (-5000) between annual precipitation and
  precipitation seasonality.

## See also

[`build_ellipsoid`](https://castanedam.github.io/nicheR/reference/build_ellipsoid.md),
[`update_ellipsoid_covariance`](https://castanedam.github.io/nicheR/reference/update_ellipsoid_covariance.md)

## Examples

``` r
data(example_sp_1)
print(example_sp_1)
#> nicheR Ellipsoid Object
#> -----------------------
#> Dimensions:        2D
#> Chi-square cutoff: 9.21
#> Centroid (mu):     26, 2375
#> 
#> Covariance matrix:
#>        bio_1   bio_12
#> bio_1      4    750.0
#> bio_12   750 293402.8
#> 
#> Covariance Limits:
#>                    min      max
#> bio_1-bio_12 -1083.333 1083.333
#> 
#> Ellipsoid semi-axis lengths:
#>   1643.885, 4.38
#> 
#> Ellipsoid axis endpoints:
#>  Axis 1:
#>           bio_1   bio_12
#> vertex_a 21.798  731.121
#> vertex_b 30.202 4018.879
#> 
#>  Axis 2:
#>          bio_1   bio_12
#> vertex_a 30.38 2374.989
#> vertex_b 21.62 2375.011
#> 
#> Ellipsoid volume:  22619.64
#> 

# Access the volume of the 2D broad warm-climate niche
example_sp_1$volume
#> [1] 22619.64

# Access the covariance matrix of the 3D seasonal precipitation niche
data(example_sp_4)
example_sp_4$cov_matrix
#>             bio_1   bio_12     bio_15
#> bio_1    2.777778    200.0     0.0000
#> bio_12 200.000000 321111.1 -5000.0000
#> bio_15   0.000000  -5000.0   136.1111
```
