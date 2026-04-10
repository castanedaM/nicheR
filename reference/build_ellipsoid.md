# Build a probabilistic ellipsoidal niche from ranges

Builds an ellipsoidal niche in multivariate environmental space using a
multivariate normal (MVN) contour defined by a constant Mahalanobis
distance. The ellipsoid is parameterized from user-provided variable
ranges by deriving a centroid and marginal standard deviations, and
assuming a diagonal covariance matrix.

## Usage

``` r
build_ellipsoid(range, cl = 0.99, verbose = TRUE)
```

## Arguments

- range:

  A 2-row `matrix` or `data.frame` of bounds, with variables as columns.
  Rows may be ordered as min/max or max/min. Column names are required
  and used as variable names.

- cl:

  Numeric confidence level in (0, 1). Used to compute the chi-square
  cutoff defining the ellipsoid contour.

- verbose:

  Logical; if `TRUE`, prints brief progress messages.

## Value

An object of class `"nicheR_ellipsoid"` produced by
[`ellipsoid_calculator`](https://castanedam.github.io/nicheR/reference/ellipsoid_calculator.md)
(via
[`new_nicheR_ellipsoid`](https://castanedam.github.io/nicheR/reference/new_nicheR_ellipsoid.md)),
containing ellipsoid geometry and associated quantities (e.g., centroid,
covariance matrix, chi-square cutoff, semi-axis lengths, axis vertex
coordinates, volume, and covariance limits).

## Details

`range` must be a 2-row `matrix` or `data.frame` with variables in
columns. Rows represent lower and upper bounds for each variable (row
order may be min/max or max/min). The centroid is computed as: \$\$\mu_i
= (m_i + M_i)/2.\$\$

Marginal standard deviations are derived assuming bounds represent
approximately \\\pm 3\\ standard deviations: \$\$\sigma_i = (M_i -
m_i)/6,\$\$ and a diagonal covariance matrix is assumed: \$\$\Sigma =
\mathrm{diag}(\sigma_1^2,\dots,\sigma_n^2).\$\$

The ellipsoid contour is defined using a chi-square cutoff \\c^2 =
\chi^2\_{n}(\mathrm{cl})\\, where \\n\\ is the number of variables.

## Examples

``` r
rng <- data.frame(bio1 = c(10, 20),
                  bio2 = c(20, 30))
ell <- build_ellipsoid(range = rng, cl = 0.95, verbose = FALSE)
print(ell)
#> nicheR Ellipsoid Object
#> -----------------------
#> Dimensions:        2D
#> Chi-square cutoff: 5.991
#> Centroid (mu):     15, 25
#> 
#> Covariance matrix:
#>       bio1  bio2
#> bio1 2.778 0.000
#> bio2 0.000 2.778
#> 
#> Ellipsoid semi-axis lengths:
#>   4.08, 4.08
#> 
#> Ellipsoid axis endpoints:
#>  Axis 1:
#>          bio1  bio2
#> vertex_a   15 20.92
#> vertex_b   15 29.08
#> 
#>  Axis 2:
#>           bio1 bio2
#> vertex_a 19.08   25
#> vertex_b 10.92   25
#> 
#> Ellipsoid volume:  52.285
#> 
```
