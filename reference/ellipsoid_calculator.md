# Calculate n-dimensional ellipsoid metrics

Computes geometric and probabilistic metrics for an n-dimensional
ellipsoid defined by a centroid and covariance matrix, including
semi-axis lengths, axis vertices, and hypervolume for a chi-square
confidence contour.

## Usage

``` r
ellipsoid_calculator(cov_matrix, centroid, cl, verbose = TRUE)
```

## Arguments

- cov_matrix:

  A square, numeric covariance matrix \\\Sigma\\. Must be SPD.
  Row/column names (if provided) are used as variable names in the
  output.

- centroid:

  Numeric vector giving the centroid \\\mu\\. Must have length equal to
  `ncol(cov_matrix)`.

- cl:

  Numeric confidence level in (0, 1). Used to compute the chi-square
  cutoff defining the ellipsoid contour.

- verbose:

  Logical; if `TRUE`, prints progress messages.

## Value

An object of class `"nicheR_ellipsoid"` created by
[`new_nicheR_ellipsoid`](https://castanedam.github.io/nicheR/reference/new_nicheR_ellipsoid.md),
containing ellipsoid geometry and associated quantities (e.g., centroid,
covariance matrix, chi-square cutoff, semi-axis lengths, axis vertex
coordinates, volume, and covariance limits).

## Details

The ellipsoid boundary is defined by the constant Mahalanobis distance
contour: \$\$(x - \mu)^\top \Sigma^{-1} (x - \mu) = c^2,\$\$ where
\\\mu\\ is the centroid, \\\Sigma\\ is the covariance matrix, and \\c^2
= \chi^2\_{n}(\mathrm{cl})\\ is the chi-square cutoff with \\n\\ degrees
of freedom.

The covariance matrix must be symmetric positive definite (SPD). The
inverse covariance is computed via the Cholesky factorization. Semi-axis
lengths are computed from covariance eigenvalues \\\lambda_i\\ as:
\$\$a_i = \sqrt{\lambda_i c^2}.\$\$

Axis vertices are computed along each eigenvector direction as \\\mu \pm
a_i v_i\\. Hypervolume is computed with
[`ellipsoid_volume`](https://castanedam.github.io/nicheR/reference/ellipsoid_volume.md),
and covariance-derived limits with
[`covariance_limits`](https://castanedam.github.io/nicheR/reference/covariance_limits.md).

## Examples

``` r
cm <- matrix(c(11.11, 0,
               0, 17777.78), nrow = 2, byrow = TRUE)
colnames(cm) <- rownames(cm) <- c("var1", "var2")
ctr <- c(20, 600)
ell <- ellipsoid_calculator(cov_matrix = cm, centroid = ctr, cl = 0.95, verbose = FALSE)
```
