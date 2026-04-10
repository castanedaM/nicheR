# nicheR_ellipsoid Class Constructor

Internal helper to create a nicheR_ellipsoid object.

## Usage

``` r
new_nicheR_ellipsoid(
  dimensions,
  var_names,
  centroid,
  cov_matrix,
  Sigma_inv,
  chol_Sigma,
  eigen,
  cl,
  chi2_cutoff,
  semi_axes_lengths,
  axes_coordinates,
  volume,
  cov_limits
)
```

## Arguments

- dimensions:

  The number of variables/dimensions.

- var_names:

  Names of the original variables.

- centroid:

  The numeric center vector.

- cov_matrix:

  The covariance matrix.

- Sigma_inv:

  The precision matrix.

- chol_Sigma:

  The Cholesky decomposition.

- eigen:

  List of eigenvectors and eigenvalues.

- cl:

  Confidence level.

- chi2_cutoff:

  Chi-squared quantile.

- semi_axes_lengths:

  Radii of the ellipsoid.

- axes_coordinates:

  List of vertex matrices.

- volume:

  Calculated hyper-volume.

- cov_limits:

  Axis-aligned limits.

## Value

An object of class `nicheR_ellipsoid` with the fields described above.
