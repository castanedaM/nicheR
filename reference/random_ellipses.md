# Generate random ellipses constrained by a point cloud and a reference ellipse

Creates n random ellipses with centroids sampled from an irregular point
cloud. Covariance matrices are built using random rotations and scaled
eigenvalues restricted by user-defined limits.

## Usage

``` r
random_ellipses(
  object,
  background,
  n = 10,
  smallest_proportion = 0.1,
  largest_proportion = 1,
  thin_background = FALSE,
  resolution = 50,
  seed = 1
)
```

## Arguments

- object:

  A nicheR_ellipsoid object used as a reference ellipse (the biggest to
  be generated), and containing at least `covariance_matrix` and `cl`.

- background:

  Matrix or Dataframe. The 2D point cloud (coordinates) used to select
  random centroids.

- n:

  Integer. Number of ellipses to generate.

- smallest_proportion:

  Numeric. Minimum scaling factor for the variance. Must be between 0
  and 1. Default = 0.1. This controls how much smaller the new ellipses
  can be compared to the reference.

- largest_proportion:

  Numeric. Maximum scaling factor for the variance relative to the
  reference. Default = 1.0. This controls how much larger the new
  ellipses can be compared to the reference. Values larger than 1 will
  result in ellipses that exceed the reference size.

- thin_background:

  Logical. If TRUE, centroids are sampled more uniformly across the
  background using a grid-based thinning approach. Default = FALSE.

- resolution:

  Integer. Number of cells per side in the grid to deal with point
  density variation across background. Default = 50.

- seed:

  Integer. Random seed for reproducibility. Default = 1. Set to NULL for
  no seeding.

## Value

An object of class `nicheR_community` containing the generated ellipses,
the reference object, and generation metadata.
