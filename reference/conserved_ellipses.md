# Generate ellipses via multivariate normal biased sampling

Creates a set of ellipses with centroids sampled from a background,
biased by their proximity to the centroid to a reference niche. Includes
an option to thin the background to reduce centroid sampling bias due to
point-density.

## Usage

``` r
conserved_ellipses(
  object,
  background,
  n = 10,
  smallest_proportion = 0.1,
  largest_proportion = 1,
  thin_background = FALSE,
  resolution = 100,
  seed = 1
)
```

## Arguments

- object:

  A nicheR_ellipsoid object used as the reference. This is will be
  considered the "largest" ellipse to be generated.

- background:

  Matrix or Dataframe. The 2D point cloud (coordinates) used to select
  centroids for the ellipses.

- n:

  Integer. Number of ellipses to generate. Default = 10.

- smallest_proportion:

  Numeric scalar in `(0, 1)`. The scale of the smallest ellipse relative
  to the original. Default is `0.1`.

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
  density variation across background. Default = 100.

- seed:

  Integer. Random seed for reproducibility. Default = 1. Set to NULL for
  no seeding.

## Value

An object of class `nicheR_community` containing the generated ellipses,
the reference object, and generation metadata.

## Details

Ellipses are generated to simulate a community of niches with varying
degrees of similarity to the reference. The distribution of the
generated ellipses is influenced by the proximity to the reference and
the density of the background points.
