# nicheR_community Class Constructor

Helper function to construct a `nicheR_community` object.

## Usage

``` r
new_nicheR_community(
  ellipse_community,
  reference,
  pattern,
  n,
  smallest_proportion,
  largest_proportion = NA,
  bias = NA,
  thin_background = NA,
  resolution = NA,
  seed = NA
)
```

## Arguments

- ellipse_community:

  A list of `nicheR_ellipsoid` objects.

- reference:

  The `nicheR_ellipsoid` object used as a template.

- pattern:

  Character. The generation pattern: "random", "nested", or "conserved".

- n:

  Integer. Number of ellipses in the community.

- smallest_proportion:

  Numeric. Minimum scaling factor used.

- largest_proportion:

  Numeric. Maximum scaling factor used (NA if not applicable).

- bias:

  Numeric. Bias exponent used (NA if not applicable).

- thin_background:

  Logical. Whether the background was thinned during generation of the
  community (NA if not applicable).

- resolution:

  Numeric. Resolution of the background grid (NA if not applicable).

- seed:

  Integer. Seed used for reproducibility (NA if not applicable).

## Value

An object of class `nicheR_community`. The object aggregates a
collection of ellipsoids, the reference niche they were derived from,
and the metadata of the generation process.
