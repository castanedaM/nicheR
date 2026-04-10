# Compute variable ranges from data or statistics with optional expansion

These functions compute the minimum and maximum values for variables,
either directly from a dataset or based on normal distribution
parameters, with the ability to expand the resulting ranges by a
percentage.

## Usage

``` r
ranges_from_data(data, expand_min = NULL, expand_max = NULL)
ranges_from_stats(mean, sd, cl = 0.95, expand_min = NULL, expand_max = NULL)

ranges_from_data(data, expand_min = NULL, expand_max = NULL)

ranges_from_stats(mean, sd, cl = 0.95, expand_min = NULL, expand_max = NULL)
```

## Arguments

- data:

  A data.frame of at least two columns. Each column should contain
  numeric values.

- expand_min:

  A named vector or list of percentages (e.g., 10 for 10 defining how
  much to expand the minimum value of specific variables.

- expand_max:

  A named vector or list of percentages defining how much to expand the
  maximum value of specific variables.

- mean:

  A named numeric vector of mean values for each variable.

- sd:

  A named numeric vector of standard deviation values for each variable.
  Names must match those in \`mean\`.

- cl:

  A numeric value indicating the confidence level (default 0.95).

## Value

A data.frame with the (potentially expanded) minimum and maximum values
of each variable.

## Examples

``` r
# From data
df <- data.frame(var1 = c(0, 10), var2 = c(50, 100))
ranges_from_data(df, expand_min = list(var1 = 10),
                 expand_max = list(var2 = 20))
#>     var1 var2
#> min   -1   50
#> max   10  110
# From statistics
m <- c(var1 = 10, var2 = 100)
s <- c(var1 = 2, var2 = 15)
ranges_from_stats(mean = m, sd = s, cl = 0.95,
                  expand_min = list(var1 = 10))
#>          var1      var2
#> min  5.296086  70.60054
#> max 13.919928 129.39946
```
