
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NicheR: Ecological Niche Modeling and Simulation in R

<!-- badges: start -->

<!-- badges: end -->

The goal of NicheR is to provide a robust set of tools for researchers
and modelers to:

Construct and define virtual niches using ellipsoid geometries in 2D or
3D environmental dimensions.

Identify and extract suitable environmental areas based on these defined
niches from various environmental data sources (rasters, data frames).

Simulate species occurrence points within suitable habitats, offering
different sampling strategies (random, biased towards center, or biased
towards edges).

Visualize niche boundaries and simulated occurrences effectively in both
environmental (E-space) and geographic (G-space).

Inspired by the methodologies in NicheA and the virtualspecies R
package, NicheR aims to streamline the process of niche
conceptualization and data generation for ecological studies.

## Authors

Mariana Castaneda-Guzman

Paanwaris Paansri

Connor Hughes

## Installation You can install the development version of NicheR from

GitHub with pak or remotes:

``` r
# 1. Install devtools if you don't have it yet
# if (!require("devtools")) {
#  install.packages("devtools")
# }

# 2. Install bean from GitHub
# devtools::install_github("castanedaM/NicheR")
```

## Dependencies NicheR relies on the following R packages, which will be

installed automatically with the package: dplyr, ggplot2, ggpubr,
plotly, RColorBrewer, and terra.

### Example This is a basic example demonstrating the core functionality

of NicheR.

### 1. Load Data and Define Environmental Background We’ll load NicheR and

prepare some example environmental raster data. For this, we’ll use
example files that should be located in your package’s inst/extdata
directory.

``` r
library(NicheR)
library(terra)
#> terra 1.8.60
library(ggplot2)
library(ggpubr)
#> 
#> Attaching package: 'ggpubr'
#> The following object is masked from 'package:terra':
#> 
#>     rotate
library(rgl)
library(plotly)
#> 
#> Attaching package: 'plotly'
#> The following object is masked from 'package:ggplot2':
#> 
#>     last_plot
#> The following object is masked from 'package:stats':
#> 
#>     filter
#> The following object is masked from 'package:graphics':
#> 
#>     layout

# Access example environmental raster files from the package
# Make sure "Bio1.tif", "Bio4.tif", "Bio12.tif" exist in inst/extdata
bio_1 <- terra::rast(system.file("extdata", "Bio1.tif", package = "NicheR"))
bio_4 <- terra::rast(system.file("extdata", "Bio4.tif", package = "NicheR"))
bio_12 <- terra::rast(system.file("extdata", "Bio12.tif", package = "NicheR"))

# Stack and name the environmental layers
env_stack <- c(bio_1, bio_12, bio_4)
names(env_stack) <- c("mean_temp", "temp_seasonality", "annual_precip")

# Convert the raster stack to a data frame for E-space plotting
env_df <- as.data.frame(env_stack, xy = TRUE, na.rm = TRUE)

# Downsample for faster plotting if the background is very large
if(nrow(env_df) > 20000){
  set.seed(1234)
  env_df_plotting <- env_df[sample(1:nrow(env_df), size = 20000, replace = FALSE), ]
} else {
  env_df_plotting <- env_df
}
```

### 2. Build a 3D Ellipsoid Niche The build_ellipsoid() function creates

the core niche object by specifying its center, semi-axis lengths, and
rotation angles.

``` r
# Define niche parameters
niche_center <- c(mean_temp = 280, temp_seasonality = 1500, annual_precip = 2000)
niche_axes <- c(100, 500, 2500) # Semi-axis lengths (breadth)
niche_angles <- c(0, 0, pi/12) # Small rotation around Z-axis

# Create the ellipsoid niche object
my_niche <- build_ellps(
  center = niche_center,
  axes = niche_axes,
  angles = niche_angles,
  n_points = 50 # Resolution for surface points
)
```

### 3. Visualize the Niche in Environmental Space plot_e_space() allows you

to visualize the environmental background and the niche ellipsoid. It
can generate both static 2D pairwise plots and interactive 3D plots.

#### 3.1. 2D Pairwise Plots

``` r
plot_e_space(
  env_bg = env_df_plotting,
  x = "mean_temp",
  y = "temp_seasonality",
  z = "annual_precip",
  labels = c("Mean Temperature (°C)", "Temp. Seasonality", "Annual Precipitation (mm)"),
  niche = my_niche,
  plot.3d = FALSE
)
#> Sampling 10000 of 20000 rows from 'env_bg' for plotting.
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

#### 3.2. 3D Interactive Plot

``` r
# This code will generate an interactive plotly plot,
# which cannot be directly embedded in static markdown but is great for exploration.
plot_e_space(
  env_bg = env_df_plotting,
  x = "mean_temp",
  y = "temp_seasonality",
  z = "annual_precip",
  labels = c("Mean Temperature (°C)", "Temp. Seasonality", "Annual Precipitation (mm)"),
  niche = my_niche,
  plot.3d = TRUE
)
```

### 4. Extract Suitable Environments and Sample Occurrences The

get_suitable_env() function identifies all environmental points that
fall within your niche, while get_sample_occ() samples a subset of these
points with a specified bias.

``` r
# Extract suitable environmental areas as a spatial raster
suitable_g_space <- get_suitable_env(
  niche = my_niche,
  env_bg = env_stack, # Use the raster stack for spatial output
  out = "spatial"
)

# Plot the suitable area on a map
plot(suitable_g_space, main = "Suitable Geographic Environment")
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

Now, let’s sample some occurrences within this suitable area.

``` r
# Sample 500 occurrences, biased towards the center of the niche
sampled_occ_center <- get_sample_occ(
  n_occ = 500,
  niche = my_niche,
  env_bg = env_stack, # Can use raster or data.frame here
  method = "center",
  seed = 42
)
#> Done sampling 500 occurrences

plot(suitable_g_space, main = "Suitable Geographic Environment")
# Overlay sampled points on the geographic suitable area map
points(sampled_occ_center$x, sampled_occ_center$y,
       pch = 20, col = "red", cex = 0.6)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

Visualize Sampled Occurrences in Environmental Space You can integrate
these simulated occurrences back into your environmental space plots.

``` r
plot_e_space(
  env_bg = env_df_plotting,
  x = "mean_temp",
  y = "temp_seasonality",
  z = "annual_precip",
  niche = my_niche,
  occ_pts = sampled_occ_center,
  show.pts.in = TRUE, # Highlight background points within niche
  show.occ.density = TRUE # Show marginal density plots for occurrences
)
#> Sampling 10000 of 20000 rows from 'env_bg' for plotting.
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

### Contributing

We welcome contributions! If you have suggestions, bug reports, or would
like to contribute code, please feel free to open an issue or submit a
pull request on the NicheR GitHub repository.
