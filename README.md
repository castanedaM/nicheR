
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/castanedaM/nicheR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/castanedaM/nicheR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# nicheR <img src='man/figures/nicheR_logo.png' align='right' height='25%' width='25%'/>

### Ellipsoid-based virtual niches in E-space and G-space

**nicheR** is an R package for building and visualizing
**ellipsoid-based ecological niches** using environmental data.

Inspired by the conceptual foundations of **NicheA** and the flexibility
of the **virtualspecies** package, **nicheR** provides a reproducible
framework that connects niche construction, prediction, sampling, and
visualization in one integrated workflow.

<br>

## Package description

**nicheR** operates across two complementary spaces:

- **E-space (Environmental Space)**: Ellipsoids represent multivariate
  environmental tolerances. This is where niches are defined,
  visualized, and compared.
- **G-space (Geographic Space)**: Predictions are projected across
  raster layers to map suitable geographic regions and generate virtual
  occurrence data.

This dual-space structure allows explicit separation between niche
definition, projection, and sampling processes, making it
straightforward to test ecological hypotheses under controlled,
reproducible conditions. The figure below shows a schematic view of how
the package works.

<br>

## Installing the package

Note: Internet connection is required to install the package.

The development version of nicheR can be installed using the code below.

``` r
# Installing and loading packages
if (!require("devtools")) install.packages("devtools")

# To install the package use
devtools::install_github("castanedaM/nicheR")

# To install the package and its vignettes use (if needed use: force = TRUE)
devtools::install_github("castanedaM/nicheR", build_vignettes = TRUE)
```

<br>

*Having problems?*

If you have any problems during installation, restart your R session,
close other RStudio sessions you may have open, and try again. If during
the installation you are asked to update packages, do so if you do not
need a specific version of one or more of the packages. If any package
gives an error when updating, install it alone using
`install.packages()`, then try installing nicheR again.

<br>

To load the package use:

``` r
library(nicheR)
```

<br>

## Workflow in nicheR

A typical nicheR workflow follows these steps:

1.  **Build an ellipsoid** — define the niche from environmental ranges
2.  **Predict suitability** — project the ellipsoid onto raster or data
    frame inputs
3.  **Prepare bias layers** *(optional)* — define and apply sampling
    biases
4.  **Sample occurrences** — draw virtual occurrence records using
    various strategies
5.  **Simulate communities** *(optional)* — generate multi-species
    virtual communities

A brief description of each step is presented below. For complete
demonstrations, see the package vignettes listed in [Checking the
vignettes](#checking-the-vignettes).

<br>

### Building an ellipsoid

The ellipsoid is the core object in nicheR. It represents the virtual
species’ niche as a multivariate ellipse in environmental space, defined
by a centroid and a covariance matrix that determine its position, size,
and orientation. Ellipsoids are built from a range data frame specifying
the minimum and maximum for each environmental variable. nicheR computes
the centroid, covariance matrix, semi-axes lengths, and niche volume
automatically.

nicheR includes base-R plotting functions for visualizing ellipsoids in
E-space, both for individual dimension pairs and all pairwise
combinations at once. Adding the environmental background contextualizes
where the niche sits relative to available conditions.

> For details on ellipsoid construction, covariance-based rotation,
> niche volume, and visualization options, see the [Building
> ellipsoids](https://castanedam.github.io/nicheR/articles/creating_ellipsoid_based_niches.html)
> vignette.

<br>

### Predicting suitability

Once the ellipsoid is built, `predict()` projects it onto environmental
data to produce Mahalanobis distance and suitability surfaces.
Suitability is highest at the niche centroid and decreases toward the
boundary. Truncated outputs set values to zero (suitability) or `NA`
(Mahalanobis) outside the ellipsoid, which is useful for enforcing
strict niche-boundary sampling downstream. Both `SpatRaster` and
`data.frame` inputs are supported and return matching output types.

> For a full guide to prediction outputs, truncation logic, raster vs.
> data frame workflows, and E-space vs. G-space visualization, see the
> [Predicting suitability and Mahalanobis
> distance](https://castanedam.github.io/nicheR/articles/predict.html)
> vignette.

<br>

### Preparing and applying bias layers

Bias layers represent external factors that influence sampling
detectability or effort, independent of the species’ actual niche. Each
layer is assigned a direction: `"direct"` means higher values increase
sampling probability, while `"inverse"` means higher values decrease it.
`prepare_bias()` standardizes all layers to \[0, 1\], applies the
directional transformations, and combines them into a single composite
bias surface. `apply_bias()` then weights the suitability prediction by
this surface to simulate non-random, biased occurrence sampling.

> For details on bias surface construction, directional effects, and how
> bias interacts with predictions and sampling, see the [Sampling bias
> data](https://castanedam.github.io/nicheR/articles/sampling_bias_data.html)
> vignette.

<br>

### Generating virtual occurrences

nicheR provides two sampling functions: `sample_data()` for unbiased
sampling directly from suitability or Mahalanobis surfaces, and
`sample_biased_data()` for sampling from a bias-weighted prediction
surface. Three sampling strategies are available: `"centroid"` samples
preferentially near the niche center, `"edge"` samples near the
boundary, and `"random"` samples uniformly across the suitable area. Two
weighting methods control sampling probability: `"suitability"` weights
by suitability score and `"mahalanobis"` weights by distance from the
centroid.

> For a full guide to sampling strategies, methods, strict
> vs. non-strict modes, and biased sampling workflows, see the [Sampling
> occurrence
> data](https://castanedam.github.io/nicheR/articles/sampling_occurrence_data.html)
> vignette.

<br>

### Simulating virtual communities

nicheR includes functions for simulating multi-species virtual
communities in environmental space. Three community types are available:
`random_ellipses()` generates species with niches distributed randomly
across E-space, `nested_ellipses()` generates nested communities where
each species niche is contained within the previous one, and
`conserved_ellipses()` generates communities where all species niches
remain similar to a reference niche, representing niche conservatism.
Community predictions produce presence-absence matrices and species
richness maps that can be used to explore macroecological patterns.

> For a full demonstration of community simulation, prediction, and
> richness mapping, see the [Virtual community
> simulation](https://castanedam.github.io/nicheR/articles/virtual_communities.html)
> vignette.

<br>

### Visualization

nicheR includes a set of plotting functions for visualizing ellipsoids,
predictions, and occurrence data in E-space. `plot_ellipsoid()` opens a
new E-space plot, `add_data()` overlays points colored by a continuous
variable, `add_ellipsoid()` overlays an ellipsoid boundary, and
`plot_ellipsoid_pairs()` produces multi-panel pairwise projections with
globally consistent axis limits.

> For a full guide to the plotting functions and all display options,
> see the [Visualizing ellipsoids in environmental
> space](https://castanedam.github.io/nicheR/articles/plotting_vignette.html)
> vignette.

<br>

## Checking the vignettes

Users can check nicheR vignettes for a full explanation of the package
functionality. The vignettes can be checked online at the [nicheR
site](https://castanedam.github.io/nicheR/) under the menu *Articles*.

To build the vignettes when installing the package from GitHub, make
sure to use the argument `build_vignettes = TRUE`.

Check each of the vignettes with the code below:

``` r
# Guide to building ellipsoid-based niches
vignette("creating_ellipsoid_based_niches")

# Guide to predicting suitability and Mahalanobis distance
vignette("predict")

# Guide to sampling occurrence data from virtual niches
vignette("sampling_occurrence_data")

# Guide to virtual occurrence sampling from simulated data
vignette("sampling_virtual_data")

# Guide to preparing and applying sampling bias
vignette("sampling_bias_data")

# Guide to simulating virtual communities
vignette("virtual_communities")

# Guide to visualizing ellipsoids in environmental space
vignette("plotting_vignette")
```

<br>

## Related work

Three tools have shaped virtual species and niche simulation in
ENM/SDMs:

- **NicheA** — a Java-based software that pioneered environmental-space
  niche visualization and explicitly linked niche theory with
  simulation. An important conceptual foundation for nicheR.
- **[virtualspecies](https://github.com/Farewe/virtualspecies)** — a
  widely adopted R package (200+ citations) for simulating virtual
  species and benchmarking SDMs and ENMs within R.
- **[evniche](https://github.com/marlonecobos/evniche/)** — a
  theoretically rigorous R package focused on environmental-space
  representations and ellipsoid-based niche concepts.

nicheR builds directly on these approaches, unifying their strengths in
a well-documented R package. A Shiny GUI for interactive niche
exploration is also under active development.

<br>

## Note on AI usage

To maintain high standards of code quality and documentation, we have
used AI LLM tools in this package. We used these tools for grammatical
polishing and exploring technical implementation strategies for
specialized functions. We manually checked and tested all code and
documentation refined with these tools.

<br>

## Contributing

We welcome contributions to improve `nicheR`. To maintain the integrity
and performance of the package, we follow a few core principles:

- **Quality over quantity**: we prioritize well-thought-out, stable
  improvements over frequent, minor changes. Please ensure your code is
  well-documented and follows the existing style of the package.
- **Minimal dependencies**: one of the goals of nicheR is to remain
  efficient. We prefer solutions that use base R or existing
  dependencies. Proposals that introduce new package dependencies will
  be evaluated for their necessity.
- **AI-assisted code**: if you use AI tools to generate code
  alternatives or improvements, please manually verify the logic and
  accuracy of the output and demonstrate the benefit in your Pull
  Request.
- **Testing**: new features should include examples, and tests should be
  performed to ensure they work as intended and do not break existing
  workflows.

If you have an idea for a major change, please open an Issue first to
discuss it with the maintainers.
