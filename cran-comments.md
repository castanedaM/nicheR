# nicheR CRAN Submission Comments

## nicheR version 0.2.0

This is an update of an existing package.

In this version the following changes were made:

- Added `run_gui()`, a Shiny interface to the package workflow. The
  application files are in `inst/shiny/` and are not evaluated during checks.
  The example is wrapped in `\dontrun{}` because the application blocks the
  console until closed and cannot run without a display.
- Added `update_ellipsoid_centroid()`, which translates an ellipsoid through
  environmental space without changing its shape, size, or orientation.
- [Changed the return value of `predict()` for `SpatRaster` input from a list
  of single-layer rasters to a single multi-layer `SpatRaster`, matching the
  layer names to the column names returned for data frame input. Code that
  indexed the result by name continues to work; code that treated it as a
  list may need updating.]
- Added the vignette `shiny_app`, documenting the interface.
- `shinydashboard` and `shinyjs` are declared in Imports. They are used by the
  application in `inst/shiny/`.

## Test environments

- Windows, R version 4.4.1 (2024-06-14 ucrt) (local)
- linux, R-devel (R-hub)
- m1-san, R-devel (R-hub)
- macos, R-devel (R-hub)
- macos-arm64, R-devel (R-hub)
- windows, R-devel (R-hub)
- macOS-latest, R release (GitHub Actions)
- Windows-latest, R release (GitHub Actions)
- Ubuntu-latest, R devel (GitHub Actions)
- Ubuntu-latest, R release (GitHub Actions)
- Ubuntu-latest, R oldrel-1 (GitHub Actions)

## R CMD check results

0 errors | 0 warnings | 2 notes

* Installed size is 9.2Mb, of which 6.5Mb is `doc/`. The vignettes document a
  workflow that is inherently visual, and the figures are needed to show what
  the ellipsoids and predictions look like.

* The note "unable to verify current time" reflects a network issue on the
  submitting machine and is unrelated to the package.

## Downstream dependencies

There are currently no downstream dependencies for this package.

---

# nicheR CRAN Submission Comments

## nicheR version 0.1.0 re-submission (second re-submission)

In this version the following changes were made:

- Added references in the DESCRIPTION file formatted as `authors (year) <doi:...>`.
- Updated vignettes to restore original `par` settings after plotting by adding `original_par <- par(no.readonly = TRUE)` at the beginning and `par(original_par)` at the end of each vignette.
- Updated documentation to add missing `\return` statements to functions with missing Rd-tags: `add_data_3d.Rd`, `add_ellipsoid_3d.Rd`, `plot_community.Rd`, and `plot_ellipsoid_3d.Rd`.
- Replaced `\dontrun{}` with `\donttest{}`. The reason for keeping `\donttest{}` is that `rgl` opens an interactive 3D device. The `requireNamespace()` guard ensures the example skips cleanly if `rgl` is not installed.

## Test environments

- Windows, R version 4.4.1 (2024-06-14 ucrt) (local)
- linux, R-devel (R-hub)
- m1-san, R-devel (R-hub)
- macos, R-devel (R-hub)
- macos-arm64, R-devel (R-hub)
- windows, R-devel (R-hub)
- macOS-latest, R release (GitHub Actions)
- Windows-latest, R release (GitHub Actions)
- Ubuntu-latest, R devel (GitHub Actions)
- Ubuntu-latest, R release (GitHub Actions)
- Ubuntu-latest, R oldrel-1 (GitHub Actions)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Downstream dependencies

There are currently no downstream dependencies for this package.

---

## nicheR version 0.1.0 re-submission (first re-submission)

In this version the following changes were made:

- Added single quotes around software names in Title and Description fields.
- No methods reference added: the associated manuscript is in preparation and not yet published.
- `NicheA`, `nicheR`, and `virtualspecies` in the DESCRIPTION are software package names, not misspellings.

## Test environments

- Windows, R version 4.4.1 (2024-06-14 ucrt) (local)
- linux, R-devel (R-hub)
- m1-san, R-devel (R-hub)
- macos, R-devel (R-hub)
- macos-arm64, R-devel (R-hub)
- windows, R-devel (R-hub)
- macOS-latest, R release (GitHub Actions)
- Windows-latest, R release (GitHub Actions)
- Ubuntu-latest, R devel (GitHub Actions)
- Ubuntu-latest, R release (GitHub Actions)
- Ubuntu-latest, R oldrel-1 (GitHub Actions)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Downstream dependencies

There are currently no downstream dependencies for this package.

---

## nicheR version 0.1.0 first submission

## Test environments

- Windows, R version 4.4.1 (2024-06-14 ucrt) (local)
- linux, R-devel (R-hub)
- m1-san, R-devel (R-hub)
- macos, R-devel (R-hub)
- macos-arm64, R-devel (R-hub)
- windows, R-devel (R-hub)
- macOS-latest, R release (GitHub Actions)
- Windows-latest, R release (GitHub Actions)
- Ubuntu-latest, R devel (GitHub Actions)
- Ubuntu-latest, R release (GitHub Actions)
- Ubuntu-latest, R oldrel-1 (GitHub Actions)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Downstream dependencies

There are currently no downstream dependencies for this package.
