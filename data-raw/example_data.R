# Title: Example data for nicheR
#
# Description: Rebuilds every example object shipped with the package, so the
# help pages and vignettes only have to load them.
#
# Not rebuilt here, raw inputs:
# inst/extdata/ma_bios.tif, WorldClim 2.1 bioclimatic variables
# inst/extdata/ma_biases.tif, species richness and nighttime lights

# SETUP ----

devtools::load_all()

bios <- terra::rast("inst/extdata/ma_bios.tif")
biases <- terra::rast("inst/extdata/ma_biases.tif")

# BACKGROUND TABLE ----

# The raster as a table, one row per cell that has data
back_data <- terra::as.data.frame(bios, xy = TRUE, na.rm = TRUE)
save(back_data, file = "data/back_data.rda", compress = "bzip2")

# REFERENCE NICHE ----

# The niche used across the vignettes. Mean annual temperature from 20 to 27
# and annual precipitation from 1000 to 2500, with a negative covariance.
ref_ranges <- data.frame(bio_1 = c(20, 27), bio_12 = c(1000, 2500))
ref_ellipse <- build_ellipsoid(ranges = ref_ranges, cl = 0.99, verbose = FALSE)
ref_ellipse <- update_ellipsoid_covariance(ref_ellipse, c("bio_1-bio_12" = -100))
save(ref_ellipse, file = "data/ref_ellipse.rda", compress = "bzip2")

# OTHER EXAMPLE NICHES ----

# Species 1 to 3 are the two-dimensional niches compared in the Build
# vignette. The ranges and covariances are the ones used there.
sp1_ranges <- data.frame(bio_1 = c(20, 32), bio_12 = c(750, 4000))
example_sp_1 <- build_ellipsoid(ranges = sp1_ranges, cl = 0.99, verbose = FALSE)
example_sp_1 <- update_ellipsoid_covariance(example_sp_1, c("bio_1-bio_12" = 750))

sp2_ranges <- data.frame(bio_1 = c(12, 24), bio_12 = c(1400, 4500))
example_sp_2 <- build_ellipsoid(ranges = sp2_ranges, cl = 0.99, verbose = FALSE)
example_sp_2 <- update_ellipsoid_covariance(example_sp_2, c("bio_1-bio_12" = 500))

sp3_ranges <- data.frame(bio_1 = c(18, 30), bio_12 = c(200, 1100))
example_sp_3 <- build_ellipsoid(ranges = sp3_ranges, cl = 0.99, verbose = FALSE)
example_sp_3 <- update_ellipsoid_covariance(example_sp_3, c("bio_1-bio_12" = 120))

# Species 4 is the three-dimensional example
sp4_ranges <- data.frame(bio_1 = c(22, 32), bio_12 = c(800, 4200),
	bio_15 = c(45, 115))
example_sp_4 <- build_ellipsoid(ranges = sp4_ranges, cl = 0.99, verbose = FALSE)
example_sp_4 <- update_ellipsoid_covariance(example_sp_4,
	c("bio_1-bio_12" = 200, "bio_1-bio_15" = 0, "bio_12-bio_15" = -5000))

save(example_sp_1, file = "data/example_sp_1.rda", compress = "bzip2")
save(example_sp_2, file = "data/example_sp_2.rda", compress = "bzip2")
save(example_sp_3, file = "data/example_sp_3.rda", compress = "bzip2")
save(example_sp_4, file = "data/example_sp_4.rda", compress = "bzip2")

# PREDICTIONS ----

# All four outputs of predict() over the raster, keeping the two variables
pred_rast <- predict(ref_ellipse, newdata = bios[[ref_ellipse$var_names]],
	include_mahalanobis = TRUE, include_suitability = TRUE,
	mahalanobis_truncated = TRUE, suitability_truncated = TRUE,
	keep_data = TRUE, verbose = FALSE)
terra::writeRaster(pred_rast, "inst/extdata/predictions_rast.tif",
	overwrite = TRUE)

# A virtual background of 1000 points drawn from the niche, then scored
virt <- virtual_data(ref_ellipse, n = 1000, truncate = FALSE,
	sampling = "centroid", seed = 1)
pred_virt <- predict(ref_ellipse, newdata = virt,
	include_mahalanobis = TRUE, include_suitability = TRUE,
	suitability_truncated = TRUE, keep_data = TRUE, verbose = FALSE)
write.csv(pred_virt, "inst/extdata/predictions_virt.csv", row.names = FALSE)

# BIASED PREDICTION ----

# Richness raises sampling effort and nighttime lights lower it
prep <- prepare_bias(bias_surface = biases,
	effect_direction = c("direct", "inverse"), verbose = FALSE)

# Suitability weighted toward the centroid, multiplied by the bias
biased <- apply_bias(prepared_bias = prep, prediction = pred_rast,
	prediction_layer = "suitability", sampling = "centroid",
	method = "suitability", verbose = FALSE)
terra::writeRaster(biased[[1]], "inst/extdata/applied_bias_rast.tif",
	overwrite = TRUE)

