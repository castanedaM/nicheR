# Title: nicheR tester functions
# Authors: Mariana Castaneda-Guzman
# Date Last Updated: 2/11/2026
# Description: Tester workflows for function in the nicheR package


# Load function from local package after any update
devtools::load_all()
set.seed(124)


# build_ellipsoid() tests -------------------------------------------------


# - test range argument inputs --------------------------------------------


# Tests for range

# Test 1: basic run
rng_t1 <- data.frame(var1 = c(10, 20),
                     var2 = c(20, 35),
                     var3 = c(67, 80))

ell_t1 <- build_ellipsoid(range = rng_t1)
plot_nicheR(list(ell_t1))


# Test 2: Wrong range input (min max switch), should fail
rng_t2 <- data.frame(var1 = c(20, 10),
                     var2 = c(20, 35),
                     var3 = c(67, 80))

ell_t2 <- build_ellipsoid(range = rng_t2)


# Test 3: Wrong range input (missing row)
rng_t3 <- data.frame(var1 = c(10, 20),
                     var2 = c(20),
                     var3 = c(67, 80))

ell_t3 <- build_ellipsoid(range = rng_t3)


# Test 4: Wrong range input (not a 2-row data.frame or matrix)
rng_t4 <- data.frame(var1 = c(10, 20, 30),
                     var2 = c(20, 35, 40),
                     var3 = c(67, 80, 85))

ell_t4 <- build_ellipsoid(range = rng_t4)


# Test 5: Range with NA values
rng_t5 <- data.frame(var1 = c(10, NA),
                     var2 = c(20, 35),
                     var3 = c(67, 80))

ell_t5 <- build_ellipsoid(range = rng_t5)


# Test 6: Range with non-numeric values, should still work
rng_t6 <- data.frame(var1 = c(10, 20),
                     var2 = c("20", "35"),
                     var3 = c(67, 80))

ell_t6 <- build_ellipsoid(range = rng_t6)
ell_t6$centroid # check if var2 is 27.5

# Test 7: Range provided as a matrix (should still work)
rng_t7 <- matrix(c(10, 20,
                   20, 35,
                   67, 80),
                 nrow = 2, byrow = FALSE)

colnames(rng_t7) <- c("var1", "var2", "var3")

ell_t7 <- build_ellipsoid(range = rng_t7)
plot_nicheR(list(ell_t7))


# Test 8: Range without colnames (should not work, variables do not auto-named)
rng_t8 <- matrix(c(10, 20,
                   20, 35,
                   67, 80),
                 nrow = 2, byrow = FALSE)

ell_t8 <- build_ellipsoid(range = rng_t8)


# Test 9: Zero-width variable (min == max)
rng_t9 <- data.frame(var1 = c(10, 10),
                     var2 = c(20, 35),
                     var3 = c(67, 80))

ell_t9 <- build_ellipsoid(range = rng_t9)


# Test 10: Very small ranges (potential numeric stability)
rng_t10 <- data.frame(var1 = c(0.000001, 0.000002),
                      var2 = c(0.000010, 0.000030),
                      var3 = c(0.000067, 0.000080))

ell_t10 <- build_ellipsoid(range = rng_t10)
plot_nicheR(list(ell_t10))


# Test 11: Negative values (should be fine if numeric)
rng_t11 <- data.frame(var1 = c(-10, -2),
                      var2 = c(-5, 5),
                      var3 = c(67, 80))

ell_t11 <- build_ellipsoid(range = rng_t11)
plot_nicheR(list(ell_t11))


# Test 12: Range provided as vector (wrong type)
rng_t12 <- c(10, 20)

ell_t12 <- build_ellipsoid(range = rng_t12)


# Test 13: Range has infinite values
rng_t13 <- data.frame(var1 = c(10, Inf),
                      var2 = c(20, 35),
                      var3 = c(67, 80))

ell_t13 <- build_ellipsoid(range = rng_t13)


# Test 14: Duplicate column names
rng_t14 <- data.frame(var1 = c(10, 20),
                      var1 = c(20, 35),
                      var3 = c(67, 80))
# R renames column automatically to var1.1

ell_t14 <- build_ellipsoid(range = rng_t14)
ell_t14$var_names

# Test 15: Large magnitude values (scale stress test)
rng_t15 <- data.frame(var1 = c(1e6, 2e6),
                      var2 = c(2e6, 3.5e6),
                      var3 = c(6.7e6, 8.0e6))

ell_t15 <- build_ellipsoid(range = rng_t15)
plot_nicheR(list(ell_t15))


# - test cov_matrix argument inputs ---------------------------------------

# Test 16: covariance with wrong dimensions
rng_t16 <- data.frame(var1 = c(10, 20),
                      var2 = c(20, 35),
                      var3 = c(67, 80))

cov_t16 <- diag(2)

ell_t16 <- build_ellipsoid(range = rng_t16, cov_matrix = cov_t16)


# Test 17: covariance not symmetric
rng_t17 <- data.frame(var1 = c(10, 20),
                      var2 = c(20, 35),
                      var3 = c(67, 80))

cov_t17 <- matrix(c(1.0, 0.9, 0.1,
                    0.2, 2.0, 0.3,
                    0.1, 0.3, 1.5),
                  nrow = 3, byrow = TRUE)

colnames(cov_t17) <- c("var1", "var2", "var3")
rownames(cov_t17) <- c("var1", "var2", "var3")

ell_t17 <- build_ellipsoid(range = rng_t17, cov_matrix = cov_t17)

# Test 18: covariance not positive definite
rng_t18 <- data.frame(var1 = c(10, 20),
                      var2 = c(20, 35),
                      var3 = c(67, 80))

cov_t18 <- matrix(c(1, 2, 3,
                    2, 4, 6,
                    3, 6, 9),
                  nrow = 3, byrow = TRUE)  # singular

colnames(cov_t18) <- c("var1", "var2", "var3")
rownames(cov_t18) <- c("var1", "var2", "var3")

ell_t18 <- build_ellipsoid(range = rng_t18, cov_matrix = cov_t18)




# - test update_covariance_ellipsoid() --------------------------------------

# Test 19: Basic run of covariance limits
rng_t19 <- data.frame(var1 = c(10, 20),
                      var2 = c(20, 35),
                      var3 = c(67, 80))

ell_t19A <- build_ellipsoid(range = rng_t19)

ell_t19A$cov_limits


ell_t19B <- update_ellipsoid_covariance(ell_t19A,
                                        covariance = c("var1-var2" = 2))

ell_t19B$cov_limits

plot_nicheR(list(ell_t19A, ell_t19B))
# TO DO: change lines

# Ecological Interpretations
# axes interpretations with change in covariance = c("var1-var2" = 2):
# axis 1: The species tolerates low Var1 only when Var2 is also low, and high Var1 only when Var2 is also high.
# axis 2: Tolerance in Var3 is independent of Var1 and Var2.
# axis 3: This axis represents combinations the species tolerates much less, because Var1 and Var2 move against each other.

# “Introducing positive covariance between Var1 and Var2 rotates the niche, such
# that tolerance is greatest when both variables increase or decrease together,
# and most restricted when they vary in opposite directions.”


# Test 20: Wrong covariance name (pair not found)
rng_t20 <- data.frame(var1 = c(10, 20),
                      var2 = c(20, 35),
                      var3 = c(67, 80))

ell_t20A <- build_ellipsoid(range = rng_t20)

ell_t20B <- update_ellipsoid_covariance(ell_t20A,
                                        covariance = c("var1-var4" = 2))

# Test 21: Wrong covariance input type (not numeric)
rng_t21 <- data.frame(var1 = c(10, 20),
                      var2 = c(20, 35),
                      var3 = c(67, 80))

ell_t21A <- build_ellipsoid(range = rng_t21)

ell_t21B <- update_ellipsoid_covariance(ell_t21A,
                                        covariance = c("var1-var2" = "high"))



# predict.nicheR_ellipsoid() tests ----------------------------------------

# TO DO: add argument keep_data = TRUE, when data.frame input
# To do: check imports of packages, to not call on the ones the package installs

# Test 22: Basic raster prediction
library(geodata)
# library(sf)

wc_t22 <- geodata::worldclim_global(var = "bio",
                                    res = 10,
                                    path = tempdir())

bios_t22 <- wc_t22[[c(1, 12, 15)]]
names(bios_t22) <- c("bio1", "bio12", "bio15")

rng_t22 <- data.frame(bio1 = c(-10, 10),
                        bio12 = c(500, 2000),
                        bio15 = c(30, 150))

ell_t22 <- build_ellipsoid(range = rng_t22)

pred_t22 <- predict(ell_t22,
                    newdata = bios_t22,
                    include_mahalanobis = TRUE,
                    include_suitability = TRUE)

plot_nicheR(list(ell_t22))
plot(pred_t22, col=rev(viridis::viridis(100)))
plot(log(pred_t22), col=rev(viridis::viridis(100)))


# Test 23: Prediction using data.frame input
bios_df_t23 <- as.data.frame(bios_t22, xy = TRUE)

pred_t23 <- predict(ell_t22,
                    newdata = bios_t22,
                    include_mahalanobis = TRUE,
                    include_suitability = TRUE,
                    mahalanobis_truncated = TRUE,
                    suitability_truncated = TRUE)

head(pred_t23)
plot(log(pred_t23), col=rev(viridis::viridis(100)))
plot(pred_t23, col=rev(viridis::viridis(100)))


# Test 24: Missing variable in newdata
bios_t24 <- bios_t22[[c("bio1", "bio12")]]  # bio15 removed

pred_t24 <- predict(ell_t22,
                    newdata = bios_t24,
                    include_mahalanobis = TRUE)


# Test 25: Extra variable in newdata
bios_df_t25 <- bios_df_t23
bios_df_t25$extra_var <- rnorm(nrow(bios_df_t25))

pred_t25 <- predict(ell_t22,
                    newdata = bios_df_t25,
                    include_mahalanobis = TRUE)

head(pred_t25)


# Test 26: Truncated suitability
pred_t26 <- predict(ell_t22,
                    newdata = bios_t22,
                    include_suitability = TRUE,
                    suitability_truncated = TRUE)

plot(pred_t26)
plot(log(pred_t26), col=rev(viridis::viridis(100)))


# Test 27: Mahalanobis only
pred_t27 <- predict(ell_t22,
                    newdata = bios_t22,
                    include_mahalanobis = TRUE,
                    include_suitability = FALSE)

plot(log(pred_t27), col=rev(viridis::viridis(100)))

# Test 28: Non-numeric input
bios_df_t28 <- bios_df_t23
bios_df_t28$bio1 <- as.character(bios_df_t28$bio1)

pred_t28 <- predict(ell_t22,
                    newdata = bios_df_t28,
                    include_mahalanobis = TRUE)

# Test 29: Wrong object class
pred_t29 <- predict(rng_t22,
                    newdata = bios_t22)

# sample_data() test ------------------------------------------------------

# Test 30: Basic run (SpatRaster, probability, center, random and edge)
library(terra)

wc_t30 <- geodata::worldclim_global(
  var = "bio",
  res = 10,
  path = tempdir()
)

bios_t30 <- wc_t30[[c(1, 12, 15)]]
names(bios_t30) <- c("bio1", "bio12", "bio15")

range_t30 <- data.frame(bio1 = c(-10, 10),
                        bio12 = c(500, 1000),
                        bio15 = c(30, 150))

ell_t30 <- build_ellipsoid(range = range_t30)
plot_nicheR(list(ell_t30))

pred_t30 <- predict(ell_t30,
                    newdata = bios_t30,
                    include_mahalanobis = TRUE,
                    include_suitability = TRUE,
                    suitability_truncated = TRUE,
                    mahalanobis_truncated = TRUE)


# Test B: check with environmental background and compare
pred_t30B <- predict(ell_t30,
                    newdata = bios_t30,
                    suitability_truncated = TRUE,
                    mahalanobis_trucated = FALSE,
                    include_mahalanobis = FALSE,
                    include_suitability = FALSE)

pred_t30B[pred_t30B == 0] <- NA

bios_t30B <- crop(mask(bios_t30, pred_t30B), pred_t30B)
plot(bios_t30B)

bios_t30B_df <- as.data.frame(bios_t30B)
head(bios_t30B_df)

plot_nicheR(list(ell_t30),
            dims = c(1, 2),
            main = "Sampled occurrences in E space",
            occ_list = list(bios_t30B_df),
            occ_col = c("grey"))

# Center
occ_t30C <- sample_data(n_occ = 200,
                       suitable_env = pred_t30,
                       sampling = "center",
                       method = "probability",
                       seed = 42)
# Edge
occ_t30D <- sample_data(n_occ = 200,
                        suitable_env = pred_t30,
                        sampling = "edge",
                        method = "probability",
                        seed = 42)
# Random
occ_t30E <- sample_data(n_occ = 200,
                        suitable_env = pred_t30,
                        sampling = "random",
                        method = "probability",
                        seed = 42)

# TO DO: allow for only e-space return, not only g-space and return with E-space. Whatever prediction came with, ends up with, and preds.
# TO DO: make it is own class object, better for plotting.
# Keep_data = NULL, data.frame = TRUE, raster = FALSE
# TO DO: add virtual data from function from the ellipsoids.
attributes(occ_t30E)
head(occ_t30E)

occ_t30C_bios <- terra::extract(bios_t30, occ_t30C)
occ_t30D_bios <- terra::extract(bios_t30, occ_t30D)
occ_t30E_bios <- terra::extract(bios_t30, occ_t30E)

plot_nicheR(list(ell_t30),
            dims = c(1, 2),
            main = "Sampled occurrences Center",
            occ_list = list(occ_t30C_bios),
            occ_col = "green4")

plot_nicheR(list(ell_t30),
            dims = c(1, 2),
            main = "Sampled occurrences Edge",
            occ_list = list(occ_t30D_bios),
            occ_col = "blue")

plot_nicheR(list(ell_t30),
            dims = c(1, 2),
            main = "Sampled occurrences Random",
            occ_list = list(occ_t30E_bios),
            occ_col = "orange")

# Test 32: Wrong input (missing suitability and suitability_trunc)
pred_t32 <- predict(ell_t30,
                    newdata = bios_t30,
                    include_mahalanobis = TRUE,
                    include_suitability = FALSE,
                    suitability_truncated = FALSE,
                    mahalanobis_truncated = TRUE)

occ_t32 <- sample_data(n_occ = 50,
                       suitable_env = pred_t32,
                       sampling = "center",
                       method = "probability",
                       seed = 42)


# TO DO: prediction = NULL, if null make virtual ell data

# Test 33: Outside sampling of the niche
occ_t33 <- sample_data(n_occ = 10000,
                        suitable_env = pred_t30$suitability,
                        sampling = "center",
                        method = "probability",
                        seed = 42)

occ_t33_bios <- terra::extract(bios_t30, occ_t33)

plot_nicheR(list(ell_t30),
            dims = c(1, 2),
            main = "Sampled occurrences Center",
            occ_list = occ_t33_bios,
            occ_col = "green4")




# apply_bias() tests -------------------------------------------------------

# Test 34: Basic run, single bias layer (terra example raster)
library(terra)

elev_t34 <- terra::rast(system.file("ex/elev.tif", package = "terra"))
names(elev_t34) <- "elev"

bias_t34 <- apply_bias(bias_surface = elev_t34,
                       out_bias = "both",
                       truncated = TRUE,
                       verbose = TRUE)

bias_t34$combination_formula
bias_t34$pooled_bias_sp
bias_t34$directional_bias_stack

plot(bias_t34$pooled_bias_sp, main = "T34 pooled bias (elev)")


# Test 35: Two layers, invert one (bias_dir = c(1, -1))
slope_t35 <- terra::terrain(elev_t34, v = "slope", unit = "degrees")
names(slope_t35) <- "slope"

bias_stack_t35 <- c(elev_t34, slope_t35)

bias_t35 <- apply_bias(bias_surface = bias_stack_t35,
                       bias_dir = c(1, -1),
                       out_bias = "both",
                       truncated = TRUE,
                       verbose = TRUE)

bias_t35$combination_formula
plot(bias_t35$directional_bias_stack, main = "T35 standardized directional layers")
plot(bias_t35$pooled_bias_sp, main = "T35 pooled bias (elev * (1-slope))")


# Test 36: Use suitable_env as mask/template (crop + align)
#create a smaller mask to force cropping
mask_t36 <- terra::crop(elev_t34, terra::ext(5.5, 6, 49.9, 50))
names(mask_t36) <- "mask_template"

bias_t36 <- apply_bias(bias_surface = bias_stack_t35,
                       bias_dir = c(1, 1),
                       suitable_env = mask_t36,
                       out_bias = "both",
                       truncated = TRUE,
                       verbose = TRUE)

plot(bias_t36$pooled_bias_sp, main = "T36 pooled bias masked to suitable_env")


# Test 37: out_bias = "standardized" (no pooled)
bias_t37 <- apply_bias(bias_surface = bias_stack_t31,
                       bias_dir = c(1, 1),
                       out_bias = "standardized",
                       truncated = TRUE,
                       verbose = TRUE)

bias_t37$pooled_bias_sp
bias_t37$directional_bias_stack
plot(bias_t37$directional_bias_stack, main = "T37 standardized only")


# Test 38: Wrong bias_dir input (should error)
bias_t38 <- apply_bias(bias_surface = bias_stack_t31,
                       bias_dir = c(1, 0),
                       out_bias = "both",
                       truncated = TRUE,
                       verbose = FALSE)


# Test 39: truncated = FALSE (no cropping step)
bias_t39 <- apply_bias(bias_surface = bias_stack_t31,
                       bias_dir = c(1, 1),
                       suitable_env = mask_t32,
                       out_bias = "both",
                       truncated = FALSE,
                       verbose = TRUE)

plot(bias_t39$pooled_bias_sp, main = "T39 pooled bias (no crop, no masked)")

