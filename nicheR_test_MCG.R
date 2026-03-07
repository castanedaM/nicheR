
# Data collatetion and preparation ----------------------------------------
wc <- geodata::worldclim_global(var = "bio",
                                res = 10,
                                path = tempdir())

bios <- wc[[c(1, 12, 15)]]
names(bios) <- c("bio1", "bio12", "bio15")

bios_df <- as.data.frame(bios, xy = TRUE)

range <- data.frame(bio1 = c(-5, 10),
                    bio12 = c(500, 750),
                    bio15 = c(30, 150))

# Bias layers
bias <- prepare_bias(bias_surface = wc[[c(6, 11)]],
                     effect_direction = "direct",
                     include_processed_layers = TRUE)


bias$combination_formula
bias$composite_surface

terra::plot(bias$composite_surface)

# Build Ellipsoid ---------------------------------------------------------

ell <- build_ellipsoid(range = range)
names(ell)
ell

plot_ellipsoid(ell)
plot_ellipsoid_pairs(ell)
# to do: plot_ellipsoid_grid or plot_all_ellipsoids

plot_ellipsoid(ell, background = bios) # has to fail

# Predict Suitability -----------------------------------------------------

bios_df_v <- virtual_data(ell, n = 1000, truncate = FALSE)
head(bios_df_v)

ell_predict_v <- predict(ell) #should fail
ell_predict_v <- predict(ell, bios_df_v)

plot_ellipsoid(ell, background = ell_predict_v)

# RASTER
ell_predict_r <- predict(ell, newdata = bios,
                         suitability_truncated   = TRUE,
                         mahalanobis_truncated   = TRUE,
                         keep_data = TRUE)
ell_predict_r #default sutibility and mahalanobis, no truncation

terra::plot(ell_predict_r)


# DATA FRAME
ell_predict_df <- predict(ell, newdata = bios_df,
                          suitability_truncated   = TRUE)
head(ell_predict_df) # default suitability and mahalanobis, no truncation

ell_predict_df_r <- as.data.frame(ell_predict_r, xy = TRUE)


# For DF
set.seed(123)

# No TRUNCATED ZOOMED OUT
plot_ellipsoid(ell,
               background = ell_predict_df,
               pch = 20,
               bg_sample = 10000)
add_data(data = ell_predict_df, x = "bio1", y = "bio12",
         col_layer = "suitability",
         rev_pal = T,
         pal = terrain.colors(100),
         bg_sample = 10000, pch = 20)
add_ellipsoid(ell, col_ell = "red")

# TRUCATED ZOOMED IN
plot_ellipsoid(ell)
add_data(data = ell_predict_df, x = "bio1", y = "bio12",
         col_layer = "suitability_trunc", rev_pal = T,
         bg_sample = 50000, pch = 20)
add_ellipsoid(ell, col_ell = "red", lwd = 2)


# Apply bias to prediction ------------------------------------------------

biased_predict_r <- apply_bias(prepared_bias = bias,
                               prediction = ell_predict_r,
                               prediction_layer = "suitability")

# Sample data -------------------------------------------------------------

# UNBIASED SAMPLING - Raster
sample_data_r <- sample_data(n_occ = 100,
                             prediction = ell_predict_r,
                             prediction_layer = "suitability",
                             sampling = "centroid",
                             method = "suitability")
head(sample_data_r)
# To do: check for suitability prediction layer with mahalanobis method

# UNBIASED SAMPLING - Data frame
sample_data_df_cn <- sample_data(n_occ = 100,
                                 prediction = ell_predict_df,
                                 prediction_layer = "suitability",
                                 sampling = "centroid",
                                 method = "mahalanobis")
head(sample_data_df_cn)


plot_ellipsoid(ell, main = "Sampled Center Unbiased")
add_data(data = ell_predict_df,
         x = "bio1", y = "bio12",
         pts_col = "grey",
         pch = 20,
         bg_sample = 8000)
add_data(data = sample_data_df_cn,
         x = "bio1", y = "bio12",
         pts_col = "black",
         pch = 4,
         lwd = 2)
add_ellipsoid(ell, col_ell = "red", lwd = 2)


# BIASED SAMPLING - Raster
sample_data_b <-
  sample_biased_data(n_occ = 100,
                     prediction = biased_predict_r,
                     prediction_layer = "suitability_biased_direct")

head(sample_data_b)


# VIRTUAL SAMPLING - virtual
sample_data_v_c <- sample_virtual_data(n_occ = 100,
                                       virtual_prediction = ell_predict_v,
                                       prediction_layer = "suitability",
                                       sampling = "centroid")

sample_data_v_e <- sample_virtual_data(n_occ = 100,
                                       virtual_prediction = ell_predict_v,
                                       prediction_layer = "suitability",
                                       sampling = "edge")

head(sample_data_v_c)


plot_ellipsoid(ell, main = "Sampled Center Unbiased")
add_data(x = ell_predict_v$bio1, y = ell_predict_v$bio12,
         pts_col = "grey", pch = 20, pts_sample = 8000)
add_data(x = sample_data_v_c$bio1, y = sample_data_v_c$bio12,
         pts_col = "black", pch = 4, lwd = 2)
add_ellipsoid(ell, col_ell = "red", lwd = 2)



# UNBIASED SAMPLING - Data frame - Edge
sample_data_df_ed <- sample_data(n_occ = 100,
                                 prediction = ell_predict_df,
                                 prediction_layer = "suitability_trunc",
                                 sampling = "edge",
                                 method = "suitability")

plot_ellipsoid(ell,
               main = "Sampled Center Unbiased",
               # background = bios_df,
               # bg_sample = 10000,
               pch = 20)
add_data(data = ell_predict_df,
         x = "bio1", y = "bio12",
         pts_col = "grey",
         pch = 20,
         bg_sample = 8000)
add_data(data = sample_data_df_ed,
         x = "bio1", y = "bio12",
         pts_col = "black",
         pch = 4,
         lwd = 2)
add_ellipsoid(ell, col_ell = "red", lwd = 2)


# PLOT FOR CENTER AND EDGE - normal data
# par(mfrow = c(2,1))
plot_ellipsoid(ell, main = "Sampled Edge Unbiased",
               background = bios_df)
add_data(x = ell_predict_df$bio1, y = ell_predict_df$bio12,
         pts_col = "grey", pch = 20, pts_sample = 8000)
add_data(x = sample_data_df_ed$bio1, y = sample_data_df_ed$bio12,
         pts_col = "black", pch = 4, lwd = 2)
add_ellipsoid(ell, col_ell = "red", lwd = 2)


plot_ellipsoid(ell, main = "Sampled Center Unbiased",
               background = bios_df[, c(3,4)])
add_data(x = ell_predict_df$bio1, y = ell_predict_df$bio12,
         pts_col = "grey", pch = 20, pts_sample = 8000)
add_data(x = sample_data_df_cn$bio1, y = sample_data_df_cn$bio12,
         pts_col = "black", pch = 4, lwd = 2)
add_ellipsoid(ell, col_ell = "red", lwd = 2)


# PLOT FOR - VIRTUAL DATA
par(mfrow = c(2,1))
plot_ellipsoid(ell, main = "Sampled Center Unbiased")
add_data(x = ell_predict_v$bio1, y = ell_predict_v$bio12,
         pts_col = "grey", pch = 20, pts_sample = 8000)
add_data(x = sample_data_v_c$bio1, y = sample_data_v_c$bio12,
         pts_col = "black", pch = 4, lwd = 2)
add_ellipsoid(ell, col_ell = "red", lwd = 2)

plot_ellipsoid(ell, main = "Sampled Edge Unbiased")
add_data(x = ell_predict_v$bio1, y = ell_predict_v$bio12,
         pts_col = "grey", pch = 20, pts_sample = 8000)
add_data(x = sample_data_v_e$bio1, y = sample_data_v_e$bio12,
         pts_col = "black", pch = 4, lwd = 2)
add_ellipsoid(ell, col_ell = "red", lwd = 2)


# To do: do an example with uniform virtual data


