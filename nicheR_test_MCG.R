devtools::load_all()



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
                     effect_direction = "direct", include_processed_layers = TRUE)
# To do: message says that it is a stack but it only has one, so remove that message

bias$combination_formula
bias$composite_surface

terra::plot(bias$composite_surface, main = "T34 pooled bias (elev)")



# Build Ellipsoid ---------------------------------------------------------

ell <- build_ellipsoid(range = range)
names(ell)
ell

plot_ellipsoid(ell)
plot_ellipsoid_pairs(ell)
# to do: plot_ellipsoid_grid or plot_all_ellipsoids

# plot_ellipsoid(ell,ell_prededict# plot_ellipsoid(ell, background = bios)
# to do: need to add checks, it does not stop when backgorung is a spatRaster

plot_ellipsoid(ell, background = bios_df[, c(3,5)], dim = c(1,3))
# to do: set the dimensions by name with the backgorund too

plot_ellipsoid(ell, background = bios_df[, c(3,4)])


# Predict Suitability -----------------------------------------------------

# VIRTUAL DATA
ell_predict_v <- predict(ell) #this will generate virtual data of 1000 by default

# to do: remove the option for creating virtual data internally, require the user to provide the newdata
newdata <- as.data.frame(virtual_data(object = ell,
                                      n = 10000,
                                      truncate = FALSE))
head(newdata)
head(ell_predict_v)

plot_ellipsoid(ell, background = ell_predict_v[,1:2])
# to do: match the same dimension used for the ellipsoid

# RASTER
ell_predict_r <- predict(ell, newdata = bios)
ell_predict_r #default sutibility and mahalanobis, no truncation

terra::plot(ell_predict_r)

# DATA FRAME
ell_predict_df <- predict(ell, newdata = bios_df)
head(ell_predict_df) #default sutibility and mahalanobis, no truncation

# For DF
set.seed(123)
plot_ellipsoid(ell)
# A different way to add the background
add_data(x = ell_predict_df$bio1, y = ell_predict_df$bio12,
         pts_sample = 8000, pts_col = "grey", pch = 20)
add_data(x = ell_predict_df$bio1, y = ell_predict_df$bio12,
         col_layer = ell_predict_df$Mahalanobis,
         rev_col = TRUE, pts_sample = 8000, pch = 16)

# For virtual data
set.seed(123)
plot_ellipsoid(ell)
# A different way to add the background
add_data(x = ell_predict_v$bio1, y = ell_predict_v$bio12,
         pts_col = "grey", pch = 20)
add_data(x = ell_predict_v$bio1,
         y = ell_predict_v$bio12,
         col_layer = ell_predict_v$suitability,
         rev_col = FALSE,
         pch = 20)
add_ellipsoid(ell, col_ell = "red", lw = 2)

# to do: add to plot_ellipsoid and argument of prediction = NULL, and if prediction is there color layer can be other than NULL


# Apply bias to prediction ------------------------------------------------

biased_predict_r <- apply_bias(prepared_bias = bias,
                               prediction = ell_predict_r)

# To do: only allow prediction to be one layer
# To do: fix space in message "tonpredictions"

# Sample data -------------------------------------------------------------

# UNBIASED SAMPLING - Raster
sample_data_r <- sample_data(n_occ = 100,
                             prediction = ell_predict_r,
                             prediction_layer = "suitability",
                             sampling = "centroid",
                             method = "suitability")
head(sample_data_r)
# To do: tolower(args) for method and sampling

# UNBIASED SAMPLING - Data frame
sample_data_df_cn <- sample_data(n_occ = 100,
                                 prediction = ell_predict_df,
                                 prediction_layer = "suitability",
                                 sampling = "centroid",
                                 method = "suitability")
head(sample_data_df_cn)


plot_ellipsoid(ell, main = "Sampled Center Unbiased")
add_data(x = ell_predict_df$bio1, y = ell_predict_df$bio12,
         pts_col = "grey", pch = 20, pts_sample = 8000)
add_data(x = sample_data_df_cn$bio1, y = sample_data_df_cn$bio12,
         pts_col = "black", pch = 4, lwd = 2)
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
# To do: Make sure this shows a stric removing the zeros outside the ellipsoid
# To do: consider some of the samples can fall outside the ellipsoid, use trucated and stric = TRUE for sampling strictly inside the ellipsoid.

# PLOT FOR CENTER AND EDGE - normal data
par(mfrow = c(2,1))
plot_ellipsoid(ell, main = "Sampled Edge Unbiased",
               background = bios_df[, c(3,4)])
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


