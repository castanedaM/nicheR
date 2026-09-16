# Title: Instructions text for nicheR Shiny app
# Description: Centralized instructional copy shown above inputs (p() blocks)
#              and tooltip text shown via icon("circle-info"). Organized by
#              the tab and section each instruction belongs to. Keys are
#              prefixed with the tab they belong to, matching the input and
#              output naming convention used in the server scripts.
# Dependencies: depends on MAX_DIMS from global.R
# Date last updated: 08/18/2026

instructions <- list(

  # ABOUT TAB --------------------------------------------------------------

  about_app = HTML("An R package for building ellipsoid-based ecological virtual niches.\nBuild niche from environmental ranges, predict suitable
areas, account for sampling bias, and generate virtual data, all in
a single reproducible workflow."),

  about_build = HTML("Define a species niche as an ellipsoid in environmental
space using background layers and user-defined variable ranges."),

  about_build_points = c(
    "Set ranges manually, from occurrence data, or from background statistics",
    "Adjust covariance to rotate the ellipsoid",
    "Move the niche centroid without changing axis tolerance ranges",
    "Save multiple named ellipsoid versions"
  ),

  about_predict = HTML("Project the ellipsoid onto raster layers or a data
frame to produce continuous or binary suitability surfaces."),

  about_predict_points = c(
    "Output: suitability, Mahalanobis distance, and truncated versions",
    "Binarization uses the ellipsoid's own confidence level",
    "E-space and G-space plots included",
    "Batch predict across multiple saved versions"
  ),

  about_bias = HTML("Upload a sampling bias raster to weight occurrence
generation toward areas with specific detection effort."),

  about_bias_points = c(
    "Optional but important for data-limited species",
    "Accepts any raster matching the background extent",
    "Higher cell values increase occurrence probability",
    "Examples: urbanization, distance to water, road density, collector coverage, any detection proxy"
  ),

  about_generate = HTML("Sample virtual data from the fitted niche,
optionally weighted by the bias layer."),

  about_generate_points = c(
    "Specify number of presences and background ratio",
    "Output is a data frame of coordinates and environmental values",
    "Designed for virtual species and simulation workflows",
    "Supports rare-species SDM validation studies"
  ),


  # BUILD TAB: INPUTS ------------------------------------------------------

  build_data_input_type = HTML("Choose how to provide your background
environmental data."),

  build_data_upload = HTML("Background data: your study area's environmental
conditions.<br>
Upload either a raster (.tif, .rds) or a CSV of the same data. Once one is
provided the other is hidden. Use Clear files to switch."),

  build_raster_file_tooltip = "Environmental conditions of the study area in raster format.
Provides both environmental and geographic space. Accepted: .tif, .rds",

  build_df_file_tooltip = "The same data in tabular form. Geographic space is
available only if the file has coordinate columns on a regular grid. Accepted: .csv, .rds",

  build_no_spatial_cols = HTML("No coordinate columns were found. The app
looks for x, lon, longitude, or easting paired with y, lat, latitude, or
northing. You can still build and inspect the ellipsoid in environmental
space, but geographic maps will not be available."),

  build_irregular_grid = HTML("Coordinate columns were found but do not form a
regular grid, so they could not be converted to a raster. Work will continue
in environmental space only."),

  build_prev_session = HTML("Upload a session file (.rds) saved from a previous
visit. The app will resume at whatever step you last left off, including
any confirmed variables and built ellipsoid if they were saved."),

  build_virtual_mode = HTML("Work entirely in environmental space (E-space),
with no geographic coordinates. You will name your variables and define their
ranges directly, without uploading any spatial data."),

  build_example_data = HTML("Use the example dataset bundled with nicheR:
WorldClim bioclimatic variables for Central America, including a bias layer
for later steps. Click Continue to select which variables to use."),


  # BUILD TAB: VARIABLES ---------------------------------------------------

  build_variable_settings = HTML(paste0("Select up to ", MAX_DIMS,
                                        " variables to use in your analysis.<br>
Uncheck a variable to exclude it. The plot beside you shows the background
distribution of your selected variables and updates as you choose.")),

  build_virtual_variables = HTML(paste0("Choose how many variables to work
with (up to ", MAX_DIMS, ") and give each one a name. Default names (var1,
var2, ...) are used if left blank.")),

  build_edit_variables_tooltip = "Editing your variables will delete the
current selection, the ellipsoid, and any covariance adjustments.",


  # BUILD TAB: RANGES ------------------------------------------------------

  build_range_choice = HTML("Ranges set the minimum and maximum value of each
variable, which together define the extent of the ellipsoid. Choose one of
the three methods below. You can switch methods and change values at any
time, then rebuild."),

  build_range_manual = HTML("Define the minimum and maximum values for each
variable.<br>
Defaults are the first and third quartiles of your background data. Lines in
the plot update as you change these values. Use Reset to Defaults to restore
them, then Initialize Ellipsoid when ready."),

  build_range_data = HTML("Upload a CSV with occurrence or other environmental
data for your species of interest.\nColumn names must match your selected variables exactly, or the upload will not be recognized. Once matched, observed min/max values appear below, and
Expand Min/Max (%) let you widen the range outward in either direction. Lines
in the plot update as you adjust them."),

  build_range_stats = HTML("Provide summary statistics for the chosen
variables. Defaults shown here are derived from your background data, but can be changed
to known values for your species of interest. Expand Min/Max (%) widen the
resulting range outward, not the mean itself. Confidence level controls how
wide the range is before any expansion."),

  build_range_manual_tooltip = "Type minimum and maximum values directly for
each variable.",

  build_range_data_tooltip = "Upload a CSV to derive observed min/max ranges,
with optional expansion.",

  build_range_stats_tooltip = "Derive ranges from mean and standard deviation,
either from your background data or entered manually.",

  build_expand_range_tooltip = "Percentage to widen the range outward from the
observed or computed bounds. Can be negative to shrink the range inward.",

  build_range_file_unreadable = "Could not read the uploaded file. Check that
it is a valid .csv or .rds.",

  build_range_file_mismatch = "The uploaded file does not contain all of your
selected variables. Column names must match exactly.",

  build_range_invalid = "Please check your range inputs before building. Every
variable needs a valid minimum below its maximum.",

  build_cl_ell_tooltip = paste0(
    "Sets the chi-square cutoff that defines the boundary of the ellipsoid. ",
    "A higher value draws a larger ellipsoid around the same centroid. ",
    "It does not change the centroid, the covariances, or the ranges, so it ",
    "applies right away and does not need a rebuild. Default is 0.95."
  ),

  build_cl_interval_tooltip = paste0(
    "Used only by From Stats. It is the normal interval around each mean ",
    "that becomes the minimum and maximum, so 0.95 gives mean +/- 1.96 SD. ",
    "This is a range input, not the ellipsoid boundary, so changing it ",
    "means rebuilding the ellipsoid. Use 0.9973 if you want the ellipsoid's ",
    "marginal SD to match the SD you typed."
  ),

  build_export_ell = paste0(
    "Saves this ellipsoid to an .rds file, including its ranges, covariances, ",
    "centroid, and confidence level. Read it back with read_nicheR(). The ",
    "name starts as the internal ID plus today's date so the file can be ",
    "traced back to the library, but you can change it."
  ),

  # BUILD TAB: COVARIANCE --------------------------------------------------

  build_covariance = HTML("Covariances describe the relationship between pairs
of variables.\nDrag a slider to rotate the ellipsoid along that pair. Moving one slider updates the valid covariance limits of the others. If a combination would make the ellipsoid invalid, an error is shown and that slider resets. Click the rotate-left icon next to a
slider to reset that pair back to zero. Click 'Set Covariances' when ready."),

  build_covariance_reset_tooltip = "Reset this pair's covariance back to zero.",

  build_cov_set = "Covariances have been set for this ellipsoid. If you have not saved this elliposid editing it will reset them back to zero or to those values of their root",


  # BUILD TAB: CENTROID ----------------------------------------------------

  build_centroid_mover = HTML("Move the ellipsoid through environmental space
without changing its shape or size. Each slider shifts the centroid along one
variable, and the ellipsoid follows. Use Reset all to return to the centroid
this version started from. Click 'Set centroid' when ready."),

  build_centroid_set = "Centroid has been set for this ellipsoid. If you have not saved this elliposid editing it will reset them back to original values or those of their root",

  # BUILD TAB: LIBRARY -----------------------------------------------------

  build_library = HTML("Every ellipsoid you save is listed here. The
highlighted row at the top is the working ellipsoid, the one used by the
plots and by every tab downstream. Each saved version can be viewed
read-only, edited in place, copied into a new version, or deleted."),

  build_library_empty = "No saved versions yet. Save the working ellipsoid to
add it here.",

  build_view_only = "You are viewing a saved version. Editing is locked. Use the 'Edit'
button to go back to editing, or the copy button in the library, to make changes.",

  build_save_ell_modal = "Give this ellipsoid a name. Use letters, numbers,
and spaces only. Spaces will be replaced with underscores.",

  build_delete_ell = "This will permanently remove the ellipsoid and any
prediction results associated with it.",

  build_reference = HTML("The ellipsoid summary reports volume and centroid
changes relative to a reference. By default this is whatever the working
ellipsoid was built or copied from. Point it at any saved version to compare
against that instead. Changing it only affects the summary, it does not move
your ellipsoid."),

  build_session_loaded = "Session loaded successfully.",

  build_session_invalid = "That file is not a nicheR session file.",

  build_session_old_version = HTML("This session was saved by an older
version of the app. It has been loaded, but bias layers may be missing and
ellipsoid ids may need checking."),

  build_range_locked = paste0(
    "These are the ranges this ellipsoid was built from. To use different ",
    "ranges, start over."
  ),

  build_start_over = paste0(
    "This clears the working ellipsoid and takes you back to choosing a ",
    "range method."
  ),

  # PREDICT TAB ------------------------------------------------------------

  predict_adjust_trunc_tooltip = "Adjust the level of truncation within the
current ellipsoid. This truncates the prediction inwards.",

  predict_library = HTML("Ellipsoids saved on the Build tab appear here.
Select one to view its settings read-only, or delete it along with any
predictions made from it. To edit an ellipsoid, go back to the Build tab."),

  predict_library_empty = "No saved ellipsoids yet. Save one on the Build tab
to predict with it.",

  predict_delete_ell = "This will permanently remove the ellipsoid, its
predictions, and any biased predictions derived from them.",

  predict_ellipsoid_select_tooltip = "Choose which saved ellipsoid to predict
with, or All versions to predict with every one at once.",

  predict_virtual_unavailable = HTML("Prediction needs raster layers, so it
is unavailable in virtual mode. Go straight to Generate, which samples
directly from the ellipsoid."),

  bias_virtual_unavailable = HTML("Sampling bias is geographic, so it is
unavailable in virtual mode."),
  # BIAS TAB ---------------------------------------------------------------

  bias = HTML("Bias input adds controlled sampling bias to a prediction
layer. Note that once bias has been applied, the prediction is no longer a
probability.<br>
In the steps below you will need one or more raster layers (.rds, .tif).
These are standardized to a 0-1 scale and can have a direct effect (increase
sampling probability) or an inverse effect (decrease sampling probability)."),

  bias_input = HTML("Upload one or more raster layers (.rds, .tif) to
represent sampling bias across the study area, or use the example layers
provided."),

  bias_skipped = "Bias skipped. Occurrences will be sampled from the
unbiased prediction.",

  bias_needs_prediction = "Run a prediction on the Predict tab before adding
sampling bias.",

  bias_prepare = HTML("Assign a direction of effect to each bias layer, then
combine them into a single composite surface. Direct means higher values
increase sampling probability, inverse means they decrease it."),

  bias_apply = HTML("Multiply a prediction by the composite bias surface.
Combinations that already exist are skipped, so you can add layers without
losing what is already applied."),

  bias_edit_upload = "This removes the current bias layers, the prepared
surface, and any applied bias. You will need to upload and prepare again.",

  bias_edit_prepare = "This removes the prepared surface and any applied
bias. The uploaded layers are kept.",

  bias_clear_apply = "This removes every biased surface. The uploaded layers
and the prepared surface are kept.",

  bias_library = HTML("Ellipsoids saved on the Build tab appear here, with
whether each has been predicted and biased. To edit an ellipsoid, go back to
the Build tab."),

  bias_library_empty = "No saved ellipsoids yet.",

  bias_delete_ell = "This will permanently remove the ellipsoid, its
predictions, its biased predictions, and any occurrences generated from them.",

  bias_mask_na_tooltip = "Union keeps any pixel with at least one valid value
across layers.
Intersection keeps only pixels valid in every layer.",

  bias_layer_stats_tooltip = "Mean and standard deviation of non-NA raster values.",

  bias_direction_tooltip = "Direct: higher values increase sampling probability.
Inverse: higher values decrease sampling probability.",

  bias_layer_tooltip = "Layer must be a suitability surface with values in [0, 1].",

  bias_ellipsoid_select_tooltip = "Choose which ellipsoid's prediction to
apply bias to, or All versions to apply to every one.",

  bias_apply_direction_tooltip = paste0(
    "Direct samples toward high prediction values, inverse away from them, ",
    "using (max + min) - x. Changes the prediction, not the bias surface."
  ),

  # GENERATE TAB -----------------------------------------------------------

  generate_needs_prediction = "Run a prediction on the Predict tab before
generating occurrences.",

  generate_intro = HTML("Sample virtual data from a prediction surface.
                        Biased layers are shown in orange if bias has been applied."),

  generate_no_surface = "Select at least one prediction layer before generating.",

  generate_mask = "If no mask is provided, sampling covers the full
prediction extent.",

  generate_library = HTML("Ellipsoids saved on the Build tab appear here,
with how many records each has generated. To edit an ellipsoid, go back
to the Build tab. To view the generated set of records press on view (eye)"),

  generate_library_empty = "No saved ellipsoids yet.",

  generate_delete_ell = "This will permanently remove the ellipsoid, its
predictions, its biased predictions, and any occurrences generated from them.",

  generate_sampling_tooltip = "Direct: higher probability near the niche centroid.
  Inverse: higher probability near the niche edge.
  Uniform: random uniform probability across all suitable cells.",

  generate_surface_tooltip = "Select one or more prediction layers to sample
from. Orange layers have bias applied.",

  generate_strict_tooltip = "When True, removes NA and zero-valued cells
before sampling. Recommended for truncated layers.",

  generate_advanced_tooltip = "Optional sampling mask and the random seed
used for reproducible draws.",

  generate_seed_tooltip = "Same seed and same settings produce the same
occurrences. Change it to draw a different replicate.",

  generate_ellipsoid_select_tooltip = "Choose which ellipsoid to generate
data from, or 'All versions' to generate for every one.",

  generate_summary = HTML("Every occurrence set you generate is listed here,
grouped by ellipsoid. Sets accumulate: generating with a new seed adds a set
rather than replacing one. Download individual sets, everything for one
ellipsoid, or all sets at once."),

  generate_summary_empty = "No occurrence sets yet. Generate one to see it here.",

  generate_max_visible = "Up to four sets can be shown at once. Hide one
before showing another.",

  generate_summary = HTML("Every occurrence set you generate is listed here.
Sets accumulate: changing any parameter adds a set rather than replacing one.
Use the eye to choose which appear in the plots, up to four at a time."),

  generate_virtual_intro = HTML("Virtual mode samples directly from the
ellipsoid's multivariate normal distribution. There is no prediction surface
and no geography, so the points are environmental values only."),

  generate_needs_ellipsoid = "Build an ellipsoid before generating points.",

  generate_truncate_tooltip = "When true, points are constrained inside the
ellipsoid's confidence limit. When false, they follow the unbounded normal
distribution and some will fall outside.",

  generate_effect_tooltip = "Direct concentrates points near the centroid.
Inverse pushes them toward the edges. Uniform spreads them evenly through
the ellipsoid volume.",

  generate_effect_needs_truncate = "Inverse and uniform require truncation.",

  # SHARED: PLOT SETTINGS --------------------------------------------------

  plot_settings = HTML("Adjust point shape, size, and colors, toggle range
lines, the ellipsoid, and centroid on or off, and zoom in manually or to
the ellipsoid's extent."),

  plot_zoom_tooltip = "Auto fits the plot to your data and ellipsoid.
Zoom to ellipsoid frames the plot tightly around the ellipsoid boundary.
Manual lets you set exact axis limits.",


  # SHARED: SESSION --------------------------------------------------------

  session_save = HTML("Save your current progress as an .rds file. This
includes confirmed variables and the built ellipsoid, if any. Load it later
from the Previous Session input option to resume where you left off.<br>
If no ellipsoid has been confirmed yet, only the most recent in-progress
values are saved.")

)
