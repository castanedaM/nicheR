# Title: Server for shiny nicheR
# Description: The server of the app
# Last Updated: 10/06/2026


function(input, output, session){

  session_data <- reactiveValues(

    # Input
    input_mode = NULL,
    file_type = NULL,
    data_source = NULL,

    bg_raster = NULL,
    bg_df = NULL,

    vars = NULL,

    # Ranges
    session_range = NULL,
    df_range = NULL,

    # Ellipsoids
    ellipsoid_list = list(),
    current_ellipsoid = NULL,
    pending_ell_delete = NULL,

    # Prediction
    ellipsoid_prediction_list = list(),
    prediction_settings = list(),

    # Bias
    bias_source = NULL,
    bias_raster = NULL,
    prepared_bias = NULL,
    ellipsoid_prediction_list_biased = list(),
    bias_settings = list(),

    # Generate
    sampling_mask = NULL,
    ellipsoid_point_sets_list = list()

  )

  source("server/helpers.R", local = TRUE)
  source("server/save_load_session.R", local = TRUE)
  source("server/report_main.R", local = TRUE)
  source("server/report_builder.R", local = TRUE)

  source("server/data_tab.R", local = TRUE)
  source("server/build_tab.R", local = TRUE)
  source("server/build_tab_plot.R", local = TRUE)

  source("server/predict_tab.R", local = TRUE)
  source("server/predict_tab_plot.R", local = TRUE)

  source("server/bias_tab.R", local = TRUE)
  source("server/bias_tab_plot.R", local = TRUE)

  source("server/generate_tab.R", local = TRUE)
  source("server/generate_tab_plot.R", local = TRUE)

}
