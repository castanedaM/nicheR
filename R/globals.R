utils::globalVariables(c(
  ".", # Often used in dplyr pipelines for the current data
  ".data", # Used by ggplot2 and dplyr for data masking
  "%", # Placeholder for the pipe operator, as it's often re-exported
  "add_trace", # From plotly
  "add_markers", # From plotly
  "aes", # From ggplot2
  "build_ellipsoid", # Your own function, but should be exported or qualified
  "geom_path", # From ggplot2
  "geom_point", # From ggplot2
  "geom_segment", # From ggplot2
  "get_suitable_environment", # Your own function, but should be exported or qualified
  # If FN_1 is an internal object, it should be treated differently.
  # If FN_1 is meant to be a dynamically generated object, its direct use
  # in plot_e_space might need reconsideration or a different approach
  # for R CMD check. Assuming FN_1 is an internal dataset/object, it should
  # be listed here.
  "FN_1",
  "%>%"
))
