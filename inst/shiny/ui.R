# Title: UI for shiny nicheR
# Description: The UI of the app
# Last Updated: 08/20/2026

dashboardPage(
  dashboardHeader(
    title = tags$span("nicheR", class = "text-app-title"),
    titleWidth = 200,
    tags$li(
      class = "dropdown",
      downloadButton("save_session_btn",
                     label = "Save",
                     icon = icon("floppy-disk"), class = "btn-primary",
                     title = "Save current session as .rds file")
    ),
    tags$li(
      class = "dropdown",
      downloadButton("create_report",
                     label = "Report",
                     icon = icon("file-lines"), class = "btn-primary",
                     title = "Save current session as an Rmarkdown with rendered HTML.")
    )
  ),

  dashboardSidebar(
    width = 200,
    sidebarMenu(
      id = "sidebar_menu",
      menuItem("About",
               tabName = "about",
               icon = icon("circle-info", style = "margin-right: 4px")),
      menuItem("1. Build ellipsoid",
               tabName = "build_tab",
               icon = icon("gear", style = "margin-right: 4px")),
      menuItem("2. Prediction",
               tabName = "predict_tab",
               icon = icon("angles-right", style = "margin-right: 4px")),
      menuItem("3. Bias",
               tabName = "bias_tab",
               icon = icon("table", style = "margin-right: 4px")),
      menuItem("4. Generate data",
               tabName = "generate_tab",
               icon = icon("eye-dropper", style = "margin-right: 4px"))
    ),

    tags$div(
      style = "position: absolute; bottom: 50px; width: 100%; text-align: center;",
      img(src = "nicheR_lg.png",
          style = "max-width: 140px;")
    ),

    tags$div(
      style = "position: absolute; bottom: 15px; width: 100%; text-align: center;",
      tags$a(
        href   = "https://github.com/castanedaM/nicheR",
        target = "_blank",
        style  = "color: #aaa; font-size: 0.8em;",
        icon("github"), " GitHub"
      )
    )
  ),

  dashboardBody(
    useShinyjs(),
    includeCSS("www/styles.css"),
    tabItems(
      tabItem(
        tabName = "about",

        fluidRow(
          column(width = 12,
                 tags$div(tags$span("About nicheR", class = "text-app-title"),
                          br(), br(),
                          tags$span(instructions$about_app,
                                    class = "text-about")
                 ),
                 br()
          )
        ),

        fluidRow(
          column(width = 12,
                 p(class = "text-about",
                   "For a walkthrough of the app itself, see the ",
                   tags$a(href   = "https://castanedam.github.io/nicheR/articles/shiny_app_vignette.html",
                          target = "_blank", rel = "noopener",
                          icon("book"), " Shiny app vignette"),
                   "."),
                 p(class = "text-about",
                   "Each step of the app corresponds to a set of functions in the ", tags$em("nicheR"),
                   " package. To understand what runs underneath, see the vignettes on ",
                   tags$a(href = "https://castanedam.github.io/nicheR/articles/creating_ellipsoid_based_niches.html",
                          target = "_blank", rel = "noopener", "building ellipsoids"), ", ",
                   tags$a(href = "https://castanedam.github.io/nicheR/articles/predict.html",
                          target = "_blank", rel = "noopener", "making predictions"), ", ",
                   tags$a(href = "https://castanedam.github.io/nicheR/articles/bias.html",
                          target = "_blank", rel = "noopener", "sampling bias"), ", and ",
                   tags$a(href = "https://castanedam.github.io/nicheR/articles/generating_occurrence.html",
                          target = "_blank", rel = "noopener", "generating data"), ".")
          )
        ),

        fluidRow(
          column(width = 2),

          column(width = 4,
                 box(title = tagList(icon("gear"),
                                     tags$span("1. Build ellipsoid",
                                               class = "text-widget-title")),
                     width = 12,
                     collapsible = TRUE,
                     collapsed = TRUE,
                     solidHeader = TRUE,
                     status = "primary",

                     p(instructions$about_build, class = "text-about"),
                     tags$ul(
                       lapply(instructions$about_build_points, function(x) {
                         tags$li(x, class = "text-about")
                       })
                     ),

                     tags$a(href   = "https://castanedam.github.io/nicheR/articles/creating_ellipsoid_based_niches.html",
                            target = "_blank",
                            icon("book"), " Build vignette")
                 )
          ),

          column(width = 4,
                 box(title = tagList(icon("angles-right"),
                                     tags$span("2. Prediction",
                                               class = "text-widget-title")),
                     width = 12,
                     collapsible = TRUE,
                     collapsed = TRUE,
                     solidHeader = TRUE,
                     status = "primary",

                     p(instructions$about_predict, class = "text-about"),
                     tags$ul(
                       lapply(instructions$about_predict_points, function(x) {
                         tags$li(x, class = "text-about")
                       })
                     ),

                     tags$a(href   = "https://castanedam.github.io/nicheR/articles/predict.html",
                            target = "_blank",
                            icon("book"), " Predict vignette")
                 )
          ),

          column(width = 2)
        ),

        br(),

        fluidRow(
          column(width = 2),

          column(width = 4,
                 box(title = tagList(icon("table"),
                                     tags$span("3. Bias",
                                               class = "text-widget-title")),
                     width = 12,
                     collapsible = TRUE,
                     collapsed = TRUE,
                     solidHeader = TRUE,
                     status = "primary",

                     p(instructions$about_bias, class = "text-about"),
                     tags$ul(
                       lapply(instructions$about_bias_points, function(x) {
                         tags$li(x, class = "text-about")
                       })
                     ),

                     tags$a(href   = "https://castanedam.github.io/nicheR/articles/bias.html",
                            target = "_blank",
                            icon("book"), " Bias vignette")
                 )
          ),

          column(width = 4,
                 box(title = tagList(icon("eye-dropper"),
                                     tags$span("4. Generate data",
                                               class = "text-widget-title")),
                     width = 12,
                     collapsible = TRUE,
                     collapsed = TRUE,
                     status = "primary",
                     solidHeader = TRUE,

                     p(instructions$about_generate, class = "text-about"),
                     tags$ul(
                       lapply(instructions$about_generate_points, function(x) {
                         tags$li(x, class = "text-about")
                       })
                     ),

                     tags$a(href   = "https://castanedam.github.io/nicheR/articles/generating_occurrence.html",
                            target = "_blank",
                            icon("book"), " Generate vignette")
                 )
          ),
          column(width = 2),
        ),

        br(),

        fluidRow(
          column(width = 12,
                 p(class = "text-about",
                   "This app allows for some basic visualization. For more specialised and customised visualization, see the ",
                   tags$a(href   = "https://castanedam.github.io/nicheR/articles/plotting_vignette.html",
                          target = "_blank",
                          icon("brush"), " visualization vignette"),
                   ".")
          )
        ),

        br(),

        fluidRow(
          column(width = 12,
                 tags$a(href   = "https://castanedaM.github.io/nicheR/authors.html",
                        target = "_blank",
                        class = "text-about",
                        icon("user-group"), "If you use nicheR in your research, please cite:"),
                 br(),
                 tags$code(
                   "Castaneda-Guzman M., Hughes C., Paansri P., Cobos M. E. (2026). ",
                   "nicheR: Ellipsoid-based ecological niche modeling. ",
                   "R package version 0.1.0. DOI: ",
                   tags$a("10.32614/CRAN.package.nicheR",
                          href = "https://cran.r-project.org/web/packages/nicheR/index.html",
                          target = "_blank",
                          rel = "noopener noreferrer",
                          style = "color: inherit; text-decoration: underline;"),
                   style = "color: grey;"
                 ),
                 br(), br()

          )
        ),

        fluidRow(
          column(width = 5),
          column(width = 2,
                 actionButton(inputId = "about_start_session_btn",
                              icon = icon("play"),
                              label = "START",
                              class = "btn-warning", width = "150px")
          ),
          column(width = 5)
        ),

        fluidRow(
          column(width = 12,
                 br(), br(),
                 tags$div(style = "font-size: 12px; color: #aaa; padding: 10px 0;
                          border-top: 0.5px solid #ddd; display: flex; gap: 12px; flex-wrap: wrap;",
                          tags$span("nicheR v0.1.0"),
                          tags$span("·"),
                          tags$a(href   = "https://github.com/castanedaM/nicheR/blob/main/LICENSE",
                                 target = "_blank", style = "color: #aaa;", "MIT license"),
                          tags$span("·"),
                          tags$a(href   = "https://github.com/castanedaM/nicheR/issues",
                                 target = "_blank", style = "color: #aaa;", "Report an issue")
                 )
          )
        )
      ),

      tabItem(
        tabName = "build_tab",
        fluidRow(class = "row-tight",
          column(width = 5,
                 tabBox(
                   id = "build_tabs",
                   width = 12,

                   tabPanel(
                     tags$span("Input", class = "text-tab-title"),
                     value = "build_data_tab",
                     fluidRow(
                       column(width = 12,
                              p(instructions$data_input_type, class = "text-instruction"),
                              radioButtons("build_data_input_type_choice",
                                           label = tags$span("Select input type:",
                                                             class = "text-widget-title"),
                                           choiceNames = list(
                                             tags$span("Background layers",
                                                       class = "text-widget-inner"),
                                             tags$span("Previous session",
                                                       class = "text-widget-inner"),
                                             tags$span("Virtual mode",
                                                       class = "text-widget-inner"),
                                             tags$span("Example data",
                                                       class = "text-widget-inner")
                                           ),
                                           choiceValues = c("bg_layers", "prev_session", "virtual_mode", "example_data"),
                                           selected = character(0),
                                           inline = TRUE),
                              uiOutput("build_data_input_type_ui")
                       )
                     )
                   ),

                   tabPanel(
                     title = tags$span("Build", class = "text-tab-title"),
                     value = "build_range_tab",
                     fluidRow(
                       uiOutput("build_variable_selector_ui")
                     ),
                     fluidRow(
                       uiOutput("build_range_method_choice_ui")
                     ),
                     fluidRow(
                       uiOutput("build_covariance_ui")
                     ),
                     fluidRow(
                       uiOutput("build_centroid_mover_ui")
                     ),
                     fluidRow(
                       uiOutput("build_ellipsoid_library_ui"),
                       uiOutput("build_reference_select_ui")
                     ),
                     fluidRow(
                       column(width = 12,
                              uiOutput("build_next_step_ui")
                       )
                     )
                   )
                 )
          ),

          column(width = 7,
                 tabBox(
                   id = "build_plot_tabs",
                   width = 12,

                   tabPanel(
                     title = tags$span("E-space", class = "text-tab-title"),
                     value = "build_espace_plot_tab",
                     uiOutput("build_espace_plot_top_options_ui"),
                     plotOutput("build_espace_plot"),
                     br(),
                     uiOutput("build_espace_plot_bottom_options_ui")
                   ),

                   tabPanel(
                     title = tags$span("G-space", class = "text-tab-title"),
                     value = "build_gspace_plot_tab",
                     uiOutput("build_gspace_plot_top_options_ui"),
                     plotOutput("build_gspace_plot"),
                     br(),
                     uiOutput("build_gspace_plot_bottom_options_ui")
                   ),

                   tabPanel(
                     title = tags$span("Combined", class = "text-tab-title"),
                     value = "build_combined_plot_tab",
                     uiOutput("build_combined_plot_top_options_ui"),
                     plotOutput("build_combined_plot"),
                     br(),
                     uiOutput("build_combined_plot_bottom_options_ui")
                   )

                 ),

                 uiOutput("build_ellipsoid_info_ui"),
                 uiOutput("build_plot_settings_ui")

          )
        )
      ),

      tabItem(
        tabName = "predict_tab",
        fluidRow(class = "row-tight",
          column(width = 5,
                 box(width = 12,
                     title = tags$span("Prediction",
                                       class = "text-tab-title"),
                     uiOutput("predict_ellipsoid_selector_ui"),

                     fluidRow(
                       column(width = 8,
                              tags$span("Prediction layers to include",
                                        class = "text-widget-title"))
                     ),

                     fluidRow(
                       column(width = 6,
                              checkboxInput("predict_suitability",
                                            label = tags$span("Suitability",
                                                              class = "text-widget-inner"),
                                            value = TRUE)),
                       column(width = 6,
                              checkboxInput("predict_suitability_trunc",
                                            label = tags$span("Suitability (truncated)",
                                                              class = "text-widget-inner"),
                                            value = FALSE))
                     ),

                     fluidRow(
                       column(width = 6,
                              checkboxInput("predict_mahalanobis",
                                            label = tags$span("Mahalanobis",
                                                              class = "text-widget-inner"),
                                            value = TRUE)),
                       column(width = 6,
                              checkboxInput("predict_mahalanobis_trunc",
                                            label = tags$span("Mahalanobis (truncated)", class = "text-widget-inner"),
                                            value = FALSE))

                     ),

                     br(),

                     fluidRow(
                       box(title = tags$span("Advanced prediction settings",
                                             class = "text-section-header"),
                           width = 12,
                           collapsible = TRUE,
                           collapsed = TRUE,
                           fluidRow(
                             column(width = 8,
                                    tagList(tags$span("Truncation level adjustment",
                                                      class = "text-widget-title"),
                                            tags$span(icon("circle-info"),
                                                      title = instructions$adjust_trunc_tooltip,
                                                      class = "tooltip-icon"))),
                             column(width = 4,
                                    numericInput(inputId = "predict_adjust_trunc",
                                                 label = NULL,
                                                 value = 0.95,
                                                 min = 0.0001,
                                                 max = 0.99999,
                                                 step = 0.05)
                             )
                           )
                       )
                     ),

                     fluidRow(
                       column(width = 12,
                              div(class = "action-btn-row",
                                  actionButton(inputId = "predict_run_btn",
                                               label = "Predict",
                                               class = "btn-continue")
                              )
                       )
                     ),

                     br(),  br(),

                     fluidRow(
                       uiOutput("predict_ellipsoid_library_ui")
                     ),

                     fluidRow(
                       column(width = 12,
                              uiOutput("predict_next_step_ui")
                       )
                     )
                 )
          ),

          column(width = 7,

                 tabBox(
                   id    = "predict_plot_tabs",
                   width = 12,

                   tabPanel(
                     title = tags$span("E-space", class = "text-tab-title"),
                     value = "predict_espace_plot_tab",
                     uiOutput("predict_espace_plot_top_options_ui"),
                     plotOutput("predict_espace_plot"),
                     br(),
                     uiOutput("predict_espace_plot_bottom_options_ui")
                   ),

                   tabPanel(
                     title = tags$span("G-space", class = "text-tab-title"),
                     value = "predict_gspace_plot_tab",
                     uiOutput("predict_gspace_plot_top_options_ui"),
                     plotOutput("predict_gspace_plot"),
                     br(),
                     uiOutput("predict_gspace_plot_bottom_options_ui")
                   ),

                   tabPanel(
                     title = tags$span("Combined", class = "text-tab-title"),
                     value = "predict_combined_plot_tab",
                     uiOutput("predict_combined_plot_top_options_ui"),
                     plotOutput("predict_combined_plot"),
                     br(),
                     uiOutput("predict_combined_plot_bottom_options_ui")
                   )
                 ),

                 br(),

                 uiOutput("predict_plot_settings_ui")
          )
        )
      ),

      tabItem(
        tabName = "bias_tab",
        fluidRow(class = "row-tight",
          column(width = 5,
                 box(title = tagList(tags$span("Bias", class = "text-tab-title")),
                     width = 12,
                     p(instructions$bias, class = "text-instruction"),
                     fluidRow(column(width = 12,
                                     uiOutput("bias_skip_ui"))),
                     br(), br(),
                     fluidRow(uiOutput("bias_upload_ui")),
                     fluidRow(uiOutput("bias_prepare_ui")),
                     fluidRow(uiOutput("bias_apply_ui")),
                     fluidRow(uiOutput("bias_ellipsoid_library_ui")),
                     fluidRow(
                       column(width = 12,
                              uiOutput("bias_next_step_ui")
                       )
                     )
                 )
          ),

          column(width = 7,

                 tabBox(
                   id = "bias_plot_tabs",
                   width = 12,

                   tabPanel(
                     title = tags$span("Bias Layers",
                                       class = "text-tab-title"),
                     value = "bias_input_layers_plot_tab",
                     plotOutput("bias_layers_plot")
                   ),

                   tabPanel(
                     title = tags$span("Bias Composite",
                                       class = "text-tab-title"),
                     value = "bias_composite_plot_tab",
                     plotOutput("bias_composite_plot")
                   ),

                   tabPanel(
                     title = tags$span("Pred and Biased G-space",
                                       class = "text-tab-title"),
                     value = "bias_gspace_plot_tab",
                     uiOutput("bias_gspace_plot_layer_select_ui"),
                     plotOutput("bias_gspace_plot"),
                     uiOutput("bias_gspace_plot_bottom_options_ui")

                   )

                 ),

                 br(),

                 uiOutput("bias_plot_settings_ui")
          )
        )
      ),

      tabItem(
        tabName = "generate_tab",
        fluidRow(class = "row-tight",
          column(width = 5,
                 fluidRow(
                   uiOutput("generate_controls_ui")
                 ),

                 fluidRow(
                   uiOutput("generate_ellipsoid_library_ui"),
                   uiOutput("generate_occurrence_summary_ui")
                  )
          ),

          column(width = 7,
                 tabBox(
                   id = "generate_plot_tabs",
                   width = 12,

                   tabPanel(
                     title = tags$span("E-space", class = "text-tab-title"),
                     value = "generate_espace_plot_tab",
                     uiOutput("generate_espace_plot_top_options_ui"),
                     plotOutput("generate_espace_plot"),
                     br(),
                     uiOutput("generate_espace_plot_bottom_options_ui")
                   ),

                   tabPanel(
                     title = tags$span("G-space", class = "text-tab-title"),
                     value = "generate_gspace_plot_tab",
                     uiOutput("generate_gspace_plot_top_options_ui"),
                     plotOutput("generate_gspace_plot"),
                     uiOutput("generate_gspace_plot_bottom_options_ui")
                   ),

                   tabPanel(
                     title = tags$span("Combined", class = "text-tab-title"),
                     value = "generate_combined_plot_tab",
                     uiOutput("generate_combined_plot_top_options_ui"),
                     plotOutput("generate_combined_plot"),
                     br(),
                     uiOutput("generate_combined_plot_bottom_options_ui")
                   )

                 ),

                 uiOutput("generate_ellipsoid_info_ui"),
                 uiOutput("generate_plot_settings_ui")

          )
        )
      )
    )
  )
)
