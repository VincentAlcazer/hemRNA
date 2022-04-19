#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
#' @import shinyFiles
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic

    dashboardPage(

      header = dashboardHeader(
        title = "hemRNA v0.82"
      ),

      #################### ==================== SIDEBAR ====================  ####################

      sidebar = dashboardSidebar(
        shinydashboard::sidebarMenu(
          # Setting id makes input$tabs give the tabName of currently-selected tab
          id = "tabs",
          menuItem("Home", tabName = "home", icon = icon("play-circle")),
          menuItem("Data loading", tabName = "data", icon = icon("spinner")),
          menuItem("Results overview", tabName = "overview", icon = icon("circle-notch")),
          menuItem("Expression", tabName = "expression", icon = icon("poll")),
          menuItem("RNAmut", tabName = "RNAmut", icon = icon("bezier-curve")),

          menuItem("Fusion", tabName = "fusion", icon = icon("exchange-alt"),
                   menuSubItem("Arriba", tabName = "arriba"),
                   menuSubItem("NF-Core", tabName = "nfcore")
                   ),
          menuItem("HotSpot variants", tabName = "hotspot",icon = icon("fire")),
          menuItem("CNV", tabName = "CNV", icon = icon("copy"))

          )
      ), # sidebar


      #################### ==================== BODY ====================  ####################

      body = dashboardBody(
        tabItems(

          ##### ===== Home

          tabItem("home", mod_home_ui("home_ui_1")),

          ##### ===== Data loading

          tabItem("data",
                             fluidPage(
                               h1("Data loading"),
                               p("hemRNA has been optimized to run with the output of its corresponding",
                                 a("RNA-seq bash pipeline.", href = "https://github.com/VincentAlcazer/hemRNA")),

                               box(title = "Folder", width = 6,
                                   ### === Load results folder
                                   HTML("<b>Select results folder</b></br>"),
                                   shinyDirButton('folder_path', label='Browse...', title='Please select a folder'),
                                   div(style = "margin-top: 20px"),

                                   ### === Load bed
                                   selectInput("panel", label = "Select a bed panel",
                                               choices = c("None"),
                                               selected = "None"),
                                   div(style = "margin-top: -20px"),
                                   p("Panels should be added to the extda/bed_panels folder as bed files to
                                     be loaded directly (the first file in alphabetical order will be
                                     loaded as default)."),

                                   div(style = "margin-top: 20px"),

                                   h3("File check:"),

                                   htmlOutput("salmon_check"),
                                   shinycssloaders::withSpinner(htmlOutput("salmon_dim"),type = 1),
                                   shinycssloaders::withSpinner(htmlOutput("panel_dim"),type = 1),

                                   htmlOutput("arriba_check"),
                                   htmlOutput("fusion_catcher_check"),

                                   htmlOutput("hotspot_check"),
                                   htmlOutput("RNAmut_check"),

                                   htmlOutput("cnv_check")

                                    # p("Preview"),
                                    # downloadButton("download_table", "Download table (.tsv)"),
                                    # shinycssloaders::withSpinner(DT::DTOutput("preview_data"),type = 6)

                               ),# box
                               box(title = "Meta", width = 6,
                                   p("/!\ First column should correspond to samples / patients ID"),
                                   fileInput("meta",
                                             label = ("Meta"),
                                             accept = c(
                                               "text/tab-separated-values",
                                               "text/comma-separated-values",
                                               "text/plain",
                                               "text/csv",
                                               ".csv",
                                               ".tsv",
                                               ".xls",
                                               ".xlsx"
                                             )
                                   ),
                                   div(style = "margin-top: -20px"),
                                   radioButtons("sep", "Separator",
                                                choices = c(
                                                  "comma-delim (.csv1)" = ",",
                                                  "semi colon-delim (.csv2)" = ";",
                                                  "tab-delim (.tsv/.txt)" = "\t",
                                                  "Excel (.xls/.xlsx)" = "xl"
                                                ),
                                                selected = "\t"
                                   ),
                                   div(style = "margin-top: -20px"),
                                   radioButtons("dec", "Decimal",
                                                choices = c(
                                                  "Comma (,)" = ",",
                                                  "Period (.)" = "."
                                                ),
                                                selected = ","
                                   ),
                                   p("Preview"),
                                   shinycssloaders::withSpinner(DT::DTOutput("preview_meta"),type = 6)
                               ), #box
                            box(title = "Missing genes", width = 6,
                                p("Genes/Features present in panel but not found in exp. dataset:"),
                                htmlOutput("missing_genes")
                                )



                             ) # fluid page
                  ), #tabItem

          ##### ===== Overview

          tabItem("overview",
                  mod_overview_ui("overview_1")
          ),#tabItem

          ##### ===== Expression

          tabItem("expression",
                  h1("Expression"),
                  HTML("<b>Methods:</b> TPM (transcripts per million) are quantified using Salmon and log-normalized (Log2 TPM+1).<br>
                (NB: TPM is a normalization accouting for both sequencing depth and transcripts size.
                When you use TPM, the sum of all TPMs in each sample are the same.
                This makes it easier to compare the proportion of reads that mapped to a gene in each sample.)"),
                  tabsetPanel(
                    id = "expression", type = "tabs",

                    tabPanel("PCA",  mod_expression_ui("expression_ui_1")),
                    tabPanel("Individual xp",mod_expression_individual_ui("expression_individual_ui_1")),
                    tabPanel("Signatures", mod_expression_signatures_ui("expression_signatures_ui_1")),
                    #tabPanel("Heatmap",  mod_expression_heatmap_ui("expression_heatmap_ui_1")),
                    tabPanel("DESEQ2", mod_expression_deseq_ui("expression_deseq_ui_1"))

                  )#TabSetPanel

        ), #tabItem
        tabItem("RNAmut", mod_RNAmut_ui("RNAmut_ui_1")),

        tabItem("arriba", mod_fusion_ui("fusion_ui_1")),
        tabItem("nfcore",
                h1("NF-Core"),
                p("Results from NF-Core"),
                tabsetPanel(
                  id = "nfcore_tabs", type = "tabs",
                  tabPanel("Fusion_catcher", mod_fusion_catcher_ui("fusion_catcher_ui_1")),
                  tabPanel("star_fusion", mod_star_fusion_ui("star_fusion_ui_1"))
                )
                ),#tabitem

        tabItem("hotspot",mod_hotspot_ui("hotspot_ui_1")),

        tabItem("CNV",
                mod_CNV_ui("CNV_ui_1"))

        )#tabItems
      ) #body

    ) #dashboard page



  ) #taglist
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){

  add_resource_path(
    'www', app_sys('app/www')
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'hemRNA'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}

