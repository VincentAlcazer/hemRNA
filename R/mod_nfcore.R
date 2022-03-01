#' nfcore UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_nfcore_ui <- function(id){
  ns <- NS(id)
  tagList(

    fluidPage(
      h1("NF-Core"),
      p("Methods: NF-Core."),

      tabsetPanel(
        id = "fusion", type = "tabs",
        tabPanel("Overall graph",

                 column(8, shinycssloaders::withSpinner(plotOutput(ns("tileplot")),type=6)
                 ),
                 column(1)


        ),#tabsetpanel
        tabPanel("Result table",
                 column(8,
                        downloadButton(ns("download_table"), "Download table (.tsv)"),
                        shinycssloaders::withSpinner(DT::DTOutput(ns("result_table")),type=6)
                 ),
                 column(1)
        )
      ),##tabsetpanel

      column(2,
             absolutePanel(
               width = 200, right = 20, draggable = T,
               style = "opacity: 0.85",
               wellPanel(
                 checkboxGroupInput(ns("confidence"),
                                    label = c("Confidence"),
                                    choices = c("High"="high",
                                                "Medium"="medium",
                                                "Low"="low"),
                                    selected = "high"),
                 sliderInput(ns("y_size"), label = "y-axis font size",
                             min = 1, max = 30, value = 12, step= 1),
                 sliderInput(ns("x_size"), label = "x-axis font size",
                             min = 1, max = 30, value = 14, step= 1),

                 selectInput(ns("legend_ext"),
                             label = ("External legend"),
                             choices = c(
                               "No" = "none",
                               "Top" = "top",
                               "Right" = "right",
                               "Left" = "left",
                               "Bottom" = "bottom"
                             ),
                             multiple = F, selected = "right"
                 )
               )
             ) # Absolutepanel
      ) # Column

    ) # Fluidpage

  )
}

#' nfcore Server Functions
#'
#' @noRd
mod_nfcore_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

  })
}

## To be copied in the UI
# mod_nfcore_ui("nfcore_ui_1")

## To be copied in the server
# mod_nfcore_server("nfcore_ui_1")
