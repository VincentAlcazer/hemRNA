#' variant UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_variant_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      tabsetPanel(
        id = "variant", type = "tabs",
        tabPanel("Overall graph",

                 column(8, shinycssloaders::withSpinner(plotOutput(ns("tileplot")),type=6)
                 ),
                 column(1)


        ),#tabsetpanel
        tabPanel("Result table",
                 column(8, shinycssloaders::withSpinner(DT::DTOutput(ns("result_table")),type=6)
                 ),
                 column(1)
        )
      ),##tabsetpanel

      column(2,
             absolutePanel(
               width = 200, right = 20, draggable = T,
               style = "opacity: 0.85",
               wellPanel(
                 radioButtons(ns("count_type"),
                                    label = c("Count"),
                                    choices = c("Raw counts"="tot_raw_count",
                                                "Norm counts"="tot_norm_count"),
                                    selected = "tot_raw_count"),
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

#' variant Server Functions
#'
#' @noRd
mod_variant_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    df_hotspot <- reactive({r$test$df_hotspot})

    # df_stat <- reactive({
    #
    #   req(df_hotspot())
    #
    #   df_stat <- df_hotspot() %>%
    #     group_by(sample_id, gene_name) %>%
    #     summarise(tot_raw_count = sum(raw_count, na.rm = T),
    #               tot_norm_count = sum(norm_count, na.rm = T))
    #
    #   return(df_stat)
    #
    # })

    plot <- reactive({

      req(df_stat())

      plot <- df_stat() %>%
        ggplot(aes_string(x="sample_id", y = "gene_name", fill = input$count_type)) +
        geom_tile() +
        labs(x="", y="")

      return(plot)

    })

    output$tileplot <- renderPlot(height = 800,{
      plot() +
        default_theme +
        scale_fill_viridis_c() +
        theme( legend.position = input$legend_ext,
               axis.text.y = element_text(size = input$y_size),
               axis.text.x = element_text(size = input$x_size, angle=90, vjust = 0.5, hjust = 1))
    })


  })
}

## To be copied in the UI
# mod_variant_ui("variant_ui_1")

## To be copied in the server
# mod_variant_server("variant_ui_1")
