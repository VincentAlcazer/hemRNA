#' fusion UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_fusion_ui <- function(id){
  ns <- NS(id)
  tagList(

    fluidPage(
      h1("Fusion detection"),
      p("Methods: arriba with all filters. High confidence are shown by default."),

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

#' fusion Server Functions
#'
#' @noRd
mod_fusion_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    df_fusion <- reactive({r$test$df_fusion})

    plot <- reactive({

      req(df_fusion())

      plot <- df_fusion() %>%
        filter(confidence %in% c(input$confidence)) %>%
        ggplot(aes(x = sample_id, y = fusion, fill = confidence)) +
        geom_tile() +
        labs(x="", y="")

      return(plot)

    })

    height <- reactive({
      if(length(unique(df_fusion()[,"fusion"]) < 20)){
        height = 800
      } else if(length(unique(df_fusion()[,"fusion"]) < 40)){
        height = 1200
      } else if(length(unique(df_fusion()[,"fusion"]) < 60)){
        height = 1600
      } else {
        height = 2000
      }


    })


    output$tileplot <- renderPlot(height = height,{
      plot() +
        default_theme +
        theme(
          legend.position = input$legend_ext,
          axis.text.y = element_text(size = input$y_size, face = "bold"),
          axis.text.x = element_text(size = input$x_size, angle=90, vjust = 0.5, hjust = 1, face = "bold"))
    })

    output$result_table <- DT::renderDT(

      df_fusion() %>%
        select(sample_id, gene1, gene2, confidence, type, reading_frame, everything()) %>%
        filter(confidence %in% c(input$confidence)), # data
      class = "display nowrap compact", # style
      filter = "top", # location of column filters
      server = T,
      rownames = FALSE,
      options = list(
        scrollX = TRUE,
        "pagelength" = 20,
        lengthChange = TRUE,
        columnDefs = list(list(className = "dt-left", targets = "_all"))
      )
    )


    output$download_table <- downloadHandler(
      filename = function() {
        paste("Arriba_fusion.tsv")
      },
      content = function(file) {
        write.table(df_fusion() %>%
                            select(sample_id, gene1, gene2, confidence, type, reading_frame, everything()) %>%
                            filter(confidence %in% c(input$confidence)), file, row.names = FALSE, sep = "\t", quote = F)
      }
    )


  })
}

## To be copied in the UI
# mod_fusion_ui("fusion_ui_1")

## To be copied in the server
# mod_fusion_server("fusion_ui_1")
