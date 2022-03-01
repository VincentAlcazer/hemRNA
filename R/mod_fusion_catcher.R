#' fusion_catcher UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_fusion_catcher_ui <- function(id){
  ns <- NS(id)
  tagList(

    fluidPage(
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
                 checkboxGroupInput(ns("tools"),
                                    label = c("Tools"),
                                    choices = c("Bowtie + STAR"="BOWTIE\\+STAR",
                                                "Bowtie"="BOWTIE",
                                                "STAR"="STAR"),
                                    selected = "BOWTIE\\+STAR"),
                 numericInput(ns("spanning_reads"), label = "Min. unique spanning reads",
                              min = 0, value = 10),
                 numericInput(ns("common_reads"), label = "Min. common mapping reads",
                             min = 0, value = 0),
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

#' fusion_catcher Server Functions
#'
#' @noRd
mod_fusion_catcher_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    df_fusion <- reactive({r$test$df_fusioncatcher})

    plot <- reactive({

      req(df_fusion())

      plot <- df_fusion() %>%
        filter(grepl(input$tools, Fusion_finding_method),
               Spanning_unique_reads >= input$spanning_reads,
               Counts_of_common_mapping_reads >= input$common_reads) %>%
        ggplot(aes(x = sample_id, y = fusion, fill = Fusion_finding_method)) +
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
        select(sample_id, gene1, gene2, Fusion_finding_method, type, reading_frame, everything()) %>%
      filter(grepl(input$tools, Fusion_finding_method),
             Spanning_unique_reads >= input$spanning_reads,
             Counts_of_common_mapping_reads >= input$common_reads), # data
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
        paste("Fusioncatcher.tsv")
      },
      content = function(file) {
        write.table(df_fusion() %>%
                            select(sample_id, gene1, gene2, Fusion_finding_method, type, reading_frame, everything()),
        file, row.names = FALSE, sep = "\t", quote = F)
      }
    )

  })
}

## To be copied in the UI
# mod_fusion_catcher_ui("fusion_catcher_ui_1")

## To be copied in the server
# mod_fusion_catcher_server("fusion_catcher_ui_1")
