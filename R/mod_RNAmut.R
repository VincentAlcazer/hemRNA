#' RNAmut UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_RNAmut_ui <- function(id){
  ns <- NS(id)
  tagList(

    fluidPage(
      h1("RNAmut"),
      p("Methods: RNAmut pipeline (Gu M. et al. Haematologica 2020).
        /!\\ Fusions VAF are provided as estimations only! (n mutated reads / total gene 1+2 WT reads)"),

      tabsetPanel(
        id = "RNAmut", type = "tabs",
        tabPanel("Graph",

                 column(8, shinycssloaders::withSpinner(plotOutput(ns("tileplot"), height = 800),type=6)
                 ),
                 column(1)


        ),#tabsetpanel
        tabPanel("Results table",
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
                 checkboxGroupInput(ns("impact"),
                               label = c("Mutation impact"),
                               choices = c("Oncogenic",
                                           "Neutral",
                                           "Sequencing_artefact"
                                           ),
                               selected = "Oncogenic"
                               ),
                 sliderInput(ns("vaf"), label = "VAF ranges",
                             min = 0, max = 1, value = c(0,1), step= 0.01),
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

#' RNAmut Server Functions
#'
#' @noRd
mod_RNAmut_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    df_RNAmut <- reactive({r$test$df_RNAmut})

    plot <- reactive({

      req(df_RNAmut())

      plot <- df_RNAmut() %>%
        filter(Info %in% input$impact) %>%
        filter(VAF2 >= input$vaf[1] & VAF2 <=input$vaf[2]) %>%
        ggplot(aes(x=sample_id, y = forcats::fct_rev(forcats::fct_reorder(gene_mut,as.numeric(type))), fill = VAF2)) +
        geom_tile(color = "white", size = 0.25) +
        labs(x="", y="", fill = "VAF")   +
        facet_grid(type~.,  space="free", scales = "free_y", switch = "y")

      return(plot)

    })

    output$tileplot <- renderPlot({
      plot() +
        default_theme +
        scale_fill_viridis_c() +
        theme( legend.position = input$legend_ext,
               strip.text.y = element_text(size = 14, face = "bold"),
               axis.text.y = element_text(size = input$y_size),
               axis.text.x = element_text(size = input$x_size, angle=90, vjust = 0.5, hjust = 1))
    })


    output$result_table <- DT::renderDT(

      df_RNAmut() %>%
        filter(Info %in% input$impact) %>%
        filter(VAF2 >= input$vaf[1] & VAF2 <=input$vaf[2]) %>%
        arrange(desc(VAF2)), # data
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
        paste("RNA_mut.tsv")
      },
      content = function(file) {
        write.table(df_RNAmut() %>%
                      filter(Info %in% input$impact) %>%
                      filter(VAF2 >= input$vaf[1] & VAF2 <=input$vaf[2]), file, row.names = FALSE, sep = "\t", quote = F)
      }
    )


  })
}

## To be copied in the UI
# mod_RNAmut_ui("RNAmut_ui_1")

## To be copied in the server
# mod_RNAmut_server("RNAmut_ui_1")
