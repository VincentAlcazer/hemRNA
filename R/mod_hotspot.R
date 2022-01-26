#' hotspot UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_hotspot_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      h1("Hotspot"),
      p("Methods: For each  unique mutations, the overall raw counts
        or the normalized count (total raw count / total Depth) is represented."),

      tabsetPanel(
        id = "hotspot", type = "tabs",
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
                 radioButtons(ns("count_type"),
                              label = c("Count"),
                              choices = c("Raw counts"="overall_count",
                                          "Norm. counts (VAF)"="overall_percent"),
                              selected = "overall_percent"),
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

#' hotspot Server Functions
#'
#' @noRd
mod_hotspot_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    df_hotspot <- reactive({r$test$df_hotspot})

    # df_stat <- reactive({
    #
    #   req(df_hotspot())
    #
    #   df_stat <- df_hotspot() %>%
    #     group_by(sample_id, gene_name) %>%
    #     summarise(tot_raw_count = sum(overall_count, na.rm = T),
    #               tot_percent = sum(overall_percent, na.rm = T))
    #
    #   return(df_stat)
    #
    # })

    plot <- reactive({

      req(df_hotspot())

      plot <- df_hotspot() %>%
        filter(overall_percent>input$vaf[1] & overall_percent<=input$vaf[2]) %>%
        ggplot(aes_string(x="sample_id", y = "gene_mut", fill = input$count_type)) +
        geom_tile() +
        labs(x="", y="", fill = "Cum. VAF")

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


    output$result_table <- DT::renderDT(

      df_hotspot() %>%
        filter(overall_percent>input$vaf[1] & overall_percent<=input$vaf[2]) %>%
         arrange(desc(overall_percent)), # data
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
        paste("Hotspot.tsv")
      },
      content = function(file) {
        write.table(df_hotspot() %>% arrange(desc(norm_count)), file, row.names = FALSE, sep = "\t", quote = F)
      }
    )


  })
}

## To be copied in the UI
# mod_hotspot_ui("hotspot_ui_1")

## To be copied in the server
# mod_hotspot_server("hotspot_ui_1")
