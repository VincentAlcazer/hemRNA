#' CNV UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_CNV_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      h1("CNV"),
      HTML("Methods: CNVkit"),
      tabsetPanel(
        id = "cnv", type = "tabs",

        tabPanel("Graph",
                 br(),

                 column(10,
                        shinycssloaders::withSpinner(plotOutput(ns("lineplot"), height = 600),type=6)
                 )
                 # column(10,
                 # ),


        ),#tabsetpanel

        tabPanel("Result table",
                 downloadButton(ns("download_table"), "Download table (.tsv)"),
                 column(8, shinycssloaders::withSpinner(DT::DTOutput(ns("preview_data")),type=6)
                 ),
                 column(1)
        )
      ),##tabsetpanel

      column(2,
             absolutePanel(
               width = 200, right = 20, draggable = T,
               style = "opacity: 0.85",
               wellPanel(
                   selectizeInput(ns("sample"),
                                label = ("Sample"),
                                multiple = F, selected = NULL,
                                ""
                 ),
                 selectizeInput(ns("method"),
                                label = ("Method"),
                                multiple = F, selected = "cbs",
                                choices=c("base","cbs", "hmm", "cnr_by_cytob")
                 ),
                 selectizeInput(ns("values"),
                                label = ("Values"),
                                multiple = F, selected = "raw",
                                choices=c("raw","mean by cytoband")
                 ),

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
                             multiple = F, selected = "top"
                 )
               )
             ) # Absolutepanel
      ) # Column

    ) # Fluidpage
  )
}

#' CNV Server Functions
#'
#' @noRd
mod_CNV_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    cnvkit_list <- reactive({r$test$cnvkit_list})

    observe({
      updateSelectizeInput(
        session,
        "sample",
        choices = c(unique(names(cnvkit_list()))),
        selected = c(unique(names(cnvkit_list())))[1],
        server = T
      )
    })

  ### === Plots

    plot <- reactive({

      req(cnvkit_list())
      req(hg19_cytoband)

      ### Points df
      datapoints <- cnvkit_list()[[input$sample]]$cnr_by_cytob

      x_breaks <- datapoints %>%
        group_by(chr_cytob3) %>%
        summarise(min_order = min(order),
                  mean_order = mean(order),
                  max_order = max(order)) %>%
        na.omit()

    ### Segmentation df

      if (input$method == "cnr_by_cytob"){

        datalines <- cnvkit_list()[[input$sample]][[input$method]] %>%
          group_by(chr_cytob3) %>%
          mutate(mean_cytob_log2 = mean(log2)) %>%
          ungroup()

        datalines$order_max <- c(datalines$order[-1],max(datapoints$order))

      } else {

        datalines <- cnvkit_list()[[input$sample]][[input$method]] %>%
          select(-gene) %>%  mutate(n_line = seq(1:nrow(.))) %>%
          as.data.frame() %>%
          left_join(select(datapoints,chromosome, start, order, chr_cytob3 ), by = "chromosome",
                    suffix = c("_line","_point")) %>%
          mutate(dist = start_line - start_point) %>% arrange(abs(dist)) %>%
          distinct(n_line, .keep_all = T) %>% arrange(n_line) %>%
          group_by(chr_cytob3) %>%
          mutate(mean_cytob_log2 = mean(log2)) %>%
          ungroup()

        datalines$order_max <- c(datalines$order[-1],max(datapoints$order))

      }



  ### Plot

      ylim = if_else(abs(max(datalines$log2)) >= 1.5, round(max(datalines$log2),1),1.5)
      message("plot")
      plot <- ggplot() +
        geom_point(data = datapoints,
                   aes(x=order, y=log2), alpha=0.2) +
        geom_linerange(data = datalines,
                       aes_string(xmin = "order", xmax = "order_max",
                                  y = if_else(input$values == "raw", "log2","mean_cytob_log2")),
                       size = 3, color = "orange") +
        ylim(-ylim, ylim) +
        geom_vline(data=datapoints,
                   aes(xintercept = vline), size = 1, linetype = 2, alpha = 0.75) +
        labs(x = "chr", y = "Copy ratio (log2)",
             title = paste0(input$sample," - ",input$method," segmentation (",input$values," values)")) +
        scale_x_continuous(breaks =x_breaks$mean_order,  labels = x_breaks$chr_cytob3)

      return(plot)

    })

    output$lineplot <- renderPlot({
      plot() +
        default_theme +
        theme(
          legend.position = input$legend_ext,
          axis.text.y = element_text(size = input$y_size, face = "bold"),
          axis.text.x = element_text(size = input$x_size,
                                     angle = 90, vjust = 0.5, hjust = 1,
                                     face = "bold")
        )
    })


    output$preview_data <- DT::renderDT(
      cnvkit_list()[[input$sample]]$cnr, # data
      class = "display nowrap compact", # style
      filter = "top", # location of column filters
      server = T,
      rownames = FALSE,
      options = list(
        scrollX = TRUE,
        lengthChange = TRUE,
        columnDefs = list(list(className = "dt-left", targets = "_all"))
      )
    )

    output$download_table <- downloadHandler(
      filename = function() {
        paste("CNV.tsv")
      },
      content = function(file) {
        write.table(cnvkit_list()[[input$sample]]$cnr, file, row.names = FALSE, sep = "\t", quote = F)
      }
    )


  })
}

## To be copied in the UI
# mod_CNV_ui("CNV_ui_1")

## To be copied in the server
# mod_CNV_server("CNV_ui_1")
