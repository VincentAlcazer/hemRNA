#' expression_individual UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_expression_individual_ui <- function(id){
  ns <- NS(id)
  tagList(
    tabsetPanel(
      id = "xp_indiv", type = "tabs",
      tabPanel("Graph (TPM)",
               br(),
              column(8,
                     shinycssloaders::withSpinner(plotOutput(ns("scatterplot")),type=6)

              ),
              column(1)
      ),#tabsetpanel
    tabPanel("Results table (TPM)",
             br(),
             column(8,
                    downloadButton(ns("download_table"), "Download table (selected genes) (.tsv)"),
                    downloadButton(ns("download_table_all"), "Download full table (all genes) (.tsv)"),
                    shinycssloaders::withSpinner(DT::DTOutput(ns("result_table")),type=6)
             ),
             column(1)
    ),
    tabPanel("Graph (ratio)",
             p("The mean log-normalized TPM value of the selected control genes is calculated as a control factor.
                 The ratio of the log-normalized TPM value of the selected genes against the control factor
                 is then represented."),
             column(4,offset = 2,
                    wellPanel(
                      p("Select one or several control genes"),
                      selectizeInput(ns("Genes_ratio"),
                                     label = ("Control genes"),
                                     multiple = T, selected = NULL,
                                     ""
                      ),
                    ) #wellpanel
             ),
             column(8,
                    shinycssloaders::withSpinner(plotOutput(ns("scatterplot_ratio")),type=6)

             ),
             column(1)
    ),#tabsetpanel
    tabPanel("Results table (ratio)",
             column(8,
                    downloadButton(ns("download_table_ratio"), "Download table (.tsv)"),
                    shinycssloaders::withSpinner(DT::DTOutput(ns("result_table_ratio")),type=6)
             ),
             column(1)
    )
  ),##tabsetpanel


    column(2,
           absolutePanel(
             width = 200, right = 20, draggable = T,
             style = "opacity: 0.85",
             wellPanel(
               selectizeInput(ns("Genes"),
                           label = ("Genes"),
                           multiple = T, selected = NULL,
                           ""
               ),
               #actionButton(ns("run"), "Run"),
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

  ) # FluidRow

}

#' expression_individual Server Functions
#'
#' @noRd
mod_expression_individual_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    df_filt <- reactive({r$test$df_filt})

    observe({
      updateSelectizeInput(
        session,
        "Genes",
        choices = c(unique(df_filt()[,1])),
        server = T
      )
    })

    observe({
      updateSelectizeInput(
        session,
        "Genes_ratio",
        choices = c(unique(df_filt()[,1])),
        server = T
      )
    })

    ### TPM graph & table

    plot_df <- reactive({

       req(df_filt())
       #req(input$run >= 1)
       input$Genes

      isolate({

        dat <- df_filt() %>%
          filter(gene_name %in% input$Genes) %>%
          pivot_longer(-gene_name, names_to = "sample_id",values_to = "Log2TPM")

        return(dat)

      })
  })

      plot <- reactive({

        req(plot_df())
        #req(input$run >= 1)
        plot <- plot_df() %>%
          ggplot(aes(x = sample_id, y = Log2TPM, color = gene_name)) +
          geom_boxplot(alpha=0.25, size = 0.75, color = "grey", outlier.shape=NA) +
          geom_jitter(width = 0.2, size = 3) +
          labs(x="", y="Log2 (TPM+1)")

        return(plot)

      })

      output$scatterplot <- renderPlot(height = 600,{
        plot() +
          default_theme +
          theme(
            legend.position = input$legend_ext,
            axis.text.y = element_text(size = input$y_size, face = "bold"),
            axis.text.x = element_text(size = input$x_size, angle=90, vjust = 0.5, hjust = 1, face = "bold")
          )
      })

      output$result_table <- DT::renderDT(
        plot_df(), # data
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
          paste("TPM_long.tsv")
        },
        content = function(file) {
          write.table(plot_df(), file, row.names = FALSE, sep = "\t", quote = F)
        }
      )

      output$download_table_all <- downloadHandler(
        filename = function() {
          paste("TPM_all.tsv")
        },
        content = function(file) {
          write.table(df_filt(), file, row.names = FALSE, sep = "\t", quote = F)
        }
      )


      ##### ===== Ratio graph & tables

      plot_df_ratio <- reactive({

        req(df_filt())
        #req(input$run >= 1)

          dat <- df_filt() %>% distinct(gene_name, .keep_all = T) %>%
            filter(gene_name %in% c(input$Genes,input$Genes_ratio)) %>%
            column_to_rownames("gene_name") %>% t() %>% data.frame(check.names = F) %>%
            rownames_to_column("sample_id") %>%
            rowwise() %>%
            mutate(control = mean(c(!!!rlang::syms(input$Genes_ratio)), na.rm=T)) %>%
            pivot_longer(input$Genes, names_to = "gene_name",values_to = "Log2TPM") %>%
            mutate(ratio = Log2TPM / control) %>%
            ungroup()

          return(dat)

        })

      plot_ratio <- reactive({

        req(plot_df_ratio())
        #req(input$run >= 1)
        plot <- plot_df_ratio() %>%
          ggplot(aes(x = sample_id, y = ratio, color = gene_name)) +
          geom_boxplot(alpha=0.25, size = 0.75, color = "grey", outlier.shape=NA) +
          geom_jitter(width = 0.2, size = 3) +
          labs(x="", y="Log2 (TPM+1) ratio")

        return(plot)

      })

      output$scatterplot_ratio <- renderPlot(height = 600,{
        #req(plot_ratio())
        plot_ratio() +
          default_theme +
          theme(
            legend.position = input$legend_ext,
            axis.text.y = element_text(size = input$y_size, face = "bold"),
            axis.text.x = element_text(size = input$x_size, angle=90, vjust = 0.5, hjust = 1, face = "bold")
          )
      })

      output$result_table_ratio <- DT::renderDT(
        plot_df_ratio(), # data
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


      output$download_table_ratio <- downloadHandler(
        filename = function() {
          paste("TPM_long.tsv")
        },
        content = function(file) {
          write.table(plot_df(), file, row.names = FALSE, sep = "\t", quote = F)
        }
      )

  })

}

## To be copied in the UI
# mod_expression_individual_ui("expression_individual_ui_1")

## To be copied in the server
# mod_expression_individual_server("expression_individual_ui_1")
