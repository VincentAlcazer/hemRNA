#' expression_signatures UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_expression_signatures_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      h1("Signatures"),
      column(8,
             textOutput(ns("missing_genes")),
             shinycssloaders::withSpinner(plotOutput(ns("barplot")),type=6)
             #shinycssloaders::withSpinner(DT::DTOutput(ns("result_table")),type=6)
      ),
      column(1),
      column(2,
             absolutePanel(
               width = 200, right = 20, draggable = T,
               style = "opacity: 0.85",
               wellPanel(
                 selectInput(ns("Signatures"),
                             label = ("Signatures"),
                             multiple = F, selected = NULL,
                             ""
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
                             multiple = F, selected = "right"
                 )
               )
             ) # Absolutepanel
      ) # Column

    ) # FluidRow

  )
}

#' expression_signatures Server Functions
#'
#' @noRd
mod_expression_signatures_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    df_filt <- reactive({r$test$df_filt})

    sig_list <- reactive({r$test$sig_list})

    observe({
      updateSelectInput(
        session,
        "Signatures",
        choices = c(names(sig_list()))
      )
    })


    sig_df <- reactive({

    req(df_filt())
    input$Signatures
      isolate({

        genes_sig <- as.character(unlist(dplyr::select(sig_list()[[input$Signatures]], 1)))
        genes_xp <- as.character(unlist(dplyr::select(df_filt(), 1)))

        intersect_genes <- intersect(genes_sig, genes_xp)

        coefs_df <- sig_list()[[input$Signatures]] %>% column_to_rownames("gene_id") %>% as.matrix()

        coefs <- as.numeric(unlist(coefs_df[intersect_genes,]))

        gene_xp <- df_filt() %>% filter(gene_name %in% intersect_genes) %>%
          group_by(gene_name) %>% summarise_all(mean) %>%
          column_to_rownames("gene_name") %>% as.matrix()

        sig_calc <- gene_xp[intersect_genes,] * coefs

        sig_final <- data.frame(Score = colSums(sig_calc)) %>%
          rownames_to_column("sample_id")
          # mutate(Median_cut = factor(cut(Score,2),labels = c("Low","High")),
          #        Tercile_cut = factor(cut(Score,3),labels = c("Low","Intermediate", "High")),
          #        Quartile_cut = factor(cut(Score,4),labels = c("Low","Intermediate1","Intermediate2", "High")))

        return(sig_final)
      })
    })


    missing_genes <- reactive({

      req(df_filt())
      isolate({

        genes_sig <- as.character(unlist(dplyr::select(sig_list()[[input$Signatures]], 1)))
        genes_xp <- as.character(unlist(dplyr::select(df_filt(), 1)))

        missing_genes <- setdiff(genes_sig, genes_xp)


        return(missing_genes)
      })
    })

    output$missing_genes <- renderText({

      c("Genes found in signature but not in gene expression dataset: ",
             missing_genes())

    })


    ##### ===== Plots


    plot <- reactive({

      req(sig_df())

      plot <- sig_df() %>%
        ggplot(aes(x = forcats::fct_reorder(sample_id,-Score), y = Score, fill = Score)) +
        geom_bar(stat = "identity", color = "black") +
        labs(x="", y="", fill = paste0(input$Signatures," Score")) +
        scale_fill_viridis_c()

      return(plot)

    })

    output$barplot <- renderPlot(height = 600,{
      plot() +
        default_theme +
        theme(
          legend.position = input$legend_ext,
          axis.text.y = element_text(size = input$y_size, face = "bold"),
          axis.text.x = element_text(size = input$x_size, angle=90, vjust = 0.5, hjust = 1, face = "bold")
        )
    })

    # output$result_table <- DT::renderDT(
    #
    #   df_filt(), # data
    #   class = "display nowrap compact", # style
    #   filter = "top", # location of column filters
    #   server = T,
    #   rownames = FALSE,
    #   options = list(
    #     scrollX = TRUE,
    #     "pagelength" = 20,
    #     lengthChange = TRUE,
    #     columnDefs = list(list(className = "dt-left", targets = "_all"))
    #   )
    # )

  })
}

## To be copied in the UI
# mod_expression_signatures_ui("expression_signatures_ui_1")

## To be copied in the server
# mod_expression_signatures_server("expression_signatures_ui_1")
