#' overview UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_overview_ui <- function(id){
  ns <- NS(id)
  tagList(

    fluidPage(
      h1("Results overview"),
      p("Methods:"),

      column(10,
             shinycssloaders::withSpinner(BioCircos::BioCircosOutput(ns("circosplot"), height = 900),type = 6)
      ),

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
                                label = ("CNV Method"),
                                multiple = F, selected = "cnr_by_cytob",
                                choices=c("cnr_by_cytob")
                 ),
                 selectizeInput(ns("values"),
                                label = ("CNV Values"),
                                multiple = F, selected = "mean by cytoband",
                                choices=c("raw","mean by cytoband")
                 )

               )
             ) # Absolutepanel
      ) # Column

    ) # Fluidpage

  )
}

#' overview Server Functions
#'
#' @noRd
mod_overview_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    ### DF

    cnvkit_list <- reactive({r$test$cnvkit_list})

    bed <- reactive({r$test$bed})

    df_xp <- reactive({r$test$df_xp})

    df_RNAmut <- reactive({r$test$df_RNAmut})

    df_arriba <- reactive({r$test$df_arriba})


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

    df_cnv <- reactive({

      req(cnvkit_list())
      req(hg19_cytoband)

      ### Points df
      datapoints <- cnvkit_list()[[input$sample]]$cnr_by_cytob

      ### Segmentation df

      datalines <- cnvkit_list()[[input$sample]][[input$method]] %>%
          group_by(chr_cytob3) %>%
          mutate(mean_cytob_log2 = mean(log2)) %>%
          ungroup()

        datalines$order_max <- c(datalines$order[-1],max(datapoints$order))



      dat <- datalines %>%
        rowwise() %>%
        mutate(CNV = if_else(input$values == "raw", log2, mean_cytob_log2),
               gene_name = gene)

      return(dat)

    })


    output$circosplot <- BioCircos::renderBioCircos({
      circos_plot(expression = df_xp() %>% select(gene_name, expression = all_of(input$sample)),
                  cnv = df_cnv(),
                  variant = df_RNAmut() %>% filter(sample_id == input$sample),
                  arriba = df_arriba() %>% filter(sample_id == input$sample),
                  bed(),
                  exp_max_rad = 0.97,exp_min_rad = 0.80,
                  CNV_max_rad = 0.75, CNV_min_rad = 0.4,
                  CNV_cutoff = 1)

    })

  })
}

## To be copied in the UI
# mod_overview_ui("overview_1")

## To be copied in the server
# mod_overview_server("overview_1")
