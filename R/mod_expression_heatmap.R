#' expression_heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_expression_heatmap_ui <- function(id, r){
  ns <- NS(id)
  tagList(

    fluidPage(
      h1("Heatmap"),

      column(8,
             shinycssloaders::withSpinner(plotOutput(ns("heatmap")),type=6)
      ),
      column(1),
      column(2,
             absolutePanel(
               width = 200, right = 20, draggable = T,
               style = "opacity: 0.85",
               wellPanel(
                 # selectInput(ns("Group"),
                 #             label = ("Groups"),
                 #             multiple = F, selected = NULL,
                 #             ""
                 # ),
                 selectInput(ns("Genes"),
                             label = ("Genes"),
                             multiple = T, selected = NULL,
                             ""
                 ),
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

#' expression_heatmap Server Functions
#'
#' @noRd
mod_expression_heatmap_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    ##### Datasets
    df_filt <- reactive({r$test$df_filt})

    # data_heat <- reactive({
    #   df_filt %>% filter(gene_name %in% input$Genes) %>% column_to_rownames("gene_name")
    # })

    ##### Variables
    # observe({
    #   req(df_filt())
    #   updateSelectInput(
    #     session,
    #     inputId = "Genes",
    #     choices = df_filt() %>% distinct(gene_name) %>% unlist() %>% as.character(),
    #
    #   )
    # })

#
#     plot <- reactive({
#
#       pheatmap(data_heat(),
#                cluster_cols = T,
#                cluster_rows = T,
#                clustering_distance_rows = "euclidean",
#                clustering_distance_cols = "euclidean",
#                clustering_method = "ward.D2")
#
#     })
#
#     output$heatmap <- renderPlot({
#       plot()
#
#     })


  })
}

## To be copied in the UI
# mod_expression_heatmap_ui("expression_heatmap_ui_1")

## To be copied in the server
# mod_expression_heatmap_server("expression_heatmap_ui_1")
