#' expression UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_expression_ui <- function(id){
  ns <- NS(id)
  tagList(

    fluidPage(
      br(),
      column(8,
             shinycssloaders::withSpinner(plotOutput(ns("pca")),type=6)
             ),
      column(1),
      column(2,
             absolutePanel(
                     width = 200, right = 20, draggable = T,
                     style = "opacity: 0.85",
                     wellPanel(
                       sliderInput(ns("var_features"),
                                    label = ("Top variable genes to show (%)"),
                                    min = 1, max = 100, step = 1, value = 25

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

#' expression Server Functions
#'
#' @noRd
mod_expression_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    df_filt <- reactive({r$test$df_filt})

    pca_plot <- reactive({
      req(df_filt())
      ##  1. Keep only features expressed in at least 5% of cells patients
      data_acp <- df_filt() %>% column_to_rownames("gene_name") %>% t()

      if(round(nrow(data_acp)*0.05) == 0){
        min_pts = 1
      } else {
        min_pts = round(nrow(data_acp)*0.05)
      }

      data_acp <- data_acp[, colSums(data_acp > 0) >= min_pts ]

      ##  2. Top 25% variables regions
      vars <- apply(data_acp, 2, var)
      top_var <- names(vars[order(vars, decreasing = T)][1:round(length(vars)*input$var_features/100)])
      data_acp <- data_acp[,top_var]


      ##Run PCA
      pca <- prcomp(data_acp, scale=T, center = T)

      ##Scree plot
      p <- factoextra::fviz_pca_ind(pca,
                               title = "PCA - Individuals",
                               legend.title = "Sample",
                               mean.point = F, #pt du centre de gravitĂ©
                               addEllipses = F,
                               ellipse.type = "confidence",
                               ellipse.level = 0.95,
                               pointsize = 4,
                               pointshape = 21,
                               col.ind = "black",
                               fill.ind = rownames(data_acp), # Colorer par groupes
                               repel = T,
                               geom.ind=c("point")
      )

      return(p)

    })

    output$pca <- renderPlot({
      pca_plot() +
        default_theme +
        theme(
          legend.position = input$legend_ext
        )
    })




  })
  }


## To be copied in the UI
# mod_expression_ui("expression_ui_1")

## To be copied in the server
# mod_expression_server("expression_ui_1")
