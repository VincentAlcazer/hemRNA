#' variant_GATK UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_variant_GATK_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      tabsetPanel(
        id = "variant", type = "tabs",
        tabPanel("Graph",

                 column(8, shinycssloaders::withSpinner(plotOutput(ns("tileplot")),type=6)
                 ),
                 column(1)
        ), #tabPanel

        tabPanel("Result table",
                 column(8, shinycssloaders::withSpinner(DT::DTOutput(ns("result_table")),type=6)
                 ),
                 column(1)
        )
      ),##tabsetpanel

      column(2,
             absolutePanel(
               width = 200, right = 20, draggable = T,
               style = "opacity: 0.85",
               wellPanel(
                 checkboxInput(ns("DB"), label = "dbSNP", value = TRUE),
                 sliderInput(ns("QD"), label = "Min. QD",
                             min = 1, max = 60, value = 2, step= 1),
                 sliderInput(ns("FS"), label = "Max. FS",
                             min = 0, max = 100, value = 60, step= 1),
                 sliderInput(ns("SOR"), label = "Max. SOR",
                             min = 0, max = 10, value = 3, step= 1),
                 sliderInput(ns("MQ"), label = "Min. MQ",
                             min = 10, max = 80, value = 40, step= 1),
                 sliderInput(ns("MQRankSum"), label = "Min. MQRankSum",
                             min = -20, max = 10, value = -13, step= 0.5),
                 sliderInput(ns("ReadPosRankSum"), label = "Min. ReadPosRankSum",
                             min = -20, max = 10, value = -8, step= 0.5),

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

#' variant_GATK Server Functions
#'
#' @noRd
mod_variant_GATK_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    vcf_path <- reactive({r$test$GATK_path})

    GATK_vcf <- reactive({

      req(vcf_path())

      vcf <- vcfR::vcfR2tidy(vcfR::read.vcfR(r$test$GATK_path, verbose = FALSE ))$fix

      return(vcf)

    })

    ## Update Parameters
    observe({
        req(GATK_vcf())
        updateSliderInput(
          session,
          "QD",
          max = round(max(GATK_vcf()$QD, na.rm=T))
      )

        updateSliderInput(
          session,
          "FS",
          max = round(max(GATK_vcf()$FS, na.rm=T))
        )
        updateSliderInput(
          session,
          "SOR",
          max = round(max(GATK_vcf()$SOR, na.rm=T))
        )
        updateSliderInput(
          session,
          "MQ",
          max = round(max(GATK_vcf()$MQ, na.rm=T))
        )
        updateSliderInput(
          session,
          "MQRankSum",
          max = round(max(GATK_vcf()$MQRankSum, na.rm=T))
        )
        updateSliderInput(
          session,
          "ReadPosRankSum",
          max = round(max(GATK_vcf()$ReadPosRankSum, na.rm=T))
        )
    })


    filt_vcf <- reactive({

      dat <- GATK_vcf() %>%
        filter(DB == input$DB,
               QD > input$QD,
               FS < input$FS,
               SOR < input$SOR,
               MQ > input$MQ,
               MQRankSum > input$MQRankSum,
               ReadPosRankSum > input$ReadPosRankSum)

      return(dat)
    })


    output$result_table <- DT::renderDT(
      filt_vcf(), # data
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
        paste("Hotspot.tsv")
      },
      content = function(file) {
        write.table(filt_vcf(), file, row.names = FALSE, sep = "\t", quote = F)
      }
    )





  })
}

## To be copied in the UI
# mod_variant_GATK_ui("variant_GATK_ui_1")

## To be copied in the server
# mod_variant_GATK_server("variant_GATK_ui_1")
