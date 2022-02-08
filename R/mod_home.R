#' home UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_home_ui <- function(id){
  ns <- NS(id)
  tagList(

    fluidPage(
      column(10,

             h1("Welcome to hemRNA"),

             p("hemRNA is a free open-source  software designed for routine RNA-sequencing analysis in Hematology.
             It has been developed with the R software, using the", a("Shiny", href = "http://shiny.rstudio.com"), "package."),

             h2("Change log"),

             HTML("
            08-02-2022: v0.4.1 - Stability improvement & bug fixes <br/>
            26-01-2022: v0.4 - CNV module & hotspot module improvement <br/>
            19-01-2022: v0.3 - Variant calling GATK module <br/>
            20-12-2021: v0.2 - Fusion & hotspot modules <br/>
            01-12-2021: v0.1 - Data loading & Expression modules <br/>
            ")


             ) #column



    ) #fluidPage

  )
}

#' home Server Functions
#'
#' @noRd
mod_home_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

  })
}

## To be copied in the UI
# mod_home_ui("home_ui_1")

## To be copied in the server
# mod_home_server("home_ui_1")
