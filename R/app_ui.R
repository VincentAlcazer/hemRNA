#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
#' @import shinyFiles
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic

    dashboardPage(

      header = dashboardHeader(
        title = "hemRNA v0.4.3"
      ),

      #################### ==================== SIDEBAR ====================  ####################

      sidebar = dashboardSidebar(
        shinydashboard::sidebarMenu(
          # Setting id makes input$tabs give the tabName of currently-selected tab
          id = "tabs",
          menuItem("Home", tabName = "home", icon = icon("play-circle")),
          menuItem("Data loading", tabName = "data", icon = icon("spinner")),
          menuItem("Expression", tabName = "expression", icon = icon("poll")),
          menuItem("Fusion", tabName = "fusion", icon = icon("exchange-alt")),
          menuItem("Variants", tabName = "variants",icon = icon("bezier-curve"),
                   menuSubItem("Hotspot", tabName = "hotspot"),
                   menuSubItem("Variant calling", tabName = "variant")
                   ),
          menuItem("CNV", tabName = "CNV", icon = icon("copy"))

          )
      ), # sidebar


      #################### ==================== BODY ====================  ####################

      body = dashboardBody(
        tabItems(

          ##### ===== Home

          tabItem("home", mod_home_ui("home_ui_1")),

          ##### ===== Data loading

          tabItem("data",
                  tabsetPanel(
                    id = "data", type = "tabs",
                    tabPanel("Full results loading",
                             fluidPage(
                               box(title = "Folder", width = 6,
                                   ### === Load results folder
                                   HTML("<b>Select results folder</b></br>"),
                                   shinyDirButton('folder_path', label='Browse...', title='Please select a folder'),
                                   div(style = "margin-top: 20px"),

                                   ### === Load bed
                                   selectInput("panel", label = "Select a bed panel",
                                               choices = c("None"),
                                               selected = "None"),
                                   div(style = "margin-top: -20px"),
                                   p("Panel should be added to the data/bed_panels folder as bed files."),

                                   div(style = "margin-top: 20px"),

                                   h3("File check:"),

                                   htmlOutput("salmon_check"),
                                   shinycssloaders::withSpinner(htmlOutput("salmon_dim"),type = 1),
                                   shinycssloaders::withSpinner(htmlOutput("panel_dim"),type = 1),

                                   htmlOutput("fusion_check"),
                                   htmlOutput("hotspot_check"),
                                   htmlOutput("cnv_check")

                                    # p("Preview"),
                                    # downloadButton("download_table", "Download table (.tsv)"),
                                    # shinycssloaders::withSpinner(DT::DTOutput("preview_data"),type = 6)

                               ),# box
                               box(title = "Meta", width = 6,
                                   p("/!\ La première colonne doit correspondre à la clé d'identification des échantillons/patients."),
                                   fileInput("meta",
                                             label = ("Meta"),
                                             accept = c(
                                               "text/tab-separated-values",
                                               "text/comma-separated-values",
                                               "text/plain",
                                               "text/csv",
                                               ".csv",
                                               ".tsv",
                                               ".xls",
                                               ".xlsx"
                                             )
                                   ),
                                   div(style = "margin-top: -20px"),
                                   radioButtons("sep", "Separator",
                                                choices = c(
                                                  "comma-delim (.csv1)" = ",",
                                                  "semi colon-delim (.csv2)" = ";",
                                                  "tab-delim (.tsv/.txt)" = "\t",
                                                  "Excel (.xls/.xlsx)" = "xl"
                                                ),
                                                selected = "\t"
                                   ),
                                   div(style = "margin-top: -20px"),
                                   radioButtons("dec", "Decimal",
                                                choices = c(
                                                  "Comma (,)" = ",",
                                                  "Period (.)" = "."
                                                ),
                                                selected = ","
                                   ),
                                   p("Preview"),
                                   shinycssloaders::withSpinner(DT::DTOutput("preview_meta"),type = 6)
                               ), #box
                            box(title = "Missing genes", width = 6,
                                p("Genes/Features present in panel but not found in exp. dataset:"),
                                htmlOutput("missing_genes")
                                )



                             ) # fluid page

                    ),# tabPanel

                    tabPanel("Individual loading",
                             fluidPage(
                               box(title = "Data", width = 6,
                                   HTML("<b>Select salmon folder</b></br>"),
                                   shinyDirButton('df_indiv_xp', label='Browse...', title='Please select a folder'),
                                   div(style = "margin-top: 20px"),
                                   fileInput("panel_indiv_xp", label = "Upload used panel (optional)"),
                                   div(style = "margin-top: -20px"),
                                   p("Panel should be uploaded as a  simple .txt file
                                     with a gene_id column corresponding to gene names"),


                                   p("Preview"),
                                   shinycssloaders::withSpinner(DT::DTOutput("preview_indiv_data"),type = 6)

                               ),# box
                              box(title = "Meta", width = 6,
                                  p("/!\ La première colonne doit correspondre à la clé d'identification des échantillons/patients."),
                                    fileInput("meta_indiv",
                                              label = ("Meta"),
                                              accept = c(
                                                "text/tab-separated-values",
                                                "text/comma-separated-values",
                                                "text/plain",
                                                "text/csv",
                                                ".csv",
                                                ".tsv",
                                                ".xls",
                                                ".xlsx"
                                              )
                                    ),
                                    div(style = "margin-top: -20px"),
                                    radioButtons("sep_indiv", "Separator",
                                                 choices = c(
                                                   "comma-delim (.csv1)" = ",",
                                                   "semi colon-delim (.csv2)" = ";",
                                                   "tab-delim (.tsv/.txt)" = "\t",
                                                   "Excel (.xls/.xlsx)" = "xl"
                                                 ),
                                                 selected = "\t"
                                    ),
                                    div(style = "margin-top: -20px"),
                                    radioButtons("dec_indiv", "Decimal",
                                                 choices = c(
                                                   "Comma (,)" = ",",
                                                   "Period (.)" = "."
                                                 ),
                                                 selected = ","
                                    ),
                                  p("Preview"),
                                  shinycssloaders::withSpinner(DT::DTOutput("preview_indiv_meta"),type = 6)
                              ) #box


                             ) # fluid page

                    )# tabPanel
                  )# tabSetPanel
                  ), # Tab Item

          ##### ===== Expression

          tabItem("expression",
                  tabsetPanel(
                    id = "expression", type = "tabs",
                    tabPanel("PCA",  mod_expression_ui("expression_ui_1")),
                    tabPanel("Individual xp",mod_expression_individual_ui("expression_individual_ui_1")),
                    tabPanel("Signatures", mod_expression_signatures_ui("expression_signatures_ui_1")),
                    #tabPanel("Heatmap",  mod_expression_heatmap_ui("expression_heatmap_ui_1")),
                    tabPanel("DESEQ2", mod_expression_deseq_ui("expression_deseq_ui_1"))

                  )#TabSetPanel

        ), #tabItem
        tabItem("fusion",mod_fusion_ui("fusion_ui_1")),

        tabItem("hotspot",mod_hotspot_ui("hotspot_ui_1")),

        tabItem("variant",
                h1("Variant calling"),
                HTML("Methods: bam were generated according to the GATK good practice pipeline,
                including: <br/>
                - Fastq reads alignement against hg19 reference genome with 1000 genome decoy using STAR (gencode v19 annotation) <br/>
                - Bam processing with markduplicates, SplitNCigarReads and baserecalibration (BQSR) <br/>
                BQSR corrected BAM were then processed with 3 different caller: <br/>
               - GATK HaplotypeCaller with default parameters and -dont-use-soft-clipped-bases -dbsnp parameters with the dbsnp138 reference.
               Basic filter was then appllied (--window 35, --cluster 3, FS > 30, QD < 2) <br/>
               - Freebayes with default parameters and --use-duplicate-reads --standard-filters <br/>
               - Samtools mpileup with -Q 30
            "),


                tabsetPanel(
                  id = "variant_sub", type = "tabs",
                  # tabPanel("Overall",
                  #          mod_variant_ui("variant_ui_1")
                  # ),#tabsetpanel
                  tabPanel("GATK",
                           mod_variant_GATK_ui("variant_GATK_ui_1")
                           )
                )##tabsetpanel
        ),

        tabItem("CNV",mod_CNV_ui("CNV_ui_1"))

        )#tabItems
      ) #body

    ) #dashboard page



  ) #taglist
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){

  add_resource_path(
    'www', app_sys('app/www')
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'hemRNA'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}

