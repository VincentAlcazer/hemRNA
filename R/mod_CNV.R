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
      tabsetPanel(
        id = "variant", type = "tabs",
        tabPanel("Individual graph",
                 p("Methods: Mean norm count per cytoband are calculated for each patient.
                    The mean norm count per patient is then divided by the mean norm count of the overall cohort
                    per cytoband.
                    The log2(ratio) is represented.
                    Important variations (fold change between -4 and +4)are outlined with dashed lines.
                   "),
                 # column(5,
                 #        radioButtons(ns("count_type"),
                 #                     label = c("Count"),
                 #                     choices = c("Raw counts"="tot_raw_count",
                 #                                 "Norm counts"="tot_norm_count"),
                 #                     selected = "tot_raw_count")
                 #        ),
                 # column(5,
                 #        sliderInput(ns("y_size"), label = "y-axis font size",
                 #                    min = 1, max = 30, value = 12, step= 1),
                 #        sliderInput(ns("x_size"), label = "x-axis font size",
                 #                    min = 1, max = 30, value = 14, step= 1)
                 #
                 #        ),

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
                 radioButtons(ns("count_type"),
                              label = c("Norm. count"),
                              choices = c("TPM","CPM"),
                              selected = "TPM"),
                 selectizeInput(ns("sample"),
                                label = ("Sample"),
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

    bed <- reactive({r$test$bed})
    ftcount <- reactive({r$test$ftcount})


    observe({
      updateSelectizeInput(
        session,
        "sample",
        choices = c(unique(gsub("_(T|C)PM","", colnames(ftcount())[-c(1:6)]))),
        selected = c(unique(gsub("_(T|C)PM","", colnames(ftcount())[-c(1:6)])))[1],
        server = T
      )
    })



    bed_mod <- reactive({

      req(bed())
      req(hg19_cytoband)

      bed_mod <- bed() %>%
        group_by(chr, gene_name) %>%
        summarise(start = min(start),
                  end = max(end)) %>%
        full_join(hg19_cytoband, by = c("chr"), suffix = c("_bed","_cytob")) %>%
        mutate(inside = if_else(start_bed >= start_cytob & end_bed <= end_cytob, TRUE, FALSE)) %>%
        filter(inside == T) %>%
        mutate(chr_n = gsub("chr","",chr)) %>%
        mutate(chr_n = gsub("X","23",chr_n),
               chr_n = gsub("Y","24",chr_n)) %>%
        mutate(chr_n = as.numeric(chr_n),
               chr_cytob = paste0(chr_n,cytoband)) %>%
        mutate(chr_cytob2 = gsub("\\.[0-9]+$","",chr_cytob)) %>%
        mutate(chr_cytob3 = gsub("[0-9]+$","",chr_cytob2)) %>%
        arrange(chr_n, start_cytob)

      bed_mod$chr_cytob <- factor(bed_mod$chr_cytob, levels = unique(bed_mod$chr_cytob))
      bed_mod$chr_cytob2 <- factor(bed_mod$chr_cytob2, levels = unique(bed_mod$chr_cytob2))
      bed_mod$chr_cytob3 <- factor(bed_mod$chr_cytob3, levels = unique(bed_mod$chr_cytob3))

      return(bed_mod)

    })

    df_CNV <- reactive({

      req(bed_mod())
      req(ftcount())

     dat <- ftcount() %>%
        left_join(bed_mod(), by = c("Geneid"="gene_name")) %>%
        mutate(order = seq(1:nrow(.))) %>%
        select(order, chr_n,chr_cytob, chr_cytob2,chr_cytob3, Geneid, ends_with(input$count_type)) %>%
        pivot_longer(ends_with(input$count_type), names_to = "sample_id", values_to = "norm_count") %>%
        mutate(log2_norm_count = log2(norm_count+1)) %>%
        group_by(chr_cytob2,sample_id) %>%
        mutate(mean_cytoband = mean(log2_norm_count)) %>%
        ungroup() %>%
        group_by(chr_cytob2) %>%
        mutate(n = n(),
               cohort_mean = mean(log2_norm_count, na.rm=T),
               cohort_sd = sd(log2_norm_count, na.rm=T)) %>%
        mutate(cohort_se = cohort_sd / sqrt(n)) %>%
        mutate(CI95 = cohort_se * qt((1 - 0.05) / 2 + 0.5, n - 1)) %>%
        mutate(CI95_sup = cohort_mean + CI95,
               CI95_inf = cohort_mean - CI95) %>%
        ungroup() %>%
        group_by(chr_cytob2,sample_id) %>%
        mutate(ratio = mean(log2_norm_count) / cohort_mean,
               ratio_inf = mean(log2_norm_count / CI95_sup),
               ratio_sup = mean(log2_norm_count / CI95_inf)
        ) %>%
        ungroup()

      return(dat)


    })





    ##### ===== Plots

    plot <- reactive({

      req(df_CNV())

      x_breaks <- df_CNV() %>%
        group_by(chr_cytob3) %>%
        summarise(mean_order = mean(order)) %>%
        na.omit()

      dat <- df_CNV() %>%
        filter(sample_id == paste0(input$sample,"_",input$count_type)) %>%
        na.omit() %>%
        ggplot(aes(x = order, y = log2(ratio), color = log2(ratio))) +
        #geom_point(alpha = 0.1, color = "grey") +
        geom_line(size = 1.5) +
        labs(x = "") +
        geom_hline(yintercept =  2, color = "black", linetype = 2) +
        geom_hline(yintercept =  - 2, color = "black", linetype = 2) +
        # geom_line(aes(x = order, y = log10(ratio_inf)), color = "black", linetype = 2) +
        # geom_line(aes(x = order, y = log10(ratio_sup)), color = "black", linetype = 2) +
        scale_x_continuous(breaks =x_breaks$mean_order,  labels = x_breaks$chr_cytob3)+
        scale_color_viridis_c() +
        coord_flip()

      return(dat)

    })

    output$lineplot <- renderPlot(height = 2400,{
      plot() +
        default_theme +
        theme(
          legend.position = input$legend_ext,
          axis.text.y = element_text(size = input$y_size, face = "bold"),
          axis.text.x = element_text(size = input$x_size, face = "bold")
        )
    })


    output$preview_data <- DT::renderDT(
      df_CNV(), # data
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
        write.table(df_CNV(), file, row.names = FALSE, sep = "\t", quote = F)
      }
    )


  })
}

## To be copied in the UI
# mod_CNV_ui("CNV_ui_1")

## To be copied in the server
# mod_CNV_server("CNV_ui_1")
