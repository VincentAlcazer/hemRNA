#' expression_deseq UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_expression_deseq_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(

      column(10,
           tabsetPanel(
             id = "DESEQ2", type = "tabs",
             tabPanel("Table",
                      column(
                        6,
                        wellPanel(
                          h3("Run new analysis..."),
                          selectInput(ns("group"),
                                      label = ("Groups to compare"),
                                      multiple = F, selected = NULL,
                                      ""
                          ),
                          selectInput(ns("covariates"),
                                      label = ("Covariates to control"),
                                      multiple = T, selected = NULL,
                                      ""
                          ),
                          #div(style = "margin-top: -20px"),
                          actionButton(ns("run_analysis"), "Run"),
                          div(style = "margin-top: 10px"),
                        )
                      ),#Column
                      column(
                        6,
                        wellPanel(

                          h3("... OR Load results from previous run"),
                          fileInput(ns("deseq_results_tsv"),
                                    label = ("DESEQ2_res.tsv"),
                                    accept = c(
                                      "text/tab-separated-values",
                                      ".tsv"
                                    )
                          )

                        )
                      ),#Column

              column(12,
                     h3("DESEQ2 Results"),
                     p("According to the number of patients/genes, DESEQ2 analysis can take a long time (~5min for 10 patients/60 000 genes)"),
                     downloadButton(ns("download_table"), "Download table (.tsv)"),
                     shinycssloaders::withSpinner(DT::DTOutput(ns("preview")),type = 6)
                     )

             ),#tabPanel
             tabPanel("Volcano plot",
                      shinycssloaders::withSpinner(plotOutput(ns("volcanoplot"), height = 800),type=6)
                      ),
             tabPanel("Heatmap",
                      shinycssloaders::withSpinner(plotOutput(ns("heatmap_plot"), height = 800),type=6),
                      column(
                        6,
                        wellPanel(
                          h3("Clustering parameters"),
                          selectInput(ns("metric"),
                                      label = ("Distance"),
                                      multiple = F, selected = "euclidean",
                                      choices = c("euclidean","pearson", "maximum","minkowski","binary")
                          ),
                          selectInput(ns("method"),
                                      label = ("Method"),
                                      multiple = F, selected = "ward.D2",
                                      choices = c("ward.D2", "single", "complete", "average","mcquitty","median","centroid")
                          ),
                          checkboxInput(ns("cluster_rows"),
                                        label = "Cluster rows",
                                        value = T),
                          checkboxInput(ns("cluster_cols"),
                                        label = "Cluster cols",
                                        value = T),
                        )
                      ),#Column
                      column(
                        6,
                        wellPanel(
                          h3("Heatmap parameters"),
                          selectInput(ns("top_anno"),
                                      label = ("Top annotation"),
                                      multiple = F, selected = NULL,
                                      ""
                          ),
                          checkboxInput(ns("row_names"),
                                        label = "Show row names",
                                        value = T),
                          sliderInput(ns("row_fontsize"),
                                      "Row fontsize",
                                      value = 12, min = 1, max = 30),
                          checkboxInput(ns("col_names"),
                                        label = "Show col names",
                                        value = T),
                          sliderInput(ns("col_fontsize"),
                                      "Column fontsize",
                                      value = 12, min = 1, max = 30)

                        )
                      )#Column

             )

           )#tabsetpanel
    ), #column
    column(1),
    column(2,
           absolutePanel(
             width = 200, right = 20, draggable = T,
             style = "opacity: 0.85",
             wellPanel(
               numericInput(ns("padj"),
                            label = ("FDR (q-value) cutoff"),
                            value = 0.05, max = 1
                            ),
               numericInput(ns("FC"),
                            label = ("Fold-change cutoff"),
                            value = 2
               ),
               numericInput(ns("top_genes"),
                            label = ("Top genes to name (Volcano)"),
                            value = 10
               ),
               # sliderInput(ns("gene_size"),
               #              label = ("Gene name size (Volcano)"),
               #              value = 10, min=1, max=30, step = 1
               # ),
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

#' expression_deseq Server Functions
#'
#' @noRd
mod_expression_deseq_server <- function(id, r){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    ##### Datasets
    meta <- reactive({r$test$meta})
    df_filt <- reactive({r$test$df_filt})

    folder_path <- reactive({r$test$folder_path})

    txi <- reactive({
      req(folder_path())
      req(meta())
      req(input$run_analysis >= 1)
      isolate({

      message("Loading salmon txi files")
      files <- list.files(paste0(folder_path(),"/salmon"), pattern = "quant.sf", recursive = T, full.names = T)

      if(length(files)==0 | is.null(files)){
        message("No expression file detected")
      } else {
        names(files) <- gsub(".salmon/quant.sf","", list.files(paste0(folder_path(),"/salmon"), pattern = "quant.sf", recursive = T))

        dat <- tximport::tximport(files[meta()$patient_id], type = "salmon", tx2gene = tx2gene, txOut = F)

        message("Salmon txi successfully loaded")
        return(dat)

      }
      })
    })

    ##### Variables
    observe({
      req(meta())
      updateSelectInput(
        session,
        "group",
        choices = colnames(meta()[-1])
      )
    })

    observe({
      req(meta())
      updateSelectInput(
        session,
        "top_anno",
        choices = colnames(meta()[-1])
      )
    })

    observe({
      req(meta())
      updateSelectInput(
        session,
        "covariates",
        choices = colnames(meta()[-1])
      )
    })

    ##### DESEQ2

    deseq2_res <- reactive({

      if(!is.null(input$deseq_results_tsv$datapath)){

        res <- data.table::fread(input$deseq_results_tsv$datapath, stringsAsFactors = F, data.table = F)
        return(res)

      } else {

      req(meta())
      req(input$run_analysis >= 1)
      isolate({


          meta_deseq <- meta()
          rownames(meta_deseq) <- meta()$patient_id

          intersect <- intersect(colnames(txi()[[1]]),rownames(meta_deseq))

          meta_deseq <- meta_deseq[intersect,]
          meta_deseq[,input$group] <- factor(meta_deseq[,input$group])

          txi_filt <- txi()

          if(length(intersect) < 2){
            message("<2 common patients found: check meta")
          }

          if(is.null(input$covariates)){
            design <- as.formula(paste0("~", input$group))
          } else {
            design <- as.formula(paste0("~",paste(input$covariates, collapse = "+")," + ",paste0(input$group)))
          }

          dds <- DESeq2::DESeqDataSetFromTximport(txi(), meta_deseq, design =  design)

          dds <- DESeq2::DESeq(dds)

          coef_n <- grep(input$group, DESeq2::resultsNames(dds))

          res <- DESeq2::lfcShrink(dds, coef=coef_n, type="ashr") %>%
            data.frame(check.names = F) %>%
            rownames_to_column("gene_id") %>%
            arrange(padj) %>%
            mutate_at(c("baseMean","log2FoldChange","lfcSE"), function(x){round(x,2)}) %>%
            mutate_at(c("pvalue","padj"), function(x){as.numeric(format.pval(x, digits = 1))}) %>%
            left_join(gene_anno) %>% filter(is.na(padj) == F)

          return(res)


      })
}
    })

    # deseq2_res <- reactive({
    #   data <- read.table("data/2021_12_15_DESEQ2_res_dev.tsv", sep = "\t", stringsAsFactors = F, header = T)
    # })

    output$preview <- DT::renderDT(
      deseq2_res(), # data
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
        paste("DESEQ2_res.tsv")
      },
      content = function(file) {
        write.table(deseq2_res(), file, row.names = FALSE, sep = "\t", quote = F)
      }
    )

    ##### ===== Volcano

    volcano <- reactive({

      padj.cutoff = input$padj
      lfc.cutoff = log2(input$FC)

      req(deseq2_res())

      volcano_df <- deseq2_res() %>%
        mutate(color = case_when(
          padj < padj.cutoff & log2FoldChange >= lfc.cutoff ~ "Upregulated",
          padj < padj.cutoff & log2FoldChange <= -lfc.cutoff ~ "Downregulated",
          padj > padj.cutoff  ~ "NS"

        ),
        genelabels="") %>%
        arrange(padj)

      n_OE <- length(which(volcano_df$log2FoldChange >= lfc.cutoff & volcano_df$padj < padj.cutoff))
      n_DE <- length(which(volcano_df$log2FoldChange <= -lfc.cutoff & volcano_df$padj < padj.cutoff))


      df_label <- data.frame(
        gene_n = c(n_DE,n_OE),
        assoc_type = c("DownExpressed","Overexpressed"),
        x_coord = c(min(volcano_df$log2FoldChange, na.rm=T)*0.8,max(volcano_df$log2FoldChange, na.rm=T)*0.8),
        y_coord = c(-log10(quantile(volcano_df$padj,1e-05, na.rm=T)),-log10(quantile(volcano_df$padj,1e-05, na.rm=T))))


      # color vector
      color_vector <- c("#9c1c1c","grey","#006aaf")
      names(color_vector) <- c("Upregulated","NS","Downregulated")


      ## Create a column to indicate which genes to label (top 10)
      if(input$top_genes > 0){
        volcano_df$genelabels[1:input$top_genes] <- volcano_df$gene_name[1:input$top_genes]
      }


      #Volcano plot
      plot <- volcano_df %>% na.omit %>%
        ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(aes(color = color), size = 2) +
        geom_text_repel(aes(label = genelabels), size = 6) +
        geom_text(data = df_label, aes(x=x_coord, y=y_coord, label = gene_n), size = 10) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 1) +
        geom_vline(xintercept = 0.58, linetype = "dashed", size = 1) +
        geom_vline(xintercept = -0.58, linetype = "dashed", size = 1) +
        scale_color_manual(values=color_vector) +
        #scale_size_manual(values = c(3,2,3)) +
        labs(title=paste0(""), color = "Sig",
             y="-log10 (Adjusted p-value)",
             x="log2 Fold Change (ashr)")

      return(plot)


    })

    output$volcanoplot <- renderPlot({
      volcano() +
        default_theme +
        theme(
          legend.position = input$legend_ext)
    })


    ##### ===== Heatmap

    heatmap <- reactive({

      req(deseq2_res())
      req(df_filt())
      req(meta())

      features <- deseq2_res() %>% filter(padj < input$padj & abs(log2FoldChange) >= input$FC) %>%
        distinct(gene_name) %>% unlist() %>% as.character()

      data_heat <-  df_filt() %>%
        filter(gene_name %in% features) %>%
        column_to_rownames("gene_name") %>% t()

      samples <- intersect(unique(meta()$patient_id),rownames(data_heat))

      data_heat <- data_heat[samples,]

      meta_heat <- meta()
      rownames(meta_heat) <- meta_heat$patient_id
      meta_heat <- meta_heat[samples,]

      plot <- fast_map(data_heat,
                       center_scale = T,
                       heat_colors="ATAC",
                       anno_pal = viridis_pal(),
                       Zlim = 3, method = input$method, metric = input$metric, title = "",
                       Data_input = "Z-score expression",
                       anno_1 = meta_heat[,input$top_anno], anno_1_name = input$top_anno, anno_1_cols = NULL,
                       anno_2 = NULL, anno_2_name = NULL, anno_2_cols = NULL,
                       cluster_rows = input$cluster_rows, cluster_cols = input$cluster_cols,
                       reorder_rows = input$cluster_rows, reorder_cols = input$cluster_cols,
                       row_names = input$row_names, col_names = input$col_names,
                       row_fontsize = input$row_fontsize, col_fontsize = input$col_fontsize,
                       show_hclust_n = 0, show_pearson = F)

      return(plot)

    })

    output$heatmap_plot <- renderPlot({
      heatmap()
    })


  })
}

## To be copied in the UI
# mod_expression_deseq_ui("expression_deseq_ui_1")

## To be copied in the server
# mod_expression_deseq_server("expression_deseq_ui_1")
