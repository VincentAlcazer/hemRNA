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

      column(8,
           tabsetPanel(
             id = "DESEQ2", type = "tabs",
             tabPanel("Parameters",
                      h1("DGE analysis parameters"),
                      column(
                        6,
                        wellPanel(
                          h3("DESEQ2 parameters"),
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
                          actionButton(ns("run_analysis"), "Run")
                        )
                      ),#Column
                      column(
                        6,
                        wellPanel(

                          selectInput(ns("test"),
                                      label = ("Test"),
                                      multiple = F, selected = NULL,
                                      ""
                          )


                        )
                      ),#Column
              column(12,
                     h3("Results table"),
                     p("According to the number of patients/genes, DESEQ2 analysis can take a long time (~5min for 10 patients/60 000 genes)"),
                     shinycssloaders::withSpinner(DT::DTOutput(ns("preview")),type = 6)
                     )

             ),#tabPanel
             tabPanel("Volcano plot",
                      h1("Volcano plot"),
                      shinycssloaders::withSpinner(plotOutput(ns("volcanoplot")),type=6)
                      )
           )#tabsetpanel
    ), #column
    column(1),
    column(2,
           absolutePanel(
             width = 200, right = 20, draggable = T,
             style = "opacity: 0.85",
             wellPanel(
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
    df_filt <- reactive({r$test$df_filt})

    meta <- reactive({r$test$meta})

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
        "covariates",
        choices = colnames(meta()[-1])
      )
    })

    ##### DESEQ2

    deseq2_res_prod <- reactive({

      req(meta())
      req(input$run_analysis >= 1)
      isolate({

        meta_deseq <- meta() %>% column_to_rownames("patient_id")

        intersect <- intersect(colnames(r$test$txi[[1]]),rownames(meta_deseq))

        meta_deseq <- meta_deseq[intersect,]
        meta_deseq[,input$group] <- factor(meta_deseq[,input$group])

        if(length(intersect) < 2){
          message("<2 common patients found: check meta")
        }

        if(is.null(input$covariates)){
          design <- as.formula(paste0("~", input$group))
        } else {
          design <- as.formula(paste0("~",paste(input$covariates, collapse = "+")," + ",paste0(input$group)))
        }

        dds <- DESeqDataSetFromTximport(r$test$txi, meta_deseq, design =  design)

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

    })

    deseq2_res <- reactive({
      data <- read.table("data/2021_12_15_DESEQ2_res_dev.tsv", sep = "\t", stringsAsFactors = F, header = T)
    })

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

    ##### ===== Volcano

    volcano <- reactive({

       padj.cutoff = 0.05
      lfc.cutoff = 0.58

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
      volcano_df$genelabels[1:10] <- volcano_df$gene_name[1:10]

      #Volcano plot
      plot <- volcano_df %>% na.omit %>%
        ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(aes(color = color, size = color)) +
        geom_text_repel(aes(label = genelabels)) +
        geom_text(data = df_label, aes(x=x_coord, y=y_coord, label = gene_n), size = 10) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 1) +
        geom_vline(xintercept = 0.58, linetype = "dashed", size = 1) +
        geom_vline(xintercept = -0.58, linetype = "dashed", size = 1) +
        scale_color_manual(values=color_vector) +
        scale_size_manual(values = c(3,2,3)) +
        labs(title=paste0(""),
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

  })
}

## To be copied in the UI
# mod_expression_deseq_ui("expression_deseq_ui_1")

## To be copied in the server
# mod_expression_deseq_server("expression_deseq_ui_1")
