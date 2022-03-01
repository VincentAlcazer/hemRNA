#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function( input, output, session ) {
  # Your application server logic


  #################### ==================== Dataframe loading  ====================  ####################

  r <- reactiveValues(
    test = reactiveValues()
  )

  ##### ===== Full results folder

  # ## For volumes parse
  volumes <- getVolumes()
  shinyDirChoose(input, 'folder_path', root=volumes, session=session)
  folder_path <- reactive({
    return(print(parseDirPath(volumes, input$folder_path)))
  })

  # Parse local folder
  # shinyDirChoose(input, 'folder_path', root=c(root='../'), session=session)
  #
  # folder_path <- reactive({
  #   return(print(parseDirPath(c(root='../'), input$folder_path)))
  # })

  ##### ===== Individual folders

  # Parse salmon folder
  # shinyDirChoose(input, 'df_indiv_xp', root=c(root='../'), session=session)
  #
  # path1 <- reactive({
  #   return(print(parseDirPath(c(root='../'), input$df_indiv_xp)))
  # })

  ##### ===== Signatures list

  sig_list <- reactive({

    sig_dir <- system.file("extdata","signatures", package = "hemRNA")
    sig_files <- list.files(sig_dir)

    sig_list <- list()

    for(sig in sig_files){

      name <- gsub(".txt$","",sig)

      sig_list[[name]] <- read.table(paste0(sig_dir,"/",sig), sep = "\t", stringsAsFactors = F, header = T)

    }

    return(sig_list)
  })



  ##### ===== Salmon expression df

  txi <- reactive({
     req(folder_path())

      files <- list.files(paste0(folder_path(),"/salmon"), pattern = "quant.sf", recursive = T, full.names = T)

      if(length(files)==0 | is.null(files)){
        message("No expression file detected")
      } else {
        names(files) <- gsub(".salmon/quant.sf","", list.files(paste0(folder_path(),"/salmon"), pattern = "quant.sf", recursive = T))
        dat <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene, txOut = F)

        # files <- list.files(path1(), pattern = "quant.sf", recursive = T, full.names = T)
        # names(files) <- gsub(".salmon/quant.sf","", list.files(path1(), pattern = "quant.sf", recursive = T))
        # data <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene, txOut = F)

        return(dat)

      }
  })

  df_xp <- reactive({
    req(r$test$txi)
    dat <- log2(r$test$txi$abundance+1) %>%
      data.frame(check.names = F) %>%
      rownames_to_column("gene_id") %>% inner_join(gene_anno) %>%
      select(gene_name,  everything(),-gene_id)

    return(dat)

  })

  #### ===== DEV
  # df_xp <- reactive({
  #
  #   dat <- read.table("data/220107_Salmon_log2TPM_dev.tsv", sep = "\t", stringsAsFactors = F, header=T) %>% distinct(gene_name, .keep_all = T)
  #   return(dat)
  # })



  ##### ===== BED file / panel df

  bed <- reactive({

    if (input$panel == "None") {
      dat <- NULL
    } else {
      dat <- bed_list[[input$panel]]

    }
    return(dat)[,1:4]

  })


  # gene_panel <- reactive({
  #
  #   req(bed())
  #   dat <- bed() %>% distinct(gene_name, .keep_all = F) %>% unlist() %>% as.character()
  #   return(dat)
  #
  # })

  ##### ===== Filtered dataset

  df_filt <- reactive({

    req(df_xp())

    if (input$panel == "None") {
      dat <- df_xp() %>% distinct(gene_name, .keep_all = T)
    } else {
      dat <- df_xp() %>% filter(gene_name %in% unique(bed()$gene_name)) %>% distinct(gene_name, .keep_all = T)
    }
    return(dat)

  })

  ##### ===== Featurecount


  ftcount <- reactive({
    req(folder_path())

    files <- list.files(paste0(folder_path(),"/featurecount"), pattern = "featurecount.cnt$", recursive = T, full.names = T)

    if(length(files)==0 | is.null(files)){
      message("No ftcount expression file detected")
    } else {

    dat <- read.delim(paste0(folder_path(),"/featurecount/featurecount.cnt"), skip=1, stringsAsFactors = F)
    colnames(dat) <- gsub("X.*.bam.(.*.)_hg19_BQSR.bam$","\\1",colnames(dat))

    TPM <- tpm_approx( dat[,7:ncol(dat)], dat$Length)
    colnames(TPM) <- paste0(colnames(TPM),"_TPM")

    CPM <- cpm_approx( dat[,7:ncol(dat)])
    colnames(CPM) <- paste0(colnames(CPM),"_CPM")

    dat <- cbind(dat,TPM,CPM)

    return(dat)
    }

  })


  ##### ===== Fusion df

  ### === arriba

  df_fusion <- reactive({

    req(folder_path())
    files <- list.files(paste0(folder_path(),"/fusion/arriba"),pattern="_arriba_fusions.tsv", recursive = T, full.names = T)
    names(files) <- gsub("_arriba_fusions.tsv$","", list.files(paste0(folder_path(),"/fusion/arriba"), pattern = "_arriba_fusions.tsv", recursive = T))

    if(length(files)==0 | is.null(files)){
      message("No fusion file detected")
    } else {

    file_list <- list()
    for(i in 1:length(files)){

      file_list[[names(files[i])]] <- data.table::fread(files[i], data.table = F)

    }
    dat <- bind_rows(file_list, .id = "sample_id") %>%
      select(sample_id, gene1 =`#gene1`, gene2, type, reading_frame, breakpoint1, site1, breakpoint2,site2, everything()) %>%
      mutate(tool = "arriba", fusion = paste0(gene1,"-",gene2))

    return(dat)
    }

  })


  ##### ===== Variant df

  ### === Hotspot

  df_hotspot <- reactive({

    req(folder_path())
    files <- list.files(paste0(folder_path(),"/variant/hotspot"), full.names = T)
    names(files) <- gsub(".hotSpot.txt$","", list.files(paste0(folder_path(),"/variant/hotspot")))


    if(length(files)==0 | is.null(files)){
      message("No hotspot file detected")
    } else {

    file_list <- list()
    for(i in 1:length(files)){

      file_list[[names(files[i])]] <- read.delim(files[i], stringsAsFactors = F, row.names = NULL, na.strings = c("","."))

      colnames(file_list[[names(files[i])]]) <- c("chr","position","gene_name","mutation","predicted_impact","Ref",
                                                  "total_Depth","ref_fw","ref_rev",
                                                  "A_fw","A_rev","C_fw","C_rev","G_fw","G_rev","T_fw","T_rev",
                                                  "Ins","Del","Temp")

    }

    dat <- bind_rows(file_list, .id = "sample_id") %>% select(-Temp) %>%
      #mutate_at(8:20, function(x){as.numeric(x)}) %>%
      rowwise() %>%
      mutate(A_percent = sum(c(A_fw,A_rev),na.rm = T) / total_Depth,
             C_percent = sum(c(C_fw,C_rev),na.rm = T) / total_Depth,
             G_percent = sum(c(G_fw,G_rev),na.rm = T) / total_Depth,
             T_percent = sum(c(T_fw,T_rev),na.rm = T) / total_Depth,
             overall_count = sum(c(A_fw,A_rev,C_fw,C_rev,G_fw,G_rev,T_fw,T_rev),na.rm = T)) %>%
      mutate(overall_percent = overall_count/total_Depth) %>%
      ungroup() %>%
      mutate(gene_mut = paste0(gene_name,"-",mutation))

    return(dat)

}
  })



  ##### ===== Meta

  meta <- reactive({

    req(input$meta)

    if(input$sep == "xl"){
      data <- readxl::read_xlsx(input$meta$datapath) %>% as.data.frame()
      colnames(data)[1] <- "patient_id"

    } else {
      data <- data.table::fread(input$meta$datapath, sep = input$sep, dec = input$dec,
                                na.strings = c("", "NA", "#N/A"), stringsAsFactors = T, data.table = F)
      colnames(data)[1] <- "patient_id"

    }
    return(data)

  })


  observe({
    updateSelectInput(
      session,
      "panel",
      choices = c("None", names(bed_list))
    )
  })


  observe({
    r$test$folder_path <- folder_path()
  })

  observe({
    r$test$sig_list <- sig_list()
  })

  observe({
    r$test$txi <- txi()
    r$test$df_xp <- df_xp()
  })

  observe({
    r$test$bed <- bed()
  })


  observe({
   r$test$df_filt <- df_filt()
  })

  observe({
    r$test$ftcount <- ftcount()
  })

  observe({
    r$test$meta <- meta()
  })

  observe({
    r$test$df_fusion <- df_fusion()
  })

  observe({
    r$test$df_hotspot <- df_hotspot()
  })


##### ===== Check if files are correctly loaded

  output$salmon_check <- renderText({
    validate(
      need(r$test$txi, "Salmon: No expression data detected" )
           )

    paste("<font color=\"#21bf88\"><b>Salmon: Expression data successfully loaded!</b></font>"  )

      })

  output$salmon_dim <- renderText({

    paste("Total features in expression dataset:",nrow(txi()[[1]]), " / n samples: ",ncol(txi()[[1]]))

  })
  output$panel_dim <- renderText({

    paste("Panel unique features: ",length(unique(bed()$gene_name)), " (",nrow(df_filt())," found in exp. dataset)")

  })

  output$missing_genes <- renderText({

    unique(setdiff(as.character(unique(bed()$gene_name)), as.character(unique(df_filt()[,1]))))

  })

  output$fusion_check <- renderText({
    validate(
      need(r$test$df_fusion, "Arriba: No fusion data detected" )
    )

    paste("<font color=\"#21bf88\"><b>Arriba: Fusion data successfully loaded!</b></font>"  )

  })

  output$hotspot_check <- renderText({
    validate(
      need(r$test$df_hotspot, "Hotspot: No mutation data detected" )
    )

    paste("<font color=\"#21bf88\"><b>Hotspot: Mutation data successfully loaded!</b></font>"  )

  })

  output$cnv_check <- renderText({
    validate(
      need(r$test$ftcount, "CNV: No ftcount data detected" )
    )

    paste("<font color=\"#21bf88\"><b>CNV: Ftcount data successfully loaded!</b></font>"  )

  })




  output$preview_data <- DT::renderDT(
    ftcount(), # data
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
      paste("Salmon_log2pTPM.tsv")
    },
    content = function(file) {
      write.table(df_filt(), file, row.names = FALSE, sep = "\t", quote = F)
    }
  )


  output$preview_meta <- DT::renderDT(
    meta(), # data
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

  #################### ==================== App modules ====================  ####################

  mod_home_server("home_ui_1")
  mod_data_server("data_ui_1")
  mod_expression_server("expression_ui_1", r=r)
  mod_expression_deseq_server("expression_deseq_ui_1", r=r)
  mod_expression_individual_server("expression_individual_ui_1", r=r)
  mod_expression_signatures_server("expression_signatures_ui_1", r=r)

  mod_fusion_server("fusion_ui_1", r=r)

  mod_variant_server("variant_ui_1", r=r)
  mod_variant_GATK_server("variant_GATK_ui_1", r=r)

  mod_hotspot_server("hotspot_ui_1", r=r)

  mod_CNV_server("CNV_ui_1", r=r)


}
