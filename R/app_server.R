#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function( input, output, session ) {
  # Your application server logic

  options(shiny.maxRequestSize = 40 * 1024^2)

  #################### ==================== Dataframe loading  ====================  ####################

  ##### ===== Set R or load full environment
  r <- reactiveValues(
    test = reactiveValues()
    )

  observe({
    req(input$env_loading)
    r$test <- readRDS(input$env_loading$datapath)$test
  })

  ##### ===== Full results folder

  ### Folder loading
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

  ### Direct environment loading


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

  tpm <- reactive({
    req(folder_path())
    message("===== Importing Salmon files...")
    files <- list.files(paste0(folder_path(),"/salmon"), pattern = "_TPM.tsv", recursive = T, full.names = T)

    if(length(files)==0 | is.null(files)){
      message("No Salmon files detected")
    } else {

      names(files) <- gsub("_TPM.tsv","", list.files(paste0(folder_path(),"/salmon"), pattern = "_TPM.tsv", recursive = T))

      file_list <- list()
      for(i in 1:length(files)){

        file_list[[names(files)[i]]] <- data.table::fread(files[i], data.table = F)
        colnames(  file_list[[names(files)[i]]])[2] <- names(files)[i]

      }

      dat <- Reduce(dplyr::full_join, file_list)

      return(dat)
      message("Salmon files successfully imported")

    }
  })


  df_xp <- reactive({
    req(tpm())
    dat <-tpm() %>%
      mutate_at(-1, function(x){log2(x+1)})%>%
      inner_join(gene_anno) %>%
      select(gene_name,  everything(),-gene_id)

    return(dat)

  })

  ##### ===== BED file / panel df

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


  bed_list <- reactive({

    bed_dir <- system.file("extdata","bed_panels", package = "hemRNA")
    bed_files <- list.files(bed_dir)

    bed_list <- list()

    for(bed in bed_files){

      name <- gsub(".bed$","",bed)

      bed_list[[name]] <- read.table(paste0(bed_dir,"/",bed), sep = "\t", stringsAsFactors = F, header = F)
      colnames(bed_list[[name]])[1:4] <- c("chr","start","end","gene_name")

    }

    return(bed_list)

  })

  bed <- reactive({


    if (input$panel == "None") {
      dat <- NULL
    } else {
      dat <- bed_list()[[input$panel]]

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



  ##### ===== Fusion df

  ### === arriba

  df_arriba <- reactive({

    req(folder_path())
    message("===== Importing Arriba files...")
    files <- list.files(paste0(folder_path(),"/fusion/arriba"),pattern="_arriba_fusions.tsv", recursive = T, full.names = T)
    names(files) <- gsub("_arriba_fusions.tsv$","", list.files(paste0(folder_path(),"/fusion/arriba"), pattern = "_arriba_fusions.tsv", recursive = T))

    if(length(files) == 0 | is.null(files)){
      message("No Arriba files detected")
    } else {
    file_list <- list()
    for(i in 1:length(files)){

      file_list[[names(files[i])]] <- data.table::fread(files[i], data.table = F)

    }
    dat <- bind_rows(file_list, .id = "sample_id") %>%
      select(sample_id, gene1 =`#gene1`, gene2, type, reading_frame, breakpoint1, site1, breakpoint2,site2, everything()) %>%
      mutate(tool = "arriba", fusion = paste0(gene1,"-",gene2))

    return(dat)
    message("Arriba files successfully imported")
    }

  })

  ### === FusionCatcher

  df_fusioncatcher <- reactive({

    req(folder_path())
    message("===== Importing NF-Core FusionCatcher files...")
    files <- list.files(paste0(folder_path(),"/nf-core/ArborescenceParEchantillon"), pattern = "_fusioncatcher.txt", recursive = T, full.names = T)

    names(files) <- gsub("/Fusioncatcher/.*.$","",
                         list.files(paste0(folder_path(),"/nf-core/ArborescenceParEchantillon"), pattern = "_fusioncatcher.txt", recursive = T))


    if(length(files)==0 | is.null(files)){
      message("No FusionCatcher files detected")
    } else {

      file_list <- list()
      for(i in 1:length(files)){

        file_list[[names(files[i])]] <- data.table::fread(files[i], data.table = F)

      }
      dat <- bind_rows(file_list, .id = "sample_id") %>%
        mutate(type = NA, site1 = NA, site2 = NA) %>%
        select(sample_id, gene1 =`Gene_1_symbol(5end_fusion_partner)`, gene2 = `Gene_2_symbol(3end_fusion_partner)`,
               type, reading_frame = Predicted_effect,
               breakpoint1 = `Fusion_point_for_gene_1(5end_fusion_partner)`, site1,
               breakpoint2 = `Fusion_point_for_gene_2(3end_fusion_partner)`, site2, everything()) %>%
        mutate(tool = "fusion_catcher", fusion = paste0(gene1,"-",gene2))

      return(dat)
      message("FusionCatcher files successfully imported")
    }

  })

  ### === Starfusion

  df_starfusion <- reactive({

    req(folder_path())
    message("===== Importing NF-Core StarFusion files...")
    files <- list.files(paste0(folder_path(),"/nf-core/ArborescenceParEchantillon"), pattern = "_star-fusion.tsv", recursive = T, full.names = T)

    names(files) <- gsub("/Star-Fusion/.*.$","",
                         list.files(paste0(folder_path(),"/nf-core/ArborescenceParEchantillon"), pattern = "_star-fusion.tsv", recursive = T))


    if(length(files)==0 | is.null(files)){
      message("No StarFusion files detected")
    } else {

      file_list <- list()
      for(i in 1:length(files)){

        file_list[[names(files[i])]] <- data.table::fread(files[i], data.table = F)

      }
      dat <- bind_rows(file_list, .id = "sample_id") %>%
        separate(`#FusionName`, sep = "--", into=c("gene1", "gene2")) %>%
        mutate(type = NA, site1 = NA, site2 = NA, reading_frame = NA) %>%
        select(sample_id, gene1, gene2,
               type, reading_frame,
               breakpoint1 = LeftBreakpoint, site1,
               breakpoint2 = RightBreakpoint, site2, everything()) %>%
        mutate(tool = "star_fusion", fusion = paste0(gene1,"-",gene2))

      return(dat)
      message("StarFusion files successfully imported")
    }

  })


  ##### ===== Variant df

  ### === Hotspot

  df_hotspot <- reactive({

    req(folder_path())
    message("===== Importing hotspot files...")
    files <- list.files(paste0(folder_path(),"/variant/hotspot"), full.names = T)
    names(files) <- gsub("_hg19_Aligned.sortedByCoord.out.hotSpot.txt$","", list.files(paste0(folder_path(),"/variant/hotspot")))


    if(length(files)==0 | is.null(files)){
      message("No hotspot files detected")
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
    message("Hotspot files successfully imported")

}
  })

  ### === RNAmut

  df_RNAmut <- reactive({

    req(folder_path())
    message("===== Importing RNAmut files...")
    files <- list.files(paste0(folder_path(),"/variant/RNAmut"),pattern="all_Sample.txt", full.names = T)
    names(files) <- gsub("_mutation-all_Sample.txt$","", list.files(paste0(folder_path(),"/variant/RNAmut"),pattern="all_Sample.txt"))


    if(length(files)==0 | is.null(files)){
      message("No RNAmut files detected")
    } else {

      file_list <- list()
      for(i in 1:length(files)){

        file_list[[names(files[i])]] <- read.delim(files[i], dec=",", stringsAsFactors = F, row.names = NULL) %>%
          mutate(WTReads = as.character(WTReads))

      }

      dat <- bind_rows(file_list, .id = "sample_id") %>%
        mutate(type = factor(if_else(grepl("-",Gene), "Fusion","Variant"),levels = c("Fusion","Variant")),
               tot_reads = sapply(strsplit(WTReads, split="-|_"),
                                  function(x){sum(unlist(as.numeric(x)))}),
               ProtMut_cut =  if_else(nchar(ProtMut) > 20,paste0(substr(ProtMut, 0, 20),"..."), ProtMut )
        ) %>%
        mutate(gene_mut = if_else(type == "Variant", paste0(Gene," [",ProtMut_cut,"]"), Gene),
               VAF2 = if_else(type == "Variant", VAF, MutReads / tot_reads))

      return(dat)

    }

  })

  ##### ===== CNV df

  ### === CNVkit

  cnvkit_list <- reactive({

    req(folder_path())
    message("===== Importing CNVkit files...")
    files <- list.files(paste0(folder_path(),"/cnvkit"), pattern = "_featurecount_genelevel_ENSID.cnr", recursive = T, full.names = T)

    if(length(files)==0 | is.null(files)){
      message("No cnvkit files detected")
    } else {

      names(files) <- gsub("_featurecount_genelevel_ENSID.cnr","",
                           list.files(paste0(folder_path(),"/cnvkit"), pattern = "_featurecount_genelevel_ENSID.cnr", recursive = T))

      file_list <- list()

      for(i in 1:length(files)){

        file_list[[names(files)[i]]] <- list(cnr = data.table::fread(files[i]) %>% mutate(order = seq(1:nrow(.))),
                                             base = data.table::fread(gsub(".cnr","_base.cns",files[i])),
                                             cbs = data.table::fread(gsub(".cnr","_cbs.cns",files[i])),
                                             hmm = data.table::fread(gsub(".cnr","_hmm.cns",files[i])))

        file_list[[names(files)[i]]]$cnr_by_cytob <- file_list[[names(files)[i]]]$cnr %>%
          mutate(chr = paste0("chr",chromosome)) %>%
          group_by(chr, gene) %>%
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

        file_list[[names(files)[i]]]$cnr_by_cytob <- file_list[[names(files)[i]]]$cnr %>%
          mutate(chr = paste0("chr",chromosome)) %>%
          group_by(chromosome) %>% mutate(vline = min(order)) %>% ungroup() %>%
          left_join(select(  file_list[[names(files)[i]]]$cnr_by_cytob, chr, gene, chr_cytob, chr_cytob2, chr_cytob3),
                    by = c("chr","gene")) %>%
          filter(is.na(chr_cytob3) == F)



      }
      return(file_list)

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
      choices = c("None", names(bed_list())),
      selected = names(bed_list())[1]
    )
  })


  observe({
    r$test$folder_path <- folder_path()
  })

  observe({
    r$test$sig_list <- sig_list()
  })

  observe({
    r$test$df_xp <- df_xp()
  })

  observe({
    r$test$bed <- bed()
  })


  observe({
   r$test$df_filt <- df_filt()
  })


  observe({
    r$test$meta <- meta()
  })

  observe({
    r$test$df_arriba <- df_arriba()
  })

  observe({
    r$test$df_fusioncatcher <- df_fusioncatcher()
  })

  observe({
    r$test$df_starfusion <- df_starfusion()
  })

  observe({
    r$test$df_hotspot <- df_hotspot()
  })

  observe({
    r$test$df_RNAmut <- df_RNAmut()
  })

  observe({
    r$test$cnvkit_list <- cnvkit_list()
  })


##### ===== Check if files are correctly loaded

  output$salmon_check <- renderText({
    validate(
      need(r$test$df_xp, "Salmon: No expression data detected" )
           )

    paste("<font color=\"#21bf88\"><b>Salmon: Expression data successfully loaded!</b></font>"  )

      })

  output$salmon_dim <- renderText({

    paste("Total features in expression dataset:",nrow(df_xp()[[1]]), " / n samples: ",ncol(df_xp()[[1]]))

  })
  output$panel_dim <- renderText({

    paste("Panel unique features: ",length(unique(bed()$gene_name)), " (",nrow(df_filt())," found in exp. dataset)")

  })

  output$missing_genes <- renderText({

    unique(setdiff(as.character(unique(bed()$gene_name)), as.character(unique(df_filt()[,1]))))

  })

  output$arriba_check <- renderText({
    validate(
      need(r$test$df_arriba, "Arriba: No fusion data detected" )
    )

    paste("<font color=\"#21bf88\"><b>Arriba: Fusion data successfully loaded!</b></font>"  )

  })

  output$fusion_catcher_check <- renderText({
    validate(
      need(r$test$df_fusioncatcher, "NF-Core: No fusion data detected" )
    )

    paste("<font color=\"#21bf88\"><b>NF-Core: Fusion catcher data successfully loaded!</b></font>"  )

  })

  output$hotspot_check <- renderText({
    validate(
      need(r$test$df_hotspot, "Hotspot: No mutation data detected" )
    )

    paste("<font color=\"#21bf88\"><b>Hotspot: Mutation data successfully loaded!</b></font>"  )

  })

  output$RNAmut_check <- renderText({
    validate(
      need(r$test$df_RNAmut, "RNAmut: No data detected" )
    )

    paste("<font color=\"#21bf88\"><b>RNAmut: Data successfully loaded!</b></font>"  )

  })

  output$cnv_check <- renderText({
    validate(
      need(r$test$cnvkit_list, "CNV: No CNVkit data detected" )
    )

    paste("<font color=\"#21bf88\"><b>CNV: CNVkit data successfully loaded!</b></font>"  )

  })


  ##### ===== Save environment


  output$save_env <- downloadHandler(
    filename = function() {
      paste("full_environment.rds")
    },
    content = function(file) {
      saveRDS(r, file)
    }
  )


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
  #mod_data_server("data_ui_1")

  mod_overview_server("overview_1", r=r)

  mod_expression_server("expression_ui_1", r=r)
  mod_expression_deseq_server("expression_deseq_ui_1", r=r)
  mod_expression_individual_server("expression_individual_ui_1", r=r)
  mod_expression_signatures_server("expression_signatures_ui_1", r=r)

  mod_fusion_server("fusion_ui_1", r=r)
  mod_fusion_catcher_server("fusion_catcher_ui_1", r=r)
  mod_star_fusion_server("star_fusion_ui_1", r=r)

  #mod_variant_server("variant_ui_1", r=r)
  #mod_variant_GATK_server("variant_GATK_ui_1", r=r)
  mod_hotspot_server("hotspot_ui_1", r=r)
  mod_RNAmut_server("RNAmut_ui_1", r=r)

  mod_CNV_server("CNV_ui_1", r=r)


}
