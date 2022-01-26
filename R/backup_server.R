
### === Fusion list
fusion_list <- reactive({
  req(folder_path())

fusion_list <- list()
### === arriba
files <- list.files(paste0(folder_path(),"/fusion/arriba"),pattern="_arriba_fusions.tsv", recursive = T, full.names = T)
names(files) <- gsub("_arriba_fusions.tsv$","", list.files(paste0(folder_path(),"/fusion/arriba"), pattern = "_arriba_fusions.tsv", recursive = T))

file_list <- list()
for(i in 1:length(files)){

  file_list[[names(files[i])]] <- data.table::fread(files[i], data.table = F)

}
fusion_list[["arriba"]] <- bind_rows(file_list, .id = "sample_id") %>%
  filter(confidence == "high") %>%
  select(sample_id, gene1 =`#gene1`, gene2, type, reading_frame, breakpoint1, site1, breakpoint2,site2) %>%
  mutate(tool = "arriba", fusion = paste0(gene1,"-",gene2))

### === Fusioncatcher
files <- list.files(paste0(folder_path(),"/fusion/fusioncatcher"),pattern="final-list_candidate-fusion-genes.txt", recursive = T, full.names = T)
names(files) <- gsub("/final-list_candidate-fusion-genes.txt$","", list.files(paste0(folder_path(),"/fusion/fusioncatcher"),pattern="final-list_candidate-fusion-genes.txt",
                                                                              recursive = T))

file_list <- list()
for(i in 1:length(files)){

  file_list[[names(files[i])]] <- data.table::fread(files[i])

}
fusion_list[["fusioncatcher"]] <- bind_rows(file_list, .id = "sample_id") %>%
  filter(!Fusion_finding_method %in% c("BOWTIE","STAR","BLAT")) %>%
  mutate(type = NA, site1 = NA, site2 = NA) %>%
  select(sample_id, gene1 =`Gene_1_symbol(5end_fusion_partner)`, gene2=`Gene_2_symbol(3end_fusion_partner)`, type,
         reading_frame = Predicted_effect,
         breakpoint1 = `Fusion_point_for_gene_1(5end_fusion_partner)`, site1,
         breakpoint2 = `Fusion_point_for_gene_2(3end_fusion_partner)` ,site2) %>%
  mutate(tool = "fusioncatcher", fusion = paste0(gene1,"-",gene2))

### === STARfusion
files <- list.files(paste0(folder_path(),"/fusion/starfusion"),pattern="_star-fusion.fusion_predictions.tsv", recursive = T, full.names = T)
names(files) <- gsub("_star-fusion.fusion_predictions.tsv$","", list.files(paste0(folder_path(),"/fusion/starfusion"),pattern="_star-fusion.fusion_predictions.tsv",
                                                                           recursive = T))

file_list <- list()
for(i in 1:length(files)){

  file_list[[names(files[i])]] <- data.table::fread(files[i])

}
fusion_list[["starfusion"]] <- bind_rows(file_list, .id = "sample_id")

fusion_list[["starfusion"]] <- fusion_list[["starfusion"]] %>%
  mutate(gene1 = lapply(strsplit(fusion_list[["starfusion"]]$`#FusionName`, "--"), function(x){x[1]}),
         gene2 = lapply(strsplit(fusion_list[["starfusion"]]$`#FusionName`, "--"), function(x){x[2]}),
         type = NA,reading_frame=NA, site1 = NA, site2 = NA) %>%
  select(sample_id, gene1, gene2, type,
         reading_frame,
         breakpoint1 = LeftBreakpoint, site1,
         breakpoint2 = RightBreakpoint ,site2) %>%
  mutate(tool = "starfusion", fusion = paste0(gene1,"-",gene2))

return(fusion_list)

})

df_fusion <- reactive({

  req(fusion_list())
  dat <- Reduce(rbind,fusion_list())

  return(dat)

})
