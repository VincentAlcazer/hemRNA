

tpm_approx <- function(data,length){

  RPK = data / (length/1000)
  TPM = RPK/(colSums(RPK)/1e6)
  return(TPM)
}

cpm_approx <- function(data){

  CPM = (data/colSums(data))* 1e6
  return(CPM)
}

depth_approx <- function(data){

  CPM = data / colSums(data)
  return(CPM)
}


#' Draw final circle plot
#'

circos_plot <- function(expression,cnv,variant, arriba, bed,
                        SNP_max_rad = 1.145, SNP_min_rad =1,
                        exp_max_rad = 0.97, exp_min_rad = 0.80,
                        CNV_max_rad = 0.75, CNV_min_rad = 0.4,
                        fusion_max_rad = 0.35,
                        CNV_cutoff = 1.5
                        ){


  ########## ========== Global parameters



  ### Set genes coords
  gene_order <- gsub("chr","",unique(bed$chr))

  genes_coords <- bed %>%
    #mutate(chr = factor(chr, levels = unique(bed$chr))) %>%
    distinct(gene_name, .keep_all = T) %>%
    group_by(chr) %>%
    mutate(order = seq(1:n()),
           chr = factor(gsub("chr","",chr), levels = gene_order))

  ### Set custom genome coords
  custom_genome <- genes_coords %>%
    group_by(chr) %>%
    summarise(length = length(unique(gene_name))) %>%
    column_to_rownames("chr") %>% t() %>% as.data.frame(check.names = F) %>% as.list()

  ########## ========== Prepare data

  ##### ===== Expression

  salmon_res <- expression %>%
    #select(gene_name, expression = all_of(input$sample)) %>% # DEV TO REMOVE
    inner_join(genes_coords, by = "gene_name") %>%
    arrange(chr)

  ##### ===== CNV

  cnv_res <- cnv %>% left_join(genes_coords, by = c("gene"="gene_name"),
                                 suffix = c("_cnv","_coords"))


  cnv_res_notsig <- cnv_res %>% filter(abs(CNV) < CNV_cutoff)
  cnv_res_sig <- cnv_res %>% filter(abs(CNV) >= CNV_cutoff)


  ##### ===== Mutations

  RNAmut_res <- variant %>% filter(Info == "Oncogenic", type == "Variant") %>%
    #filter(sample_id == input$sample) %>% # DEV TO REMOVE
    arrange(desc(VAF2)) %>%
    distinct(gene_mut, .keep_all = T) %>%
    left_join(genes_coords, by = c("Gene"="gene_name"))

  ##### ===== Fusions

  fusion_res <- arriba %>% select(-read_identifiers, -fusion_transcript) %>%
    #filter(sample_id == input$sample) %>% # DEV TO REMOVE
    filter(confidence == "high") %>%
    left_join(genes_coords, by = c("gene1"="gene_name")) %>%
    left_join(genes_coords, by = c("gene2"="gene_name"), suffix = c("_g1","_g2")) %>%
    na.omit()
    # separate(breakpoint1, into = c("chr1","linkpos1")) %>%
    # separate(breakpoint2, into = c("chr2","linkpos2"))



  ########## ========== Draw tracks


  ##### ===== SNP tracks

  if(nrow(RNAmut_res) == 0){

    tracks = BioCircosTextTrack("text1",
                                         'No variants found',
                                         weight = "normal",
                                         x = -0.2, y = -1.07)
    tracks = tracks + BioCircosTextTrack("text2",
                                         'by RNAMut',
                                         weight = "normal",
                                         x = -0.1, y = -1.02)

  } else {

    tracks =  BioCircosArcTrack("RNAmut_variants",
                                        chromosome = as.character(RNAmut_res$chr),
                                        starts = RNAmut_res$order,
                                        ends = RNAmut_res$order+4,
                                        labels = RNAmut_res$gene_mut,

                                        maxRadius = SNP_max_rad, minRadius = SNP_min_rad,
                                        color = "black")

  }


  ##### ===== Expression track
  tracks = tracks + BioCircosHeatmapTrack("heatmap1",
                                          as.character(salmon_res$chr),
                                          salmon_res$order-1,
                                          salmon_res$order,
                                          round(salmon_res$expression,2),
                                          labels = salmon_res$gene_name,
                                          minRadius = exp_min_rad,
                                          maxRadius = exp_max_rad)

  ##### ===== CNV tracks

  tracks = tracks + BioCircosCNVTrack("cnv_backm2",
                                      as.character(cnv_res$chromosome),
                                      cnv_res$order_coords-1,
                                      cnv_res$order_coords,
                                      rep(-2,length(as.character(cnv_res$chromosome))),
                                      color = "black",
                                      range = c(-2,2),
                                      minRadius = CNV_min_rad,
                                      maxRadius = CNV_max_rad)

  tracks = tracks + BioCircosCNVTrack("cnv_backm1",
                                      as.character(cnv_res$chromosome),
                                      cnv_res$order_coords-1,
                                      cnv_res$order_coords,
                                      rep(-1,length(as.character(cnv_res$chromosome))),
                                      color = "lightgrey",
                                      range = c(-2,2),
                                      minRadius = CNV_min_rad,
                                      maxRadius = CNV_max_rad)

  tracks = tracks + BioCircosCNVTrack("cnv_back0",
                                      as.character(cnv_res$chromosome),
                                      cnv_res$order_coords-1,
                                      cnv_res$order_coords,
                                      rep(0,length(as.character(cnv_res$chromosome))),
                                      color = "grey",
                                      range = c(-2,2),
                                      minRadius = CNV_min_rad,
                                      maxRadius = CNV_max_rad)

  tracks = tracks + BioCircosCNVTrack("cnv_back1",
                                      as.character(cnv_res$chromosome),
                                      cnv_res$order_coords-1,
                                      cnv_res$order_coords,
                                      rep(1,length(as.character(cnv_res$chromosome))),
                                      color = "lightgrey",
                                      range = c(-2,2),
                                      minRadius = CNV_min_rad,
                                      maxRadius = CNV_max_rad)

  tracks = tracks + BioCircosCNVTrack("cnv_back2",
                                      as.character(cnv_res$chromosome),
                                      cnv_res$order_coords-1,
                                      cnv_res$order_coords,
                                      rep(2,length(as.character(cnv_res$chromosome))),
                                      color = "black",
                                      range = c(-2,2),
                                      minRadius = CNV_min_rad,
                                      maxRadius = CNV_max_rad)

  tracks = tracks + BioCircosCNVTrack("cnv_notsig",width = 2,
                             as.character(cnv_res_notsig$chromosome),
                             cnv_res_notsig$order_coords-1,
                             cnv_res_notsig$order_coords,
                             cnv_res_notsig$CNV,
                             color = "cornflowerblue",
                             range = c(-2,2),
                             minRadius = CNV_min_rad,
                             maxRadius = CNV_max_rad)

  tracks = tracks + BioCircosCNVTrack("cnv_sig", width = 5,
                                      as.character(cnv_res_sig$chromosome),
                                      cnv_res_sig$order_coords-1,
                                      cnv_res_sig$order_coords,
                                      cnv_res_sig$CNV,
                                      color = "firebrick",
                                      range = c(-2,2),
                                      minRadius = CNV_min_rad,
                                      maxRadius = CNV_max_rad)





  ##### ===== Fusions track

  if(nrow(fusion_res) == 0){

    tracks = tracks + BioCircosTextTrack("text3",
                                        'No fusions found',
                                        weight = "normal",
                                        x = -0.2, y = 0)
    tracks = tracks + BioCircosTextTrack("text4",
                                         'by Arriba',
                                         weight = "normal",
                                         x = -0.10, y = 0.1)

  } else {

    tracks = tracks + BioCircosLinkTrack("fusions",
                                         gene1Chromosomes = as.character(fusion_res$chr_g1),
                                         gene1Starts = as.numeric(fusion_res$order_g1),
                                         gene1Ends = as.numeric(fusion_res$order_g1)+4,
                                         gene2Chromosomes = as.character(fusion_res$chr_g2),
                                         gene2Starts = as.numeric(fusion_res$order_g2),
                                         gene2Ends = as.numeric(fusion_res$order_g2)+4,
                                         labels = fusion_res$fusion,
                                         displayAxis = F, axisPadding = 6,
                                         color = "firebrick", width = "0.5em",
                                         displayLabel = F,
                                         maxRadius = fusion_max_rad)
  }

  ##### ===== Backgrounds
#
  # tracks = tracks + BioCircosBackgroundTrack("BG_SNP",
  #                                            borderColors = "black",
  #                                            fillColor = "white", borderSize = 0.6,
  #                                            minRadius = CNV_max_rad/2, maxRadius = CNV_max_rad)
  # tracks = tracks + BioCircosBackgroundTrack("BG_SNP_2",
  #                                            borderColors = "black",
  #                                            fillColor = "white", borderSize = 0.6,
  #                                            minRadius = CNV_min_rad, maxRadius = CNV_max_rad/2)
  #


  ##### ===== Final viz
  p <- BioCircos(tracks,
                 genome = custom_genome,
                 genomeFillColor = "Spectral",
                 yChr = F,
                 chrPad = 0,
                 displayGenomeBorder = F,
            genomeTicksLen = 3,
            genomeTicksTextSize = 0,
            genomeTicksScale = 50000000,
            genomeLabelTextSize = 18,
            genomeLabelDy = 0)

  return(p)



}




#' Input = df with samples as rownames and features as cols
#' Annotation = vector of annotation (in the same order as input's rownames!!)
#' anno_pal = viridis_pal, hue_pal, grey_pal...
#'


fast_map <- function(df, top_features_p = 1, center_scale = T,
                     heat_colors="ATAC",
                     anno_pal = viridis_pal(),
                     Zlim = 3, method = "ward.D2", metric = "euclidean", title = "title",
                     Data_input = "Z-score expression",
                     anno_1 = NULL, anno_1_name = NULL, anno_1_cols = NULL,
                     anno_2 = NULL, anno_2_name = NULL, anno_2_cols = NULL,
                     cluster_rows = T, cluster_cols = T,
                     reorder_rows = T, reorder_cols = T,
                     row_names = T,col_names=T,
                     row_fontsize = 12, col_fontsize = 12,
                     show_hclust_n = 0, show_pearson = F ){

  if(heat_colors == "Hiroshige"){
    heat_colors=rev(colorRampPalette(MetBrewer::MetPalettes$Hiroshige[[1]])(100))
  } else if(heat_colors == "ATAC"){
    heat_colors=colorRampPalette(c('#3361A5','#248AF3','#14B3FF','#88CEEF','#C1D5DC',
                                   '#EAD397','#FDB31A','#E42A2A','#A31D1D'))(100)
  }

  if(show_pearson == T){

    Data_input = "Pearson's R"
    heat_col <- circlize::colorRamp2(seq(-1, 1, by =2/99),heat_colors)

  } else {

    heat_col <- circlize::colorRamp2(seq(-Zlim, +Zlim, by =(2*Zlim)/99),heat_colors)

  }

  sup_leg = list(legend_1 = Legend(col_fun = heat_col,
                                   title = Data_input,
                                   labels_gp = gpar(fontsize = 12),
                                   direction = c("horizontal"),
                                   title_gp = gpar(fontsize = 12, fontface="bold"))
  )



  if(is.null(anno_1)){

    sup_ha = NULL

  } else if(is.null(anno_2)){

    if(is.numeric(anno_1) | is.double(anno_1)){

      ##### ===== 1 unique continuous annotation

      col_list <- list(supha_1 = circlize::colorRamp2(c(min(anno_1, na.rm=T),
                                                        median(anno_1,na.rm=T),
                                                        max(anno_1,na.rm=T)),
                                                      c(anno_pal(100)[1],
                                                        anno_pal(100)[50],
                                                        anno_pal(100)[100])))

      sup_ha =  HeatmapAnnotation(supha_1 = anno_1,
                                  show_annotation_name = F,
                                  col = col_list,
                                  show_legend = F,
                                  na_col = "grey")

      sup_leg$legend_2 = Legend(col_fun = col_list$supha_1,
                                title=anno_1_name,
                                direction = c("horizontal"), ncol=2,
                                labels_gp = gpar(fontsize = 12),
                                title_gp = gpar(fontsize = 12, fontface="bold"))


    } else {

      ##### ===== 1 unique categorical annotation

      if(is.null(anno_1_cols)){
        ha_1 <- anno_pal(length(unique(na.omit(anno_1))))
      } else {
        ha_1 <- anno_1_cols
      }

      names(ha_1) <- unique(na.omit(anno_1))

      col_list <- list(supha_1 = ha_1)

      sup_ha =  HeatmapAnnotation(supha_1 = anno_1,
                                  show_annotation_name = F,
                                  col = col_list,
                                  show_legend = F)

      sup_leg$legend_2 = Legend(labels=names(col_list[[1]]),
                                title=anno_1_name,
                                direction = c("horizontal"), ncol=2,
                                labels_gp = gpar(fontsize = 12),
                                title_gp = gpar(fontsize = 12, fontface="bold"),
                                legend_gp = gpar(fill = as.vector(col_list[[1]])))

    }



  } else {

    ##### ===== 2 annotations

    ### CHeck if anno 1 is num
    if(is.numeric(anno_1) | is.double(anno_1)){

      ha_1 = circlize::colorRamp2(c(min(anno_1, na.rm=T),
                                    median(anno_1,na.rm=T),
                                    max(anno_1,na.rm=T)),
                                  c(anno_pal(100)[1],
                                    anno_pal(100)[50],
                                    anno_pal(100)[100]))

      sup_leg$legend_2 = Legend(col_fun = ha_1,
                                title=anno_1_name,
                                direction = c("horizontal"), ncol=2,
                                labels_gp = gpar(fontsize = 12),
                                title_gp = gpar(fontsize = 12, fontface="bold"))

    } else {

      if(is.null(anno_1_cols)){
        ha_1 <- anno_pal(length(unique(na.omit(anno_1))))
      } else {
        ha_1 <- anno_1_cols
      }


      names(ha_1) <- unique(na.omit(anno_1))

      sup_leg$legend_2 = Legend(labels=names(ha_1),
                                title=anno_1_name,
                                direction = c("horizontal"), ncol=2,
                                labels_gp = gpar(fontsize = 12),
                                title_gp = gpar(fontsize = 12, fontface="bold"),
                                legend_gp = gpar(fill = as.vector(ha_1)))


    }

    ### CHeck if anno 2 is num

    if(is.numeric(anno_2) | is.double(anno_2)){

      ha_2 = circlize::colorRamp2(c(min(anno_2, na.rm=T),
                                    median(anno_2,na.rm=T),
                                    max(anno_2,na.rm=T)),
                                  c(anno_pal(100)[1],
                                    anno_pal(100)[50],
                                    anno_pal(100)[100]))


      sup_leg$legend_3 = Legend(col_fun = ha_2,
                                title=anno_2_name,
                                direction = c("horizontal"), ncol=2,
                                labels_gp = gpar(fontsize = 12),
                                title_gp = gpar(fontsize = 12, fontface="bold"))



    } else {

      if(is.null(anno_2_cols)){
        ha_2 <- RColorBrewer::brewer.pal(length(unique(na.omit(anno_2))), "Dark2")[1:length(unique(na.omit(anno_2)))]
      } else {
        ha_2 <- anno_2_cols
      }

      names(ha_2) <- unique(na.omit(anno_2))

      sup_leg$legend_3 = Legend(labels=names(ha_2),
                                title=anno_2_name,
                                direction = c("horizontal"), ncol=2,
                                labels_gp = gpar(fontsize = 12),
                                title_gp = gpar(fontsize = 12, fontface="bold"),
                                legend_gp = gpar(fill = as.vector(ha_2)))



    }

    sup_ha =  HeatmapAnnotation(supha_1 = anno_1,
                                supha_2 = anno_2,
                                show_annotation_name = F,
                                col = list(supha_1 = ha_1,
                                           supha_2 = ha_2),
                                show_legend = F)






  }


  ##### ===== Prepare df

  matrix <- df
  matrix <- matrix[,order(apply(matrix, 2, var), decreasing = T)[1:round(top_features_p*ncol(matrix))]]

  if(center_scale == T){
    matrix <- scale(matrix)
  }

  # if(show_pearson == T){
  #
  #   matrix <- Hmisc::rcorr(t(matrix), type = "pearson")$r
  #
  #   # sup_leg$legend_1 <- Legend(col_fun = heat_col,
  #   #                            title = "Pearson's R",
  #   #                            labels_gp = gpar(fontsize = 12),
  #   #                            direction = c("horizontal"),
  #   #                            title_gp = gpar(fontsize = 12, fontface="bold"))
  #
  #   message("Max R: ",round(max(matrix[matrix < 1]),2),
  #           " - Min R: ", round(min(matrix[matrix > -1]),2))
  #
  #
  #
  # }


  ##### ===== Show Hclust n

  if(show_hclust_n > 0){

    model <- factoextra::eclust(matrix, "hclust", k = show_hclust_n, hc_metric = metric,
                                stand = FALSE, hc_method = method, graph = FALSE)

    if(show_hclust_n == 2){
      clusters <- RColorBrewer::brewer.pal(3, "Set1")[c(1,3)]
    } else {
      clusters <- RColorBrewer::brewer.pal(length(unique(model$cluster)), "Set1")
    }

    names(clusters) <- unique(model$cluster)

    bot_ha =  HeatmapAnnotation(botha1 = model$cluster,
                                show_annotation_name = F,
                                col = list(botha1 = clusters),
                                show_legend = F,
                                gp = gpar(fontsize = 12))

    bot_leg = Legend(
      labels =  names(clusters),
      title = "Clusters",
      direction = c("vertical"), ncol = 2,
      labels_gp = gpar(fontsize = 12),
      title_gp = gpar(fontsize = 12, fontface="bold"),
      legend_gp = gpar(fill = as.vector(clusters))
    )


  } else {

    bot_ha = NULL
    bot_leg = NULL
  }


  ht <-   ComplexHeatmap::Heatmap(t(matrix),
                  name = "heat",
                  col = heat_col,
                  show_column_names = col_names,
                  show_row_names = row_names,
                  #column_split = k,
                  cluster_columns = cluster_cols,
                  show_column_dend = cluster_cols,
                  column_dend_reorder =   reorder_cols,
                  clustering_distance_columns = metric,
                  clustering_method_column = method,

                  show_heatmap_legend = F,

                  show_row_dend = cluster_rows,
                  cluster_rows = cluster_rows,
                  row_dend_reorder = reorder_rows,
                  clustering_distance_rows = metric,
                  clustering_method_rows = method,
                  use_raster = TRUE,
                  raster_device = c("png"),
                  raster_quality = 2,
                  top_annotation = sup_ha,
                  bottom_annotation = bot_ha,
                  column_title = title,
                  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                  row_title_gp = gpar(fontsize = 12),

                  row_names_gp = gpar(fontsize = row_fontsize),
                  column_names_gp = gpar(fontsize = col_fontsize)
                  )


  draw(ht,annotation_legend_list = c(sup_leg,bot_leg), annotation_legend_side = "bottom", newpage = T)



}

