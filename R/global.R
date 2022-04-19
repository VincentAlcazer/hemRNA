
#' Global parameters

#  library(BiocManager)
# options(repos = BiocManager::repositories())

#' @import shinyFiles
#' @import ggplot2
#' @import dplyr
#' @import tibble
#' @import DESeq2
#' @import ggrepel
#' @import viridis
#' @import tidyr
#' @importFrom data.table fread
#' @import BioCircos
#' @import rmarkdown
#' @import ComplexHeatmap
#' @import grid
#' @import circlize
#'

#library(ggplot2)

########## ========== Datasets

tx2gene <- read.delim("data/gencode_v38_tx2Gene.tsv")
gene_anno <- read.delim("data/gencode_v38_gene_names.tsv")

# sig_files <- list.files("inst/extdata/signatures/")
#
# sig_list <- list()
# for(sig in sig_files){
#
#     name <- gsub(".txt$","",sig)
#
#     sig_list[[name]] <- read.table(paste0("inst/extdata/signatures/",sig), sep = "\t", stringsAsFactors = F, header = T)
#
# }

# bed_files <- list.files("inst/extdata/bed_panels/")
#
# bed_list <- list()
# for(bed in bed_files){
#
#   name <- gsub(".bed$","",bed)
#
#   bed_list[[name]] <- read.table(paste0("inst/extdata/bed_panels/",bed), sep = "\t", stringsAsFactors = F, header = F)
#   colnames(bed_list[[name]])[1:4] <- c("chr","start","end","gene_name")
#
# }


hg19_cytoband <- read.table("data/UCSC_hg19_cytoBand.txt.gz",
                            sep = "\t", stringsAsFactors = F,
                            col.names = c("chr","start","end","cytoband","gstat"))

########## ========== Parameters

### Ggplot 2 default theme

default_theme <- theme_bw() + theme(
  plot.title = element_text(size = 18, face = "bold"),
  axis.text = element_text(size = 12, color = "black"),
  axis.title = element_text(size = 14, face = "bold"),
  legend.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12),
  legend.position = "none"
)
