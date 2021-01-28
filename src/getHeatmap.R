library(vroom)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

figdir = "fig"

args <- commandArgs(trailingOnly = T)

if (length(args) < 3 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  tissue = args[1]
  biomart = args[2]
  matrix_file = args[3]
  
  if(!is.null(args[4])) {
    figdir = args[4]
  }
}

cnv_matrix <- vroom(matrix_file)
genes <- unlist(cnv_matrix[, 1], use.names = F)
cnv_matrix <- as.matrix(cnv_matrix[,-1])
rownames(cnv_matrix) <- genes

annot <- vroom(biomart, 
               col_names = c("id", "chr", "start", "end", "gc_content", "type"),
               skip = 1)
chrs <- c(as.character(1:22), "X")

annot <- annot %>% 
  filter(id %in% genes, chr %in% chrs) %>%  
  mutate(chr = factor(chr, levels = chrs)) %>%  
  arrange(chr, start) 
 
genes <- annot$id
cnv_matrix <- cnv_matrix[genes, ]

col_fun = colorRamp2(c(0, 2, 6), c("blue", "white", "red"))
chromosomes.pal <- c("#D909D1", "#0492EE", "#5DA0CB", "#106F35", "#5BD2AE", "#199F41", 
                     "#FE0F43", "#00FFCC", "#F495C5", "#E1BF5D", "#5F166F", "#088ACA",
                     "#41CFE0", "#0F0A71", "#FFFF99", "#B06645", "#800092", "#B925AE",
                      "#B1B719", "#CB97E8", "#130B9E", "#E12B29", "#79A5B9")

names(chromosomes.pal) <- c("22","11","12","13","14","15","16","17","18","19","1" ,"2" ,"3" ,"4" ,"5" ,
                            "6" ,"7" ,"X" ,"8" ,"9" ,"20","10","21")

chrs <- as.data.frame(annot %>% select(chr) %>% 
                        rename(Chr = chr))

ha <- HeatmapAnnotation(df = chrs, name = "Chr", show_annotation_name = F,
                        col = list(Chr = chromosomes.pal),
                        which = "row", width = unit(0.5, "cm"))

ht <- Heatmap(cnv_matrix, cluster_rows = F, cluster_columns = F, 
        col = col_fun, name = "CN", na_col = "#000000",
        show_column_names = F, show_row_names = F)

png(filename = paste0(figdir,"/", tissue, "-ascat-heatmap.png"), width = 1200, height = 600)
ht + ha
dev.off()
