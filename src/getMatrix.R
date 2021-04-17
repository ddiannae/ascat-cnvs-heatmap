library(readr)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = T)

if (length(args) < 3 ) {
  stop("Incorrect number of arguments", call.=FALSE)
} else {
  tissue = args[1]
  type = args[2]
  datadir = args[3]
}
datadir <- paste(datadir, tissue, sep="/")
rawdir <- paste(datadir, "raw", paste0(tissue,"-", type, "-ascat"), sep="/")

files_to_read <- list.files(path = rawdir, 
                            pattern = "\\.tsv$", full.names = T, recursive = T)

all_files <- lapply(files_to_read, function(file) {
  cnvs <- read_tsv(file)
  cnvs <- cnvs %>% mutate(gene_id = as.character(lapply(str_split(gene_id, "\\."), "[[", 1)))
  return(cnvs %>% select(gene_id, copy_number))
})
  
size <- unique(do.call(rbind,lapply(all_files, dim)))
stopifnot(nrow(size)==1)

genes <- do.call(cbind, lapply(all_files, function(x) select(x, "gene_id")))
genes <- t(unique(t(genes)))
stopifnot(dim(genes)==c(size[1,1], 1))

targets <- data.frame(id = paste(tissue, type, 1:length(files_to_read), sep = "_"), 
                      file = unlist(lapply(str_split(files_to_read, "/"), "[[", 6)))


matrix <- bind_cols(lapply(all_files, function(x) select(x, "copy_number")))
colnames(matrix) <- targets$id
matrix <- matrix %>% mutate(ensembl_id = genes[,1]) %>% 
  select(ensembl_id, everything())

write_tsv(matrix, paste0(datadir, "/", tissue, "-", type, "-ascat-matrix.tsv"))
write_tsv(targets, paste0(datadir,"/", tissue, "-", type, "-ascat-files.tsv"))
