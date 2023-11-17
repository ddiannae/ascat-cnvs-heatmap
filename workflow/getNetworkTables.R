################################################################################
## Script to get interactions and vertices tables to build networks of a 
## top MI cutoff. The annotations are added from a biomart file 
## It requires the mi matrix as input.
###############################################################################
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vroom)
library(tidyr)
library(dplyr)
library(janitor)

CUTOFF <- as.numeric(snakemake@params[["cutoff"]])
COND <- snakemake@params[["cond"]]

cat("Loading files\n")
MImatrix <- vroom::vroom(snakemake@input[["mi_matrix"]], col_types = c(.default = "d"))

if(nrow(MImatrix) == 0) {
  vroom_write(tibble(),snakemake@output[["vertices"]])
  vroom_write(tibble(),snakemake@output[["interactions"]])
  quit()
}

genes <- colnames(MImatrix)

cat("MI matrix with ", nrow(MImatrix), " rows and ", ncol(MImatrix), " columns loaded \n")
MImatrix <- MImatrix %>% 
  as.matrix()

MImatrix[is.na(MImatrix)] <- 0
MImatrix[lower.tri(MImatrix, diag = T)] <- NA

MImatrix <- as_tibble(MImatrix)
MImatrix$source <- genes

cat("Annotating interactions\n")
MIvals <- MImatrix %>% pivot_longer(cols = starts_with("ENSG"), 
                                     names_to = "target",
                                     values_to = "mi",
                                     values_drop_na = TRUE) %>% 
  filter(source != target) %>%
  arrange(desc(mi)) %>% 
  mutate(row_num = row_number()) %>% 
  filter(row_num <= CUTOFF)
MIvals$cond <- COND

annot <-  vroom::vroom(snakemake@params[["biomart"]],  
                       col_names = c("ensembl_id", "chr", "start", "end", "gc", "type", "gene_name"), 
                         skip = 1) %>%
  select(ensembl_id, chr, start, end, gene_name)

cat("Merging annotations\n")
colnames(annot) <-  c("source", "source_chr", "source_start", "source_end","source_name")
MIvals <- merge(annot, MIvals)
colnames(annot) <-  c("target", "target_chr", "target_start", "target_end",  "target_name")
MIvals <- merge(annot, MIvals)
MIvals <- MIvals %>% mutate(inter = if_else(source_chr == target_chr,  F, T), 
                      interaction_type = if_else(inter == T, "Inter", "Intra"),
                      distance = if_else(inter == F, as.integer(pmax(source_start, target_start) - 
                                           pmin(source_start, target_start)), as.integer(-1)))

targets <- MIvals %>% select(target, target_chr, target_start, target_end, target_name)
sources <- MIvals %>% select(source, source_chr, source_start, source_end, source_name)
colnames(targets) <- c("ensembl", "chr", "start", "end",  "symbol")
colnames(sources)  <-  c("ensembl", "chr", "start", "end", "symbol")
vertices <- bind_rows(targets, sources)
vertices <- vertices[!duplicated(vertices$ensembl), ] 

MIvals <- MIvals %>% 
  select(source, target, mi, distance, row_num, interaction_type, cond) %>%
  arrange(row_num) 

cat("Saving files\n")
vroom_write(vertices, file = snakemake@output[["vertices"]])
vroom_write(MIvals, file = snakemake@output[["interactions"]])
