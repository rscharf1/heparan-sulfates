library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)

healthy <- readRDS("inputs/Healthy.rds")

hs_enzymes <- c("Hs2st1", "Hs3st1", "Hs3st3b1", "Hs6st1", "Hs6st2", "Ndst3", 
"Ndst2", "Ndst1", "Glce", "Ext2", "Extl2", "Extl1", "Extl3", 
"Ext1", "Gpc4", "Tgfbr3", "Agrn")

genes <- healthy[["RNA"]]@counts %>% rownames()

length(genes)

hs_enzymes %in% genes