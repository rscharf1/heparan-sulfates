library(data.table)
library(dplyr)
library(ggplot2)

dt <- fread("outputs/active_enhancers_with_genes.bed")

colnames(dt) <- c(
	"other_chr", 
	"other_start", 
	"other_end", 
	"target_gene", 
	"active_chr", 
	"active_start", 
	"active_end", 
	"active_score_max", 
	"active_score_mean"
)

dt[, gene_name := sub("-\\d+$", "", strsplit(target_gene, ",")[[1]][1]), by = .I]

