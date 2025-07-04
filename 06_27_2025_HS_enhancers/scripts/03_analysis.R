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

dt[grep("Lmo2", target_gene)][order(-active_score_mean)]

##############
# PCHi-C data
	# Bait is the target gene 
	# Other is the "enhancer"
dt <- fread("inputs/HPC7_Promoter_Capture_Interactions.ibed")

dt$enhancer_length <- dt$otherEnd_end - dt$otherEnd_start

dt[grep("Lmo2", bait_name)]


tmp <- dt[grep("Lmo2", bait_name)][, c(6,7), with = FALSE]
tmp[[2]] <- tmp[[2]] - min(tmp[[2]])
tmp[[1]] <- tmp[[1]] - min(tmp[[1]])

##############
# HS genes
hs <- fread("inputs/HS_relevant_genes.csv")

dt <- fread("inputs/HPC7_Promoter_Capture_Interactions.ibed")

dt <- dt[, c("bait_name", "otherEnd_chr", "otherEnd_start", "otherEnd_end", "N_reads", "score")]

dt[, gene_name := sub("-\\d+$", "", strsplit(bait_name, ",")[[1]][1]), by = .I]

dt[tolower(gene_name) %in% tolower(hs$hgnc_symbol)]

hs$hgnc_symbol[!(tolower(hs$hgnc_symbol) %in% tolower(dt$gene_name))]

hs$hgnc_symbol[tolower(hs$hgnc_symbol) %in% tolower(dt$gene_name)]

##############
# Plot Gardener
library(plotgardener)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)





