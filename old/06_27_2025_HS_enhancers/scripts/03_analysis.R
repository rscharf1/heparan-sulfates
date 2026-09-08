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
library(dplyr)
library(ggplot2)

pchic <- read.table("inputs/HPC7_Promoter_Capture_Interactions.ibed", header = TRUE, sep = "\t")
# hk27ac <- read.table("inputs/GSM1329815_ChIPseq_H3K27Ac_HPC7.cleaned.bedgraph", header = FALSE)

gata2_bait <- pchic[grep("Gata2", pchic$bait_name), ]

gata2_bedpe <- gata2_bait[, c("bait_chr", "bait_start", "bait_end",
                              "otherEnd_chr", "otherEnd_start", "otherEnd_end")]

gata2_bedpe$bait_chr <- paste0("chr", gata2_bedpe$bait_chr)
gata2_bedpe$otherEnd_chr <- paste0("chr", gata2_bedpe$otherEnd_chr)

write.table(gata2_bedpe, "gata2_interactions.bedpe", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

pageCreate(width = 7, height = 3, default.units = "inches")

chrom <- "chr6"
region_start <- 88120000
region_end   <- 88480000

params <- pgParams(
  chrom = chrom,
  chromstart = region_start,
  chromend = region_end,
  assembly = "mm9",
  x = 0.5,
  y = 0.5,
  width = 6,
  height = 1,
  just = c("left", "top"),
  default.units = "inches"
)

pchic_plot <- plotPairs(
  data = "gata2_interactions.bedpe",
  params = params,
  fill = "darkred",
  linecolor = "darkred",
  lwd = 1.5,
  opacity = 0.8
)

annoGenomeLabel(
  plot = pchic_plot,
  params = params,
  scale = "Kb",
  fontsize = 8,
  y = 1.55
)

ggsave("gata2_pchic_interactions.pdf")

#################################
library(plotgardenerData)
data("IMR90_DNAloops_pairs")

pdf("test.pdf")

## Set the coordinates
params <- pgParams(
    chrom = "chr21",
    chromstart = 27900000, chromend = 30700000,
    assembly = "hg19",
    width = 7
)

## Create a page
pageCreate(width = 7.5, height = 2.1, default.units = "inches")

IMR90_DNAloops_pairs$length <- 
        (IMR90_DNAloops_pairs$start2 - IMR90_DNAloops_pairs$start1) / 1000

bedpePlot <- plotPairs(
    data = IMR90_DNAloops_pairs, params = params,
    fill = colorby("length", palette = 
                colorRampPalette(c("dodgerblue2", "firebrick2"))),
    lwd = 2, spaceHeight = .7,
    x = 0.25, y = 0.25, height = 1.5,
    just = c("left", "top"), default.units = "inches"
)

annoGenomeLabel(plot = bedpePlot, x = 0.25, y = 1.78, scale = "Mb")

annoHeatmapLegend(
    plot = bedpePlot, fontcolor = "black",
    x = 7.0, y = 0.25,
    width = 0.10, height = 1, fontsize = 10
)

plotText(
    label = "Kb", rot = 90, x = 6.9, y = 0.75,
    just = c("center", "center"), fontsize = 10
)

dev.off()

#####################
library(plotgardener)
library(dplyr)
library(ggplot2)

pchic <- read.table("inputs/HPC7_Promoter_Capture_Interactions.ibed", header = TRUE, sep = "\t")
# hk27ac <- read.table("inputs/GSM1329815_ChIPseq_H3K27Ac_HPC7.cleaned.bedgraph", header = FALSE)

gata2_bait <- pchic[grep("Gata2", pchic$bait_name), ]

gata2_bedpe <- gata2_bait[, c("bait_chr", "bait_start", "bait_end",
                              "otherEnd_chr", "otherEnd_start", "otherEnd_end")]

gata2_bedpe$bait_chr <- paste0("chr", gata2_bedpe$bait_chr)
gata2_bedpe$otherEnd_chr <- paste0("chr", gata2_bedpe$otherEnd_chr)

gata2_bedpe$length <- (gata2_bedpe$otherEnd_start - gata2_bedpe$bait_start) / 1000
gata2_bedpe$other_len <- (gata2_bedpe$otherEnd_end - gata2_bedpe$otherEnd_start) / 1000

pdf("test2.pdf")

chrom <- "chr6"
region_start <- 88120000
region_end   <- 88480000

params <- pgParams(
  chrom = chrom,
  chromstart = region_start,
  chromend = region_end,
  assembly = "mm9",
  x = 0.5,
  y = 0.5,
  width = 6,
  height = 1,
  just = c("left", "top"),
  default.units = "inches"
)

pageCreate(width = 7, height = 2.1, default.units = "inches")

bedpePlot <- plotPairs(
    data = gata2_bedpe, params = params,
    fill = colorby("length", palette = 
                colorRampPalette(c("dodgerblue2", "firebrick2"))),
    lwd = 2, spaceHeight = .7,
    x = 0.5, y = 0.15, height = 1.5,
    just = c("left", "top"), default.units = "inches"
)

annoGenomeLabel(plot = bedpePlot, x = 0.5, y = 1.78, scale = "Mb")

annoHeatmapLegend(
    plot = bedpePlot, fontcolor = "black",
    x = 6.5, y = 0.15,
    width = 0.10, height = 1, fontsize = 10
)

plotText(
    label = "Kb", rot = 90, x = 6.4, y = 0.65,
    just = c("center", "center"), fontsize = 10
)

gene_track <- plotGenes(
  params = params, 
  x = 0.5,
  y = 2.35,
  width = 6,
  height = 0.3,
  fontsize = 6,
  strandLabels = FALSE,
  just = "left"
)

annoGenomeLabel(
  plot = gene_track,
  params = params,
  x = 0.5,
  y = 2.85,
  scale = "Kb",
  fontsize = 7,
  just = "left"
)

gata2_bait_range <- data.frame(
  chrom = "chr6",
  start = 88137173,
  end = 88146476
)

plotRanges(
  data = gata2_bait_range,
  chrom = "chr6",
  chromstart = 88120000,
  chromend = 88480000,
  assembly = "mm9",
  x = 0.5, y = 0.15,
  width = 6,
  height = 0.1,
  fill = "black",
  just = c("left", "top"),
  default.units = "inches"
)

plotText(
  label = "Gata2 bait",
  params = params,
  x = 0.5,
  y = 0.075,
  fontsize = 8,
  just = c("left", "top")
)

dev.off()





