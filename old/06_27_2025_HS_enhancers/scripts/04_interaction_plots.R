library(plotgardener)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(data.table)

pchic <- read.table("inputs/HPC7_Promoter_Capture_Interactions.ibed", header = TRUE, sep = "\t")
# hk27ac <- read.table("inputs/GSM1329815_ChIPseq_H3K27Ac_HPC7.cleaned.bedgraph", header = FALSE)


# GATA2 EXAMPLE

gata2_bait <- pchic[grep("Gata2", pchic$bait_name), ]

gata2_bedpe <- gata2_bait[, c("bait_chr", "bait_start", "bait_end",
                              "otherEnd_chr", "otherEnd_start", "otherEnd_end")]

gata2_bedpe$bait_chr <- paste0("chr", gata2_bedpe$bait_chr)
gata2_bedpe$otherEnd_chr <- paste0("chr", gata2_bedpe$otherEnd_chr)

gata2_bedpe$length <- (gata2_bedpe$otherEnd_start - gata2_bedpe$bait_start) / 1000
gata2_bedpe$other_len <- (gata2_bedpe$otherEnd_end - gata2_bedpe$otherEnd_start) / 1000

enhancers <- fread("inputs/known_enhancers.csv")
enhancers[gene == "Gata2"]
enhancers$start <- enhancers$middle - 500
enhancers$end <- enhancers$middle + 500

pdf("test2.pdf")

chrom <- "chr6"
region_start <- 88050000
region_end   <- 88480000

params <- pgParams(
  chrom = chrom,
  chromstart = region_start,
  chromend = region_end,
  assembly = "mm9",
  x = 0.5,
  y = 0.5,
  width = 4,
  height = 5,
  just = c("left", "top"),
  default.units = "inches"
)

pageCreate(width = 5, height = 4.1, default.units = "inches", showGuides=TRUE)

# PCHi-C interactions
bedpePlot <- plotPairs(
    data = gata2_bedpe, params = params,
    fill = colorby("length", palette = 
                colorRampPalette(c("dodgerblue2", "firebrick2"))),
    lwd = 2, spaceHeight = .7,
    x = 0.5, y = 0.1, height = 1.5,
    just = c("left", "top"), default.units = "inches"
)

# annoGenomeLabel(plot = bedpePlot, x = 0.5, y = 1.78, scale = "Mb")

annoHeatmapLegend(
    plot = bedpePlot, fontcolor = "black",
    x = 4.5, y = 0.15,
    width = 0.10, height = 1, fontsize = 10
)

plotText(
    label = "Kb", rot = 90, x = 4.4, y = 0.65,
    just = c("center", "center"), fontsize = 10
)

# Enhancers
enhancer_gr <- GRanges(
  seqnames = enhancers$chrom,
  ranges = IRanges(start = enhancers$start, end = enhancers$end)
)

plotRanges(
  data = enhancer_gr,
  chrom = chrom,
  chromstart = region_start,
  chromend = region_end,
  assembly = "mm9",
  x = 0.5,
  y = 1.9,
  width = 4,
  height = 0.08,
  fill = "firebrick3",
  just = c("left", "top"),
  default.units = "inches"
)

# Add enhancer label
plotText(
  label = "-77 kb enhancer",
  params = params,
  x = enhancers$start,
  y = 1.65,
  fontsize = 7,
  just = c("left", "top")
)

# Genes
gene_track <- plotGenes(
  params = params, 
  x = 0.5,
  y = 2.35,
  width = 4,
  height = 0.3,
  fontsize = 8.5,
  strandLabels = TRUE,
  just = "left"
)

annoGenomeLabel(
  plot = gene_track,
  params = params,
  x = 0.5,
  y = 2.65,
  scale = "Kb",
  # fontsize = 7,
  just = "left"
)

# bait
gata2_bait_range <- data.frame(
  chrom = "chr6",
  start = 88137173,
  end = 88146476
)

plotRanges(
  data = gata2_bait_range,
  chrom = "chr6",
  chromstart = 88050000,
  chromend = 88480000,
  assembly = "mm9",
  x = 0.5, y = 0.15,
  width = 4,
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

plotText(
  label = "PCHi-C Interactions",
  params = params,
  x = 0.5,
  y = 0.3,
  fontsize = 8,
  just = c("left", "top")
)

plotText(
  label = "Known Enhancers",
  params = params,
  x = 0.5,
  y = 1.65,
  fontsize = 8,
  just = c("left", "top")
)

dev.off()

######################################
library(plotgardener)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(data.table)

pchic <- read.table("inputs/HPC7_Promoter_Capture_Interactions.ibed", header = TRUE, sep = "\t")

genes <- c("gata2", "sfpi1", "runx1", "myc-")

for(g in genes) {
  bait <- pchic[grep(g, pchic$bait_name, ignore.case = TRUE), ]

  g <- gsub("-", "", g)

  # bait <- bait[1:3, ]

  bedpe <- bait[, c("bait_chr", "bait_start", "bait_end",
                              "otherEnd_chr", "otherEnd_start", "otherEnd_end")]

  bedpe$bait_chr <- paste0("chr", bedpe$bait_chr)
  bedpe$otherEnd_chr <- paste0("chr", bedpe$otherEnd_chr)

  bedpe$length <- (bedpe$otherEnd_start - bedpe$bait_start) / 1000
  bedpe$other_len <- (bedpe$otherEnd_end - bedpe$otherEnd_start) / 1000

  enhancers <- fread("inputs/known_enhancers.csv")
  enhancers <- enhancers[grep(g, gene, ignore.case=TRUE)]

  enhancers$start <- enhancers$middle - 500
  enhancers$end <- enhancers$middle + 500

  pdf(paste0("outputs/", g, ".pdf"))

  chrom <- bedpe$bait_chr[1]
  region_start <- c(bedpe$bait_start, bedpe$otherEnd_start, enhancers$start) %>% min() %>% (function(x) x/10000) %>% floor() %>% (function(x) x*10000)
  region_start <- region_start - 1e5
  region_end <- c(bedpe$bait_end, bedpe$otherEnd_end, enhancers$end) %>% max() %>% (function(x) x/10000) %>% ceiling() %>% (function(x) x*10000)
  region_end <- region_end + 1e5

  params <- pgParams(
    chrom = chrom,
    chromstart = region_start,
    chromend = region_end,
    assembly = "mm9",
    x = 0.5,
    y = 0.5,
    width = 4,
    height = 5,
    just = c("left", "top"),
    default.units = "inches"
  )

  pageCreate(width = 5, height = 1.5 + nrow(bait) / 6, default.units = "inches", showGuides=TRUE)

  # bait
  bait_range <- data.frame(
    chrom = chrom,
    start = bedpe$bait_start[1],
    end = bedpe$bait_end[1]
  )

  cur_y <- 0.075

  plotText(
    label = "Bait",
    params = params,
    x = 0.5,
    y = cur_y, #0.075
    fontsize = 8,
    just = c("left", "top")
  )

  cur_y <- cur_y + 0.075

  plotRanges(
    data = bait_range,
    chrom = chrom,
    chromstart = region_start,
    chromend = region_end,
    assembly = "mm9",
    x = 0.5, y = cur_y,
    width = 4,
    height = 0.1,
    fill = "black",
    just = c("left", "top"),
    default.units = "inches"
  )

  cur_y <- cur_y + 0.15

  # PCHi-C interactions
  plotText(
    label = "PCHi-C Interactions",
    params = params,
    x = 0.5,
    y = cur_y,
    fontsize = 8,
    just = c("left", "top")
  )

  bedpePlot <- plotPairs(
      data = bedpe, params = params,
      fill = colorby("length", palette = 
                  colorRampPalette(c("dodgerblue2", "firebrick2"))),
      lwd = 2, spaceHeight = .7,
      x = 0.5, y = cur_y, height = nrow(bedpe) / 7, #2
      just = c("left", "top"), default.units = "inches"
  )
  annoHeatmapLegend(
      plot = bedpePlot, fontcolor = "black",
      x = 4.5, y = cur_y,
      width = 0.10, height = 1, fontsize = 10
  )

  plotText(
      label = "Kb", rot = 90, x = 4.4, y = cur_y + 0.4,
      just = c("center", "center"), fontsize = 10
  )

  cur_y <- cur_y + nrow(bedpe) / 7 + 0.2

  # Enhancers

  plotText(
    label = "Known Enhancers",
    params = params,
    x = 0.5,
    y = cur_y,
    fontsize = 8,
    just = c("left", "top")
  )

  cur_y <- cur_y + 0.15

  enhancer_gr <- GRanges(
    seqnames = enhancers$chrom,
    ranges = IRanges(start = enhancers$start, end = enhancers$end)
  )

  plotRanges(
    data = enhancer_gr,
    chrom = chrom,
    chromstart = region_start,
    chromend = region_end,
    assembly = "mm9",
    x = 0.5,
    y = cur_y,
    width = 4,
    height = 0.08,
    fill = "firebrick3",
    just = c("left", "top"),
    default.units = "inches"
  )

  cur_y <- cur_y + 0.3

  # Genes
  gene_track <- plotGenes(
    params = params, 
    x = 0.5,
    y = cur_y,
    width = 4,
    height = 0.3,
    fontsize = 8.5,
    strandLabels = TRUE,
    just = "left"
  )

  cur_y <- cur_y + 0.35

  annoGenomeLabel(
    plot = gene_track,
    params = params,
    x = 0.5,
    y = cur_y,
    scale = "Kb",
    # fontsize = 7,
    just = "left"
  )

  dev.off()

}























