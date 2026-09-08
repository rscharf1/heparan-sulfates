library(data.table)
library(dplyr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)

# Use the path to your narrowPeak file

peak_files <- list.files(path = "outputs", pattern = "\\.trimmed.bed$", full.names = TRUE)

names(peak_files) <- sub(".*_(AML[0-9]+-[A-Z]+)_.*", "\\1", peak_files)

sample_loc <- peak_files[x]
sample_name <- names(peak_files[x])[1]

parallel::mclapply(seq_along(peak_files), function(x) {
  go(peak_files[x], names(peak_files[x])[1])
}, mc.cores = 4)

go <- function(sample_loc, sample_name) {

  # Load peaks
  peak <- readPeakFile(sample_loc)

  # Load genome annotation
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

  # Annotate peaks
  peak_anno <- annotatePeak(peak, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")

  # Save the annotated table 
  peak_dt <- data.frame(peak_anno) %>% data.table()

  # Ways to find significant peaks
    # Score in V5
    # Where it is binding relative to the gene 
      # Promoter (priority)
      # Enhancer
      # Proximity to TSS 
  # PLOT B
  # Peaks in the promoter region
  promoter_peaks <- peak_dt[grepl("Promoter", annotation)]

  promoter_peaks$TSS_bin <- cut(promoter_peaks$distanceToTSS,
                                 breaks = c(-3000, -1000, -500, 0, 500, 1000, 3000),
                                 labels = c("-3kb to -1kb", "-1kb to -500", "-500 to 0",
                                            "0 to 500", "500 to 1kb", "1kb to 3kb"))


  b <- ggplot(promoter_peaks[V5 > 100], aes(x = TSS_bin, y = as.numeric(V5))) +
    geom_violin(fill = "skyblue", alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    labs(
      x = "Distance to TSS (binned)",
      y = "PU.1 Peak Score (V5)",
      title = "PU.1 Binding Strength Across TSS Proximity Bins"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # pdf("outputs/b.pdf")
  # print(b)
  # dev.off()


  # PLOT C
  controls_pos <- c("SPI1", "Csf1r", "Itgam", "Gfi1", "Cebpa", "Mafb", "Il1b", "Cd14", "Runx1", "Irf8") %>% toupper()
  controls_neg <- c("Ins", "Pax6", "MyoD1", "Otx2", "Tnnt2", "Foxg1", "Fabp4", "Nkx2-5") %>% toupper()

  promoter_peaks$my_label <- "None"
  promoter_peaks[SYMBOL %in% controls_pos]$my_label <- "Pos"
  promoter_peaks[SYMBOL %in% controls_neg]$my_label <- "Neg"

  # top 10k
  promoter_peaks <- promoter_peaks[order(-V5)][1:10000]



  c <- ggplot(promoter_peaks, aes(
      x = V5,
      y = distanceToTSS,
      color = my_label,
      size = my_label,
      alpha = my_label
  )) +
    # geom_point() +
    geom_point(data = m[my_label == "None"], aes(x = V5, y = distanceToTSS), color = "gray80", size = 1, alpha = 0.2) +
    geom_point(data = m[my_label != "None"], aes(x = V5, y = distanceToTSS, color = my_label, size = my_label, alpha = my_label)) + 
    scale_color_manual(
      values = c("None" = "gray80", "Pos" = "red3", "Neg" = "blue3"),
      name = "Control Type"
    ) +
    scale_size_manual(
      values = c("None" = 1, "Pos" = 3, "Neg" = 3),
      guide = "none"
    ) +
    scale_alpha_manual(
      values = c("None" = 0.2, "Pos" = 1, "Neg" = 1),
      guide = "none"
    ) +
    scale_x_continuous(trans = "log10", limits = c(min(promoter_peaks$V5), 1000)) +
    geom_hline(yintercept = c(-1000, 1000), linetype = "dashed", color = "black") +
    theme_minimal() +
    theme(legend.position = "top") + 
    geom_text_repel(
      data = promoter_peaks[my_label != "None"],
      aes(label = SYMBOL),
      size = 3,
      max.overlaps = 10
    ) +
    labs(
      x = "Relative Signal Strength", 
      y = "Distance to TSS"
    )

  # pdf("outputs/c.pdf")
  # print(c)
  # dev.off()

  # PLOT #
  # Plot where HS relevant genes fall on that plot 
  hs_genes <- fread("inputs/HS_relevant_genes.csv")

  m <- merge(promoter_peaks, hs_genes, by.x = "ENSEMBL", by.y = "ensembl_gene_id", all.x=TRUE)

  m[hgnc_symbol %in% hs_genes$hgnc_symbol]$my_label <- "HS"

  d <- ggplot(m, aes(
      x = V5,
      y = distanceToTSS,
      color = my_label,
      size = my_label,
      alpha = my_label
  )) +
    # geom_point() +
    geom_point(data = m[my_label == "None"], aes(x = V5, y = distanceToTSS), color = "gray80", size = 1, alpha = 0.2) +
    geom_point(data = m[my_label != "None"], aes(x = V5, y = distanceToTSS, color = my_label, size = my_label, alpha = my_label)) + 
    scale_color_manual(
      values = c("None" = "gray80", "Pos" = "red3", "Neg" = "blue3", "HS" = "black"),
      name = "Control Type"
    ) +
    scale_size_manual(
      values = c("None" = 1, "Pos" = 3, "Neg" = 3, "HS" = 3),
      guide = "none"
    ) +
    scale_alpha_manual(
      values = c("None" = 0.2, "Pos" = 1, "Neg" = 1, "HS" = 1),
      guide = "none"
    ) +
    scale_x_continuous(trans = "log10", limits = c(min(promoter_peaks$V5), 1000)) +
    geom_hline(yintercept = c(-1000, 1000), linetype = "dashed", color = "black") +
    theme_minimal() +
    theme(legend.position = "top") + 
    geom_text_repel(
      data = m[my_label != "None"],
      aes(label = SYMBOL),
      size = 3,
      max.overlaps = 10
    ) +
    labs(
      x = "Relative Signal Strength", 
      y = "Distance to TSS",
      title = sample_name)

  # pdf("outputs/d.pdf")
  # print(d)
  # dev.off()

  pdf(paste0("outputs/", sample_name, ".pdf"))
    print(b)
    print(c)
    print(d)
  dev.off()  
}


# NEXT
# Visualize in IGV