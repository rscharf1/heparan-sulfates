library(GenomicRanges)
library(rtracklayer)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)

# 1. Load PU.1 peaks
pu1_peaks <- import("inputs/PU.1_IP_VEH_URE_peaks.narrowPeak", format = "narrowPeak")

# 2. Load gene annotations
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)

# 3. Method 1: Find genes whose **promoters** overlap PU.1 peaks
promoters_gr <- promoters(genes_gr, upstream = 1000, downstream = 500)
promoter_hits <- findOverlaps(promoters_gr, pu1_peaks)
promoter_gene_ids <- names(promoters_gr)[queryHits(promoter_hits)]

# 4. Method 2: Find the **nearest gene** to each peak
nearest_hits <- nearest(pu1_peaks, genes_gr)
nearest_gene_ids <- names(genes_gr)[nearest_hits]

# 5. Map Entrez IDs to gene symbols
gene_symbols_promoter <- mapIds(org.Hs.eg.db, keys = promoter_gene_ids, column = "SYMBOL", keytype = "ENTREZID")
gene_symbols_nearest <- mapIds(org.Hs.eg.db, keys = nearest_gene_ids, column = "SYMBOL", keytype = "ENTREZID") %>% as.character()

nearest_gene_ids <- names(genes_gr)[nearest_hits]
nearest_gene_symbols <- mapIds(org.Hs.eg.db, keys = nearest_gene_ids, 
                               column = "SYMBOL", keytype = "ENTREZID", multiVals = "first") %>% as.character()

has_promoter_hit <- logical(length(pu1_peaks))
has_promoter_hit[subjectHits(promoter_hits)] <- TRUE

bound_genes_df <- data.frame(
  Peak_ID = pu1_peaks$name,
  Chrom = as.character(seqnames(pu1_peaks)),
  Start = start(pu1_peaks),
  End = end(pu1_peaks),
  NearestGene = as.character(unname(nearest_gene_symbols)),
  OverlapsPromoter = has_promoter_hit,
  SignalValue = pu1_peaks$signalValue,
  PValue = pu1_peaks$pValue,
  QValue = pu1_peaks$qValue,
  PeakSummitOffset = pu1_peaks$peak
)

controls_pos <- c("SPI1", "Csf1r", "Itgam", "Gfi1", "Cebpa", "Mafb", "Il1b", "Cd14", "Runx1", "Irf8") %>% toupper()
controls_neg <- c("Ins", "Pax6", "MyoD1", "Otx2", "Tnnt2", "Foxg1", "Fabp4", "Nkx2-5") %>% toupper()

table(toupper(controls_pos) %in% bound_genes_df$NearestGene)
table(toupper(controls_neg) %in% bound_genes_df$NearestGene)


dt <- bound_genes_df %>% data.table()

dt$my_label <- "None"
dt[NearestGene %in% toupper(controls_pos)]$my_label <- "Pos"
dt[NearestGene %in% toupper(controls_neg)]$my_label <- "Neg"

dt_top_peaks <- dt[ , .SD[which.max(SignalValue)], by = NearestGene]

a <- ggplot(dt_top_peaks, aes(x = PValue, y = SignalValue)) +
  geom_point(alpha = 0.5, size = 1.5) +
  labs(
    title = "Peak Strength vs. Significance",
    x = expression(-log[10](p~value)),
    y = "Signal Value (Fold Enrichment)"
  ) +
  theme_minimal()

b <- ggplot() +
  # First plot "None" points (background)
  geom_point(
    data = dt_top_peaks[my_label == "None"],
    aes(x = PValue, y = SignalValue, color = my_label, size = my_label, alpha = my_label)
  ) +
  # Then plot "Pos" and "Neg" (foreground)
  geom_point(
    data = dt_top_peaks[my_label != "None"],
    aes(x = PValue, y = SignalValue, color = my_label, size = my_label, alpha = my_label)
  ) +
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
  theme_minimal() +
  theme(legend.position = "top") +
  labs(
    x = expression(-log[10](p~value)), 
    y = "Signal Value (Fold Enrichment)"
  ) + 
  geom_text_repel(
    data = dt_top_peaks[my_label != "None"],
    aes(x = PValue, y = SignalValue, label = NearestGene),
    size = 3,
    max.overlaps = 10
  )

hs_genes <- fread("inputs/HS_relevant_genes.csv")


table(hs_genes$hgnc_symbol %in% dt_top_peaks$NearestGene)

dt_top_peaks[NearestGene %in% hs_genes$hgnc_symbol]$my_label <- "HS"

c <- ggplot() +
  # First plot "None" points (background)
  geom_point(
    data = dt_top_peaks[my_label == "None"],
    aes(x = PValue, y = SignalValue, color = my_label, size = my_label, alpha = my_label)
  ) +
  # Then plot "Pos" and "Neg" (foreground)
  geom_point(
    data = dt_top_peaks[my_label != "None"],
    aes(x = PValue, y = SignalValue, color = my_label, size = my_label, alpha = my_label)
  ) +
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
  theme_minimal() +
  theme(legend.position = "top") +
  labs(
    x = expression(-log[10](p~value)), 
    y = "Signal Value (Fold Enrichment)"
  ) + 
  geom_text_repel(
    data = dt_top_peaks[my_label != "None"],
    aes(x = PValue, y = SignalValue, label = NearestGene),
    size = 3,
    max.overlaps = 10
  )

pdf("outputs/out.pdf")
  print(a)
  print(b)
  print(c)
dev.off()


