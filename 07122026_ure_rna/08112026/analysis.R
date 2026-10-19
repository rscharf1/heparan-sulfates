library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

# cat GSE5654_series_matrix.txt | grep -v "^!" > mat_clean.txt

counts <- fread("inputs/mat_clean.txt")
map <- fread("inputs/probe_gene_map.txt") %>% .[, c("ID","Gene Symbol")]

counts <- map[counts, on = c(ID = "ID_REF")]

HS_genes <- fread("inputs/HS_relevant_genes.csv") %>% .$hgnc_symbol

counts[, WT := rowMeans(.SD),
       .SDcols = c("GSM132382", "GSM132383", "GSM132384")]

counts[, PU1_KD := rowMeans(.SD),
       .SDcols = c("GSM132385", "GSM132386", "GSM132387")]

plot_dt <- counts[
  WT > 0 & PU1_KD > 0
]

plot_dt[, log2_WT := log2(WT+1)]
plot_dt[, log2_PU1_KD := log2(PU1_KD+1)]

pdf("08112026/out1.pdf")

  ggplot(plot_dt, aes(x = log2_WT, y = log2_PU1_KD)) +
    geom_point(
      alpha = 0.25,
      size = 0.8
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2
    ) +
    coord_equal() +
    theme_classic() +
    labs(
      x = "WT HSC expression (log2)",
      y = "PU.1 KD HSC expression (log2)"
    )

dev.off()

plot_dt[grepl("Xylt1", `Gene Symbol`, ignore.case=TRUE)]

to_highlight <- c(
  "PU1"  = "Spi1",
  "EXT1" = "Ext1",
  "EXT2" = "Ext2",
  "EXTL2" = "Extl2",
  "EXTL3" = "Extl3",
  "B3GAT1" = "B3gat1",
  "B3GAT3" = "B3gat3",
  "B4GALT1" = "B4galt1",
  "B4GALT3" = "B4galt3",
  "B4GALT7" = "B4galt7",
  "FAM20B" = "Fam20b",
  "GLCE" = "Glce",
  "XYLT2" = "Xylt2"
)

hl <- data.table(
  label = names(to_highlight),
  `Gene Symbol` = unname(to_highlight)
)

hl <- merge(
  hl,
  plot_dt,
  by = "Gene Symbol"
)


pdf("08112026/out2.pdf")

  ggplot(
    plot_dt,
    aes(x = log2_WT, y = log2_PU1_KD)
  ) +
    
    # All probes
    geom_point(
      alpha = 0.25,
      size = 0.8
    ) +
    
    # Identity line
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2
    ) +
    
    # Genes of interest ON TOP
    geom_point(
      data = hl,
      color = "red",
      size = 3
    ) +
    
    # Labels
    geom_text_repel(
      data = hl,
      aes(label = label),
      color = "red",
      fontface = "bold",
      size = 4,
      max.overlaps = Inf
    ) +
    
    coord_equal() +
    theme_classic() +
    labs(
      x = "WT HSC expression (log2)",
      y = "PU.1 KD HSC expression (log2)"
    )

dev.off()









