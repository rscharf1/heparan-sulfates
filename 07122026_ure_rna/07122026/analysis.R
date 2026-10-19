library(data.table)
library(dplyr)
library(limma)
library(ggplot2)
library(ggrepel)

# cat GSE5654_series_matrix.txt | grep -v "^!" > mat_clean.txt

counts <- fread("inputs/mat_clean.txt")
map <- fread("inputs/probe_gene_map.txt") %>% .[, c("ID","Gene Symbol")]

counts <- map[counts, on = c(ID = "ID_REF")]

HS_genes <- fread("inputs/HS_relevant_genes.csv") %>% .$hgnc_symbol

table(tolower(HS_genes) %in% tolower(counts2$`Gene Symbol`))


counts_clean <- counts[
  !is.na(`Gene Symbol`) &
  trimws(`Gene Symbol`) != ""
]

sample_cols <- c(
  "GSM132382", "GSM132383", "GSM132384",
  "GSM132385", "GSM132386", "GSM132387"
)

# Expression matrix
expr <- as.matrix(counts_clean[, ..sample_cols])
mode(expr) <- "numeric"

# Shift because dChip values include negatives, then log2-transform
offset <- 1 - min(expr, na.rm = TRUE)
expr_log2 <- log2(expr + offset)

# Collapse probes with the same gene symbol by averaging
expr_gene <- avereps(
  expr_log2,
  ID = counts_clean$`Gene Symbol`
)

## Experimental groups
group <- factor(c(
  "Healthy", "Healthy", "Healthy",
  "Mutant", "Mutant", "Mutant"
))

## Design matrix
design <- model.matrix(~ group)

## Differential expression
fit <- lmFit(expr_gene, design)
fit <- eBayes(fit)

res <- topTable(
  fit,
  coef = "groupMutant",
  number = Inf,
  adjust.method = "BH"
)

res <- as.data.table(res, keep.rownames = "Gene")



# Volcano

HS_core <- c(
  "XYLT1", "XYLT2",
  "B4GALT7",
  "B3GALT6",
  "B3GAT1", "B3GAT2", "B3GAT3",
  "EXT1", "EXT2",
  "EXTL1", "EXTL2", "EXTL3",
  "FAM20B",
  "GLCE",
  "NDST1", "NDST2", "NDST3", "NDST4",
  "HS2ST1",
  "HS6ST1", "HS6ST2", "HS6ST3",
  "HS3ST1", "HS3ST2", "HS3ST3A1", "HS3ST3B1",
  "HS3ST4", "HS3ST5", "HS3ST6"
)

HS_sulfation <- c(
  "PAPSS1",
  "PAPSS2",
  "SLC35B2",
  "SLC35B3",
  "NDST1", "NDST2", "NDST3", "NDST4",
  "HS2ST1",
  "HS6ST1", "HS6ST2", "HS6ST3",
  "HS3ST1", "HS3ST2", "HS3ST3A1",
  "HS3ST3B1", "HS3ST4", "HS3ST5", "HS3ST6",
  "SULF1",
  "SULF2"
)

HS_remodeling <- c(
  "HPSE",
  "HPSE2",
  "SULF1",
  "SULF2"
)

HS_proteoglycans <- c(
  "SDC1", "SDC2", "SDC3", "SDC4",
  "GPC1", "GPC2", "GPC3",
  "GPC4", "GPC5", "GPC6",
  "HSPG2",
  "AGRN",
  "COL18A1"
)



res[, Gene_upper := toupper(Gene)]

# res[, Gene_set := fcase(
#   Gene_upper %in% toupper(HS_core),         "Core biosynthesis",
#   Gene_upper %in% toupper(HS_sulfation),    "Sulfation",
#   Gene_upper %in% toupper(HS_remodeling),   "Remodeling",
#   Gene_upper %in% toupper(HS_proteoglycans),"Proteoglycans",
#   default = "Other"
# )]



make_volcano <- function(res, genes, plot_title) {

  plot_genes <- res[
    Gene_upper %in% toupper(genes)
  ]

  spi1 <- res[
    Gene_upper == "SPI1"
  ]

  ggplot(res, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(
      color = "grey75",
      alpha = 0.7,
      size = 1.2
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    geom_point(
      data = plot_genes,
      color = "red",
      size = 2.3
    ) +
    geom_point(
      data = spi1,
      color = "blue",
      size = 3
    ) +
    geom_text_repel(
      data = plot_genes,
      aes(label = Gene),
      color = "red",
      size = 3.5,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.2
    ) +
    geom_text_repel(
      data = spi1,
      aes(label = Gene),
      color = "blue",
      fontface = "bold",
      size = 4,
      max.overlaps = Inf,
      box.padding = 0.5
    ) +
    theme_classic() +
    labs(
      title = plot_title,
      subtitle = paste0(nrow(plot_genes), " highlighted genes"),
      x = expression(log[2]~fold~change),
      y = expression(-log[10]~italic(P))
    )
}

p_core <- make_volcano(
  res,
  HS_core,
  "Core HS biosynthesis genes"
)

p_sulfation <- make_volcano(
  res,
  HS_sulfation,
  "HS sulfation genes"
)

p_remodeling <- make_volcano(
  res,
  HS_remodeling,
  "HS remodeling genes"
)

p_proteoglycans <- make_volcano(
  res,
  HS_proteoglycans,
  "Heparan sulfate proteoglycans"
)

pdf("07122026/volcano2.pdf", width = 5, height = 4)

p_core
p_sulfation
p_remodeling
p_proteoglycans

dev.off()

# Where these genes fall expression-wise 

counts_clean[tolower(`Gene Symbol`) %in% tolower(HS_core)]

counts_clean[`Gene Symbol` == "Spi1"]

healthy_samples <- c(
  "GSM132382",
  "GSM132383",
  "GSM132384"
)

# Average expression across healthy samples
baseline <- data.table(
  Gene = rownames(expr_gene),
  MeanExpression = rowMeans(
    expr_gene[, healthy_samples, drop = FALSE],
    na.rm = TRUE
  )
)

baseline[, Gene_upper := toupper(Gene)]

genes_to_show <- unique(c(HS_core, "SPI1"))

highlight_baseline <- baseline[
  Gene_upper %in% toupper(genes_to_show)
]

baseline <- baseline[order(MeanExpression)]
baseline[, Rank := .I]

pdf("07122026/expr.pdf")

ggplot(baseline, aes(Rank, MeanExpression)) +
    geom_line(color="grey70") +
    geom_point(
        data=highlight_baseline,
        color="red",
        size=2.5
    ) +
    geom_text_repel(
        data=highlight_baseline,
        aes(label=Gene),
        size=3
    ) +
    theme_classic()

dev.off()




healthy_samples <- c(
  "GSM132382",
  "GSM132383",
  "GSM132384"
)

# Average baseline expression across healthy samples
baseline <- data.table(
  Gene = rownames(expr_gene),
  MeanExpression = rowMeans(
    expr_gene[, healthy_samples, drop = FALSE],
    na.rm = TRUE
  )
)

baseline[, Gene_upper := toupper(Gene)]

# Rank all genes from lowest to highest baseline expression
setorder(baseline, MeanExpression)
baseline[, Rank := .I]

make_expression_rank_plot <- function(baseline, genes, plot_title) {

  plot_genes <- baseline[
    Gene_upper %in% toupper(genes)
  ]

  spi1 <- baseline[
    Gene_upper == "SPI1"
  ]

  ggplot(baseline, aes(x = Rank, y = MeanExpression)) +
    geom_line(
      color = "grey70",
      linewidth = 0.7
    ) +
    geom_point(
      data = plot_genes,
      color = "red",
      size = 2.5
    ) +
    geom_point(
      data = spi1,
      color = "blue",
      size = 3
    ) +
    geom_text_repel(
      data = plot_genes,
      aes(label = Gene),
      color = "red",
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.2,
      min.segment.length = 0
    ) +
    geom_text_repel(
      data = spi1,
      aes(label = Gene),
      color = "blue",
      fontface = "bold",
      size = 3.5,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.2,
      min.segment.length = 0
    ) +
    theme_classic() +
    labs(
      title = plot_title,
      subtitle = paste0(
        nrow(plot_genes),
        " genes detected; average expression across healthy samples"
      ),
      x = "Gene expression rank",
      y = "Mean log2 expression"
    )
}

p_core <- make_expression_rank_plot(
  baseline,
  HS_core,
  "Baseline expression of core HS biosynthesis genes"
)

p_sulfation <- make_expression_rank_plot(
  baseline,
  HS_sulfation,
  "Baseline expression of HS sulfation genes"
)

p_remodeling <- make_expression_rank_plot(
  baseline,
  HS_remodeling,
  "Baseline expression of HS remodeling genes"
)

p_proteoglycans <- make_expression_rank_plot(
  baseline,
  HS_proteoglycans,
  "Baseline expression of heparan sulfate proteoglycans"
)

pdf("07122026/expr2.pdf", width = 5, height = 4)

p_core
p_sulfation
p_remodeling
p_proteoglycans

dev.off()



