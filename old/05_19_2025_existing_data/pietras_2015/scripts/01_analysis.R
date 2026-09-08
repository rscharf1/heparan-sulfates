suppressPackageStartupMessages( {
  library(forcats)
  library(ggplot2)
  library(GEOquery)
  library(pheatmap)
  library(preprocessCore) # quantile normalization
  library(tidyverse)
  library(viridis)
  # library(AnnotationDbi)
  # library(mogene10sttranscriptcluster.db)
  library(data.table)
} )

gse <- getGEO("GSE68529",GSEMatrix=TRUE)

length(gse)
class(gse)

expr <- exprs(gse[[1]])
pheno <- pData(gse[[1]]) %>% data.table()
features <- fData(gse[[1]]) %>% data.table()

# NORMALIZATION
probe_ids <- rownames(expr)
sample_ids <- colnames(expr)
colmeans <- colSums(expr)/nrow(expr)

expr <- data.matrix(expr)
expr <- normalize.quantiles(expr, copy = TRUE)

colnames(expr) <- sample_ids
rownames(expr) <- probe_ids

# SAMPLES
samples <- pheno[, c("geo_accession", "title", "description.1"), with = FALSE]

name_map <- setNames(samples$description.1, samples$geo_accession)
colnames(expr) <- name_map[colnames(expr)]

# PROBE TO GENE MAP 
gpl <- getGEO("GPL6246", AnnotGPL = TRUE)
probe_gene_dt <- Table(gpl) %>% data.table()
probe_gene_dt <- probe_gene_dt[, c(1:3), with = FALSE]
colnames(probe_gene_dt) <- c("ID", "title", "symbol")
probe_gene_dt <- probe_gene_dt[symbol != ""]
probe_gene_dt <- probe_gene_dt[!grepl("///", symbol)]

probe_gene_dt[ID %in% rownames(expr)]

expr[1:4,1:4]

# rownames(expr)

# expr[rownames(expr) %in% probe_gene_dt$ID] %>% dim()

# sample <- expr[head(intersect(probe_gene_dt$ID, rownames(expr))), ]
# probe_gene_dt[ID %in% rownames(sample)]
# genes <- probe_gene_dt[ID %in% rownames(sample)]$symbol

# LOOK AT HS GENES 
hs_genes <- fread("inputs/HS_relevant_genes.csv")

probe_gene_dt[toupper(symbol) %in% hs_genes$hgnc_symbol]$ID

expr_sub <- expr[as.character(probe_gene_dt[toupper(symbol) %in% hs_genes$hgnc_symbol]$ID), ]
rownames(expr_sub) <- probe_gene_dt[toupper(symbol) %in% hs_genes$hgnc_symbol]$symbol

# Convert wide matrix to long format
expr_long <- melt(
  as.data.table(expr_sub, keep.rownames = "Gene"),
  id.vars = "Gene",
  variable.name = "Sample",
  value.name = "Expression"
)

# Add "Cell" group by parsing Sample names
expr_long[, Cell := sub("-\\d+$", "", Sample)]  # Removes "-1", "-2", etc.

gene_to_plot <- "B3gat1"

cell_order <- c("HSCLT", "HSCST", "MPP2", "MPP3", "MPP4", "CMP", "GMP", "PreGr", "Gr")

expr_long[, Cell := factor(Cell, levels = cell_order)]

plot_data <- expr_long[Gene == gene_to_plot]

pdf("outputs/out.pdf")
ggplot(expr_long[Gene == "B3gat1"], aes(x = Cell, y = Expression)) +
  geom_boxplot(aes(fill = Cell), alpha = 0.4) +
  geom_jitter(width = 0.2, size = 2) +
  labs(x = "Cell", y = "B3gat1", title = "B3gat1 Expression by Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

genes <- expr_long$Gene %>% unique()

pdf("outputs/all_genes.pdf")

lapply(seq_along(genes), function(x) {
  ggplot(expr_long[Gene == genes[x]], aes(x = Cell, y = Expression)) +
    geom_boxplot(aes(fill = Cell), alpha = 0.4) +
    geom_jitter(width = 0.2, size = 2) +
    labs(x = "Cell", y = genes[x], title = paste(genes[x], "Expression by Cell Type")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  })

dev.off()



















