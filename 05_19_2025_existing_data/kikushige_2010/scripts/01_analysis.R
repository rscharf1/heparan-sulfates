suppressPackageStartupMessages( {
  library(forcats)
  library(ggplot2)
  library(GEOquery)
  library(preprocessCore) # quantile normalization
  library(tidyverse)
  library(viridis)
  library(data.table)
} )

gse <- getGEO("GSE24395",GSEMatrix=TRUE)

expr_mat <- exprs(gse[[1]])
sample_info <- pData(gse[[1]]) %>% data.table()
probe_info <- fData(gse[[1]]) %>% data.table()

# NORMALIZE
expr_mat[1:4,1:3]

probe_ids <- rownames(expr_mat)
sample_ids <- colnames(expr_mat)

expr_mat_norm <- normalize.quantiles(expr_mat, copy = TRUE)
expr_mat_norm[1:4,1:3]

colnames(expr_mat_norm) <- sample_ids
rownames(expr_mat_norm) <- probe_ids

# SAMPLES
sample_info[, !c("extract_protocol_ch1", "label_protocol_ch1", "data_processing"), with = FALSE]
sample_info <- sample_info[, c("geo_accession", "title", "disease:ch1"), with = FALSE]

name_map <- setNames(sample_info$title, sample_info$geo_accession)
name_map <- name_map %>% gsub("normal", "", .) %>% trimws(.)
colnames(expr_mat_norm) <- name_map[colnames(expr_mat_norm)]
expr_mat_norm[1:4,1:3]



# PROBE TO GENE MAP 
gpl <- getGEO("GPL6106", AnnotGPL = TRUE)
probe_gene_dt <- Table(gpl) %>% data.table()
probe_gene_dt <- probe_gene_dt[, c("ID", "Symbol"), with = FALSE]
probe_gene_dt <- probe_gene_dt[Symbol != ""]
probe_gene_dt$ID <- as.character(probe_gene_dt$ID)

# Control genes? 
  # Pos: genes we know will differ between HSCs and LSCs
  # Neg: genes that shouldn't change much 

pos <- c(
  "AKR1C3", "CD34", "GPR56", "SOCS2", "DNMT3B", "SPINK2", "NGFRAP1", "MMRN1",
  "KIAA0125", "EMP1", "LAPTM4B", "CPXM1", "NYNRIN", "MEF2C", "FLT3", "PRDM16"
)

neg <- c("ACTB", "GAPDH", "B2M", "RPL13A", "HPRT1", "PPIA", "TBP", "GUSB", "SDHA")

hs_genes <- fread("inputs/HS_relevant_genes.csv")

table(pos %in% probe_gene_dt$Symbol)
table(neg %in% probe_gene_dt$Symbol)

expr_mat_norm[c("4070341", "6290634"), ]

expr_long <- expr_mat_norm[probe_gene_dt[Symbol %in% c(pos, neg, hs_genes$hgnc_symbol)]$ID %>% as.character(),] %>% 
  as.data.table(., keep.rownames = "Gene") %>% 
  melt(
    .,
    id.vars = "Gene",
    variable.name = "Sample",
    value.name = "Expression"
  )

expr_long <- merge(expr_long, probe_gene_dt, by.x="Gene", by.y="ID")

expr_long[, Group := ifelse(grepl("^LSC", Sample), "LSC", "HSC")]

# POS
probes <- expr_long[Symbol %in% pos]$Gene %>% unique()

pdf("outputs/pos-ctrl.pdf")

lapply(seq_along(probes), function(x) {
  gene <- expr_long[Gene == probes[x]]$Symbol[1]

  ggplot(expr_long[Gene == probes[x]], aes(x = Group, y = Expression)) +
    geom_boxplot(aes(fill = Group), alpha = 0.4) +
    geom_jitter(width = 0.2, size = 2) +
    labs(x = "Cell", y = gene, title = paste(gene, "Expression by Cell Type")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  })

dev.off()

# NEG
probes <- expr_long[Symbol %in% neg]$Gene %>% unique()

pdf("outputs/neg-ctrl.pdf")

lapply(seq_along(probes), function(x) {
  gene <- expr_long[Gene == probes[x]]$Symbol[1]

  ggplot(expr_long[Gene == probes[x]], aes(x = Group, y = Expression)) +
    geom_boxplot(aes(fill = Group), alpha = 0.4) +
    geom_jitter(width = 0.2, size = 2) +
    labs(x = "Cell", y = gene, title = paste(gene, "Expression by Cell Type")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  })

dev.off()

# HS
probes <- expr_long[Symbol %in% hs_genes$hgnc_symbol]$Gene %>% unique()

pdf("outputs/hs-genes.pdf")

lapply(seq_along(probes), function(x) {
  gene <- expr_long[Gene == probes[x]]$Symbol[1]

  ggplot(expr_long[Gene == probes[x]], aes(x = Group, y = Expression)) +
    geom_boxplot(aes(fill = Group), alpha = 0.4) +
    geom_jitter(width = 0.2, size = 2) +
    labs(x = "Cell", y = gene, title = paste(gene, "Expression by Cell Type")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  })

dev.off()






