suppressPackageStartupMessages( {
  library(forcats)
  library(ggplot2)
  library(GEOquery)
  library(preprocessCore) # quantile normalization
  library(tidyverse)
  library(viridis)
  library(data.table)
  library(limma)
  library(ggrepel)
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

# 1. Sample metadata
group <- ifelse(grepl("^LSC", colnames(expr_mat_norm)), "LSC", "HSC")
group <- factor(group, levels = c("HSC", "LSC"))  # HSC as reference

# 2. Design matrix for limma
design <- model.matrix(~ group)

# 3. Fit model
fit <- lmFit(expr_mat_norm, design)

# 4. Apply empirical Bayes moderation
fit <- eBayes(fit)

topTable_results <- topTable(fit, coef = "groupLSC", number = Inf, adjust.method = "BH") %>%
  as.data.table(., keep.rownames = "probe_id")

topTable_results[adj.P.Val < 0.05][order(-logFC)]

# PROBE TO GENE MAP 
gpl <- getGEO("GPL6106", AnnotGPL = TRUE)
probe_gene_dt <- Table(gpl) %>% data.table()
probe_gene_dt <- probe_gene_dt[, c("ID", "Symbol"), with = FALSE]
probe_gene_dt <- probe_gene_dt[Symbol != ""]
probe_gene_dt$ID <- as.character(probe_gene_dt$ID)

diff <- merge(topTable_results, probe_gene_dt, by.x = "probe_id", by.y = "ID")
diff[, neg_log10_p := -log10(P.Value)]
diff[, sig := ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Significant", "Not Significant")]


a <- ggplot(diff, aes(x = logFC, y = neg_log10_p)) +
  geom_point(aes(color = sig), alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: LSC vs HSC",
    x = "Log2 Fold Change",
    y = "-log10(p-value)",
    color = "Significance"
  )

pos <- c(
  "AKR1C3", "CD34", "GPR56", "SOCS2", "DNMT3B", "SPINK2", "NGFRAP1", "MMRN1",
  "KIAA0125", "EMP1", "LAPTM4B", "CPXM1", "NYNRIN", "MEF2C", "FLT3", "PRDM16"
)

neg <- c("ACTB", "GAPDH", "B2M", "RPL13A", "HPRT1", "PPIA", "TBP", "GUSB", "SDHA")

diff[, control_type := fcase(
  Symbol %in% pos, "Positive Control",
  Symbol %in% neg, "Negative Control",
  default = "Other"
)]


b <- ggplot() +
  geom_point(data = diff[control_type == "Other"], aes(x = logFC, y = neg_log10_p), 
             color = "gray", alpha = 0.5) +
  geom_point(data = diff[control_type != "Other"], 
             aes(x = logFC, y = neg_log10_p, color = control_type),
             size = 2, alpha = 0.9) +
  scale_color_manual(values = c(
    "Positive Control" = "red",
    "Negative Control" = "blue"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: LSC vs HSC",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Control Type"
  )

hs_genes <- fread("inputs/HS_relevant_genes.csv")

diff[, control_type := fcase(
  Symbol %in% pos, "Positive Control",
  Symbol %in% neg, "Negative Control",
  Symbol %in% hs_genes$hgnc_symbol, "HS Gene",
  default = "Other"
)]

c <- ggplot() +
  geom_point(data = diff[control_type == "Other"], aes(x = logFC, y = neg_log10_p), 
             color = "gray", alpha = 0.5) +
  geom_point(data = diff[control_type != "Other"], 
             aes(x = logFC, y = neg_log10_p, color = control_type),
             size = 2, alpha = 0.9) +
  scale_color_manual(values = c(
    "Positive Control" = "red",
    "Negative Control" = "blue",
    "HS Gene" = "black"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: LSC vs HSC",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Control Type"
  ) +
  geom_text_repel(
    data = diff[control_type != "HS Gene"],
    aes(label = Symbol),
    size = 3,
    max.overlaps = 10
  )

pdf("outputs/diff-exp.pdf")
  print(a)
  print(b)
  print(c)
dev.off()






