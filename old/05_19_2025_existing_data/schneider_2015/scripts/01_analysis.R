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

gse <- getGEO("GSE68172",GSEMatrix=TRUE)

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
expr_mat_norm[1:4,1:3]

#SAMPLES 
sample_info <- sample_info[, c("geo_accession", "disease status:ch1", "cd34 status:ch1", "cd38 status:ch1"), with = FALSE]
colnames(sample_info) <- colnames(sample_info) %>% gsub("status:ch1", "", .) %>% trimws()

sample_info$concat <- paste0(sample_info$disease, "_", sample_info$cd34, "_", sample_info$cd38) %>%
  gsub("\\+|-| ", "", .) %>% make.unique()

sample_info[, disease := factor(disease, levels = c("AML", "healthy"))]

name_map <- setNames(sample_info$concat, sample_info$geo_accession)

colnames(expr_mat_norm) <- name_map[colnames(expr_mat_norm)]

expr_mat_norm[1:4,1:4]

# sample_labels <- colnames(expr_mat_norm)

# sample_meta <- data.table(
#   sample = sample_labels,
#   disease  = sub("_.*$", "", sample_labels),
#   cell_type = sub("^[^_]+_", "", sample_labels)
# )

# # Convert to factor with order
# sample_meta[, disease := factor(disease, levels = c("AML", "healthy"))]
# sample_meta[, cell_type := factor(cell_type, levels = c("plus_minus", "minus_minus", "plus_plus", "minus_plus"))]

gpl <- getGEO("GPL570", AnnotGPL = TRUE)
probe_gene_dt <- Table(gpl) %>% data.table()
probe_gene_dt <- probe_gene_dt[, c("ID", "Gene symbol"), with = FALSE]
probe_gene_dt <- probe_gene_dt[`Gene symbol` != ""]
probe_gene_dt[, `Gene symbol` := tstrsplit(`Gene symbol`, "///", fixed = TRUE)[[1]]]
probe_gene_dt$ID <- as.character(probe_gene_dt$ID)

pos <- c(
  "AKR1C3", "CD34", "GPR56", "SOCS2", "DNMT3B", "SPINK2", "NGFRAP1", "MMRN1",
  "KIAA0125", "EMP1", "LAPTM4B", "CPXM1", "NYNRIN", "MEF2C", "FLT3", "PRDM16"
)

neg <- c("ACTB", "GAPDH", "B2M", "RPL13A", "HPRT1", "PPIA", "TBP", "GUSB", "SDHA")

hs_genes <- fread("inputs/HS_relevant_genes.csv")

table(pos %in% probe_gene_dt$`Gene symbol`)
table(neg %in% probe_gene_dt$`Gene symbol`)
table(hs_genes$hgnc_symbol %in% probe_gene_dt$`Gene symbol`)

# DIFF EXP
samples <- sample_info[grepl("plus", cd34) & grepl("minus", cd38), concat]
expr_samples <- expr_mat_norm[, samples, drop = FALSE]
group <- sample_info[grepl("plus", cd34) & grepl("minus", cd38), disease]

# Design matrix
design <- model.matrix(~ group)  # intercept + leukemia effect
colnames(design)  # "(Intercept)" and "groupLeukemia"

fit <- lmFit(expr_samples, design)
fit <- eBayes(fit)
topTable_results <- topTable(fit, coef = "grouphealthy", number = Inf, adjust.method = "BH") %>%
  as.data.table(., keep.rownames = "probe_id")

diff <- merge(topTable_results, probe_gene_dt, by.x = "probe_id", by.y = "ID")
diff[, neg_log10_p := -log10(P.Value)]
diff[, sig := ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Significant", "Not Significant")]

mytitle <- "Volcano Plot: CD34+/CD38- | healthy / AML"

a <- ggplot(diff, aes(x = logFC, y = neg_log10_p)) +
  geom_point(aes(color = sig), alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    title = mytitle,
    x = "Log2 Fold Change",
    y = "-log10(p-value)",
    color = "Significance"
  )

diff[, control_type := fcase(
  `Gene symbol` %in% pos, "Positive Control",
  `Gene symbol` %in% neg, "Negative Control",
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
    title = mytitle,
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Control Type"
  )

diff[, control_type := fcase(
  `Gene symbol` %in% pos, "Positive Control",
  `Gene symbol` %in% neg, "Negative Control",
  `Gene symbol` %in% hs_genes$hgnc_symbol, "HS Gene",
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
    title = mytitle,
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Control Type"
  ) +
  geom_text_repel(
    data = diff[control_type == "HS Gene"],
    aes(x = logFC, y = neg_log10_p, label = `Gene symbol`),
    size = 3,
    max.overlaps = 10
  )

pdf("outputs/34pos_38neg.pdf")
  print(a)
  print(b)
  print(c)
dev.off()

