library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)

# Load data
expr <- fread("data/rsem_gene_tmp.gz", data.table = FALSE)
pheno <- fread("data/metadata.txt.gz", data.table = FALSE) %>% data.table()

colnames(pheno) <- c("sample", "detailed_site", "primary_disease", "primary_site", "sample_type", "sex", "study")

rownames(expr) <- expr[,1]
expr <- expr[,-1]

hs_genes <- fread("~/steidl_lab/heparan-sulfates/useful_stuff/HS_relevant_genes.csv")
gene_map <- fread("data/gene_id_map.tsv")

length(intersect(hs_genes$hgnc_symbol, gene_map$gene))

setdiff(
	hs_genes$hgnc_symbol,
	intersect(hs_genes$hgnc_symbol, gene_map$gene)
)

# Start with one cancer type and then expand 
	# Plot 1: In tumor and healthy samples, show where these enzymes fall expression wise relative to other genes 
all_gene_expression <- function(tmp_mat, tmp_title) {
	gene_means <- rowMeans(tmp_mat)
	gene_means <- sort(gene_means, decreasing = TRUE)

	df <- data.frame(
	  gene = names(gene_means),
	  mean_expression = unname(gene_means),
	  stringsAsFactors = FALSE
	)

	df$gene <- factor(df$gene, levels = df$gene)

	df$rank <- seq_len(nrow(df))

	df$highlight <- ifelse(df$gene %in% gene_map[gene %in% hs_genes$hgnc_symbol]$id, "HS", "Other")

	df <- df %>%
	  left_join(
	    gene_map[gene %in% hs_genes$hgnc_symbol] %>% select(gene = id, gene_symbol = gene),
	    by = "gene"
	  )

	df %>% filter(!is.na(gene_symbol))


	p <- ggplot(df, aes(x = rank, y = mean_expression)) +
	  geom_point(
	    data = subset(df, highlight == "Other"),
	    color = "lightgray",
	    size = 1
	  ) +
	  geom_point(
	    data = subset(df, highlight == "HS"),
	    color = "red",
	    size = 2
	  ) +
	  geom_text_repel(
	    data = df %>% filter(!is.na(gene_symbol)),
	    aes(label = gene_symbol),
	    size = 3,
	    max.overlaps = 50,
	    box.padding = 0.3,
	    point.padding = 0.2,
	    segment.color = "gray50"
	  ) +
	  theme_bw() +
	  labs(
	    x = "Genes ranked by mean expression",
	    y = "Mean expression",
	    title = tmp_title
	  )

	p
}

nl <- pheno[study == "GTEX"][primary_disease == "Pancreas"]
expr_normal <- expr[, colnames(expr) %in% nl$sample] %>% as.matrix()
mode(expr_normal) <- "numeric"


tm <- pheno[study == "TCGA"][detailed_site == "Pancreatic Adenocarcinoma"][sample_type == "Primary Tumor"]
expr_tumor <- expr[, colnames(expr) %in% tm$sample] %>% as.matrix()
mode(expr_tumor) <- "numeric"

p1 <- all_gene_expression(expr_normal, "Pancreas - Healthy")
p2 <- all_gene_expression(expr_tumor, "Pancreas - Cancer")

pdf("out/01.pdf")
	p1
	p2
dev.off()

	# Plot 2: Volcano plot, highlighting enzymes 
volcano <- function(expr_tumor, expr_normal, tmp_title = "empty title") {

	expr_all <- cbind(expr_tumor, expr_normal)

	keep <- rowSums(expr_all > -2) >= 0.1 * ncol(expr_all)

	expr_tumor <- expr_tumor[keep, ]
	expr_normal <- expr_normal[keep, ]

	logFC <- rowMeans(expr_tumor) - rowMeans(expr_normal)

	pvals <- sapply(seq_len(nrow(expr_tumor)), function(i) {
	  wilcox.test(expr_tumor[i, ], expr_normal[i, ])$p.value
	})

	padj <- p.adjust(pvals, method = "BH")

	res <- data.frame(
	  gene = rownames(expr_tumor),
	  logFC = logFC,
	  pval = pvals,
	  padj = padj,
	  stringsAsFactors = FALSE
	)

	res$group <- "NS"
	res$group[res$padj < 0.05 & res$logFC > 2]  <- "Up"
	res$group[res$padj < 0.05 & res$logFC < -2] <- "Down"

	res$highlight <- ifelse(res$gene %in% gene_map[gene %in% hs_genes$hgnc_symbol]$id, "HS", "Other")

	res <- res %>%
	  left_join(
	    gene_map[gene %in% hs_genes$hgnc_symbol] %>% select(gene = id, gene_symbol = gene),
	    by = "gene"
	  )

	plot_df <- res %>%
	  mutate(
	    padj_plot = padj,
	    padj_plot = ifelse(is.na(padj_plot), 1, padj_plot),
	    padj_plot = ifelse(padj_plot == 0, 1e-300, padj_plot),
	    neglog10_padj = -log10(padj_plot)
	  )

	label_df <- plot_df %>%
	  filter(!is.na(gene_symbol), gene_symbol != "")

	p1 <- ggplot(plot_df, aes(x = logFC, y = neglog10_padj)) +
	  geom_point(
	    data = subset(plot_df, group == "NS"),
	    color = "gray80",
	    size = 1.2,
	    alpha = 0.7
	  ) +
	  geom_point(
	    data = subset(plot_df, group == "Up"),
	    color = "red",
	    size = 1.8,
	    alpha = 0.85
	  ) +
	  geom_point(
	    data = subset(plot_df, group == "Down"),
	    color = "blue",
	    size = 1.8,
	    alpha = 0.85
	  ) +
	  geom_vline(xintercept = c(-2, 2), linetype = "dashed", linewidth = 0.4) +
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
	  geom_text_repel(
	    data = label_df,
	    aes(label = gene_symbol),
	    size = 3,
	    max.overlaps = 50,
	    box.padding = 0.3,
	    point.padding = 0.2,
	    segment.color = "gray50"
	  ) +
	  theme_bw(base_size = 12) +
	  theme(
	    panel.grid.minor = element_blank(),
	    panel.grid.major = element_blank(),
	    plot.title = element_text(face = "bold"),
	    axis.title = element_text(face = "bold")
	  ) +
	  labs(
	    x = "Log2 fold-change (Tumor - Normal)",
	    y = expression(-log[10]("adjusted p-value")),
	    title = tmp_title
	  )

	p2 <- ggplot(plot_df, aes(x = logFC, y = neglog10_padj)) +
	  geom_point(
	    data = subset(plot_df, group == "NS"),
	    color = "gray80",
	    size = 1.2,
	    alpha = 0.7
	  ) +
	  geom_point(
	    data = subset(plot_df, group == "Up"),
	    color = "red",
	    size = 1.8,
	    alpha = 0.85
	  ) +
	  geom_point(
	    data = subset(plot_df, group == "Down"),
	    color = "blue",
	    size = 1.8,
	    alpha = 0.85
	  ) +
	  geom_vline(xintercept = c(-2, 2), linetype = "dashed", linewidth = 0.4) +
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
	  geom_text_repel(
	    data = label_df %>% filter(group != "NS"),
	    aes(label = gene_symbol),
	    size = 3,
	    max.overlaps = 50,
	    box.padding = 0.3,
	    point.padding = 0.2,
	    segment.color = "gray50"
	  ) +
	  theme_bw(base_size = 12) +
	  theme(
	    panel.grid.minor = element_blank(),
	    panel.grid.major = element_blank(),
	    plot.title = element_text(face = "bold"),
	    axis.title = element_text(face = "bold")
	  ) +
	  labs(
	    x = "Log2 fold-change (Tumor - Normal)",
	    y = expression(-log[10]("adjusted p-value")),
	    title = tmp_title
	  )

	# Paired box plots of the DE HS genes 
	genes <- res %>% filter(!is.na(gene_symbol) & group != "NS") %>% select(gene, gene_symbol)

  tumor_df <- expr_tumor[genes$gene, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "expr") %>%
    mutate(group = "tumor")

  normal_df <- expr_normal[genes$gene, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "expr") %>%
    mutate(group = "normal")

  combined <- rbind(tumor_df, normal_df)

  combined <- combined %>%
  	left_join(genes, by = "gene")

  p3 <- ggplot(combined, aes(x = group, y = expr, fill = group)) +
		  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
		  stat_summary(
		    fun = median,
		    geom = "point",
		    size = 2,
		    color = "black"
		  ) +
		  facet_wrap(~ gene_symbol, scales = "free_y") +
		  theme_bw(base_size = 12) +
		  theme(
		    strip.text = element_text(face = "bold"),
		    axis.title = element_text(face = "bold"),
		    legend.position = "none",
		    panel.grid.major.x = element_blank(),
		    panel.grid.minor = element_blank()
		  ) +
		  labs(
		    x = NULL,
		    y = "log(TPM)",
		    title = "Expression of selected genes"
		  )


  list(volc1 = p1, volc2 = p2, volc3 = p3)
}


# Loop through cancer types
	# First make a lookup table that maps normal <> tumor 
pheno[study == "GTEX",
      .N,
      by = .(detailed_site, primary_disease, primary_site, sample_type)
][order(primary_site)][N>5] %>% fwrite(., "data/gtex_unique_combinations.csv")

pheno[study == "TCGA",
      .N,
      by = .(detailed_site, primary_disease, primary_site, sample_type)
][order(primary_site)][N>5] %>% fwrite(., "data/tcga_unique_combinations.csv")

tcga_gtex_map <- fread("data/TCGA_GTEX_map.csv") # field is primary_disease 

parallel::mclapply(seq_len(nrow(tcga_gtex_map)), function(x) {
	print(x)

	cancer_tissue <- tcga_gtex_map$TCGA[x]
	healthy_tissue <- tcga_gtex_map$GTEX[x]

	nl <- pheno[study == "GTEX"][primary_disease == healthy_tissue]
	expr_normal <- expr[, colnames(expr) %in% nl$sample] %>% as.matrix()
	mode(expr_normal) <- "numeric"

	tm <- pheno[study == "TCGA"][primary_disease == cancer_tissue]
	expr_tumor <- expr[, colnames(expr) %in% tm$sample] %>% as.matrix()
	mode(expr_tumor) <- "numeric"

	p1 <- all_gene_expression(expr_normal, healthy_tissue)
	p2 <- all_gene_expression(expr_tumor, cancer_tissue)

	volcanos <- volcano(expr_tumor, expr_normal, paste0(cancer_tissue, " v. ", healthy_tissue))

	pdf(paste0("out/", cancer_tissue, ".pdf"), width = 10)
		print(p1 | p2)
		print(volcanos$volc1)
		print(volcanos$volc2)
		print(volcanos$volc3)
	dev.off()

}, mc.preschedule = FALSE)




