library(data.table)
library(dplyr)
library(ggplot2)

dt <- fread("inputs/clean_counts.txt")

# Clean up gene table 
gene_table <- fread("inputs/annotation_table_clean.txt") %>% .[, c("ID", "Gene Symbol"), with = FALSE]

colnames(gene_table) <- c("id", "gene1")

gene_table[, gene2 :=
  trimws(sub("^[^/]*///", "", gene1))
]

gene_table[, gene1 :=
  trimws(sub("[[:space:]]*///.*$", "", gene1))
]

gene_table[gene2 == gene1, gene2 := NA_character_]
gene_table[gene2 == "", gene2 := NA_character_]

gene_table <- gene_table[nchar(gene1) > 0]

gene_table[, gene1_probe := {
  if (.N > 1L && !is.na(gene1[1L]) && nzchar(gene1[1L])) {
    paste0(gene1, "_probe", seq_len(.N))
  } else {
    gene1
  }
}, by = gene1]

# Join expression and gene data
annot <- gene_table[, .(ID_REF = id, gene = gene1_probe)]

expr_annotated <- merge(dt, annot, by = "ID_REF")

sample_cols <- setdiff(names(dt), "ID_REF")
expr_matrix <- as.matrix(expr_annotated[, ..sample_cols])
rownames(expr_matrix) <- expr_annotated$gene

colnames(expr_matrix) <- c(
  paste0("MDS_", seq_len(183)),
  paste0("Healthy_", seq_len(17))
)

# Transformation Data
transformation <- fread("inputs/transformation_data.csv") %>% .[, c(1,4)]
colnames(transformation) <- c("patient", "aml")
transformation[, patient_original := patient]
transformation[, patient := paste0("MDS_", seq_len(.N))]

transformation$aml %>% table()

# Correlative testing 
expr_matrix[grep("ext1", rownames(expr_matrix), ignore.case = TRUE),, drop = FALSE]

gene_row <- "CTSL"

stopifnot(gene_row %in% rownames(expr_matrix))

expr_dt <- data.table(
  patient = colnames(expr_matrix),
  expression = as.numeric(expr_matrix[gene_row, ])
)

outcome <- transformation[, .(patient, aml)][aml %chin% c("YES", "NO")]

analysis_dt <- merge(expr_dt, outcome, by = "patient")
analysis_dt[, progressed := as.integer(aml == "YES")]

# Group sizes and expression summaries
analysis_dt[, .(
  n = .N,
  median_expression = median(expression, na.rm = TRUE),
  mean_expression = mean(expression, na.rm = TRUE)
), by = aml]

# Quick comparison between the two groups
wilcox.test(expression ~ aml, data = analysis_dt, exact = FALSE)

# Odds of progression per 1 SD increase in expression
analysis_dt[, expression_z := as.numeric(scale(expression))]
fit <- glm(progressed ~ expression_z, data = analysis_dt, family = binomial())

summary(fit)
exp(cbind(OR = coef(fit), confint.default(fit)))

test_one_gene <- function(gene_row) {
  expr_dt <- data.table(
    patient = colnames(expr_matrix),
    expression = as.numeric(expr_matrix[gene_row, ])
  )

  outcome <- copy(transformation[, .(patient, aml)])
  outcome[, aml := toupper(trimws(as.character(aml)))]
  outcome <- outcome[aml %chin% c("YES", "NO")]

  healthy <- data.table(
    patient = grep("^Healthy_[0-9]+$", colnames(expr_matrix), value = TRUE),
    aml = "Healthy"
  )

  outcome <- rbindlist(list(outcome, healthy), use.names = TRUE)

  analysis_dt <- merge(expr_dt, outcome, by = "patient")

  analysis_dt[, .(
    probe = gene_row,
    n = .N,
    median_expression = median(expression, na.rm = TRUE),
    mean_expression = mean(expression, na.rm = TRUE)
  ), by = aml]
}

out <- lapply(matched_rows, function(g) {
	test_one_gene(g)
}) %>% rbindlist()

test_one_gene("CTSL")






list_of_genes <- c("ext1", "extl3", "hpse")
row_names <- rownames(expr_matrix)
rows_lower <- tolower(row_names)

matched_rows <- unique(unlist(lapply(list_of_genes, function(g) {
  row_names[rows_lower == g | startsWith(rows_lower, paste0(g, "_probe"))]
}), use.names = FALSE))

hs_genes <- c(
	fread("~/steidl_lab/heparan-sulfates/useful_stuff/HS_relevant_genes.csv") %>% 
	.$hgnc_symbol,
	fread("~/steidl_lab/heparan-sulfates/useful_stuff/HS_degradation.csv", header=FALSE) %>% 
		.$V1
) %>% tolower() %>% unique()

matched_rows <- unique(unlist(lapply(hs_genes, function(g) {
  row_names[rows_lower == g | startsWith(rows_lower, paste0(g, "_probe"))]
}), use.names = FALSE))

out <- lapply(matched_rows, function(g) {
	test_one_gene(g)
}) %>% rbindlist()

plot_dt <- copy(out)

wide <- dcast(plot_dt, probe ~ aml, value.var = "median_expression")
wide[, log2_ratio := log2(YES / NO)]
setorder(wide, log2_ratio)

# plot_dt[, probe := factor(probe, levels = wide$probe)]
plot_dt[, probe := factor(probe, levels = sort(unique(probe)))]

pdf("10162026/out1.pdf", height = 12)

ggplot(plot_dt, aes(
  x = median_expression,
  y = probe,
  color = aml,
  group = probe
)) +
  geom_line(color = "grey75") +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    x = "Median expression (log10 display scale)",
    y = NULL,
    color = "Progressed to AML"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6))

dev.off()

plot_dt <- dcast(out, probe ~ aml, value.var = "median_expression")

pdf("10162026/out2.pdf")

ggplot(plot_dt, aes(x = YES, y = NO)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point() +
  coord_equal() +
  labs(
    x = "Median expression in MDS patients who progressed to AML",
    y = "Median expression in MDS patients who did not progress"
  ) +
  theme_bw()

dev.off()



