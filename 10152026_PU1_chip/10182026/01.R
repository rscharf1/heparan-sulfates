library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

metadata <- fread("metadata/TSS_window_metadata.tsv")

summary_files <- list.files(
  "MOLM13_summaries",
  pattern = "\\.tsv$",
  full.names = TRUE
)

signal_tables <- lapply(summary_files, function(file) {
  x <- fread(
    file,
    header = FALSE,
    col.names = c("region_id", "size", "covered", "sum", "mean0", "mean")
  )

  sample_name <- sub("\\.tsv$", "", basename(file))
  result <- x[, .(region_id, mean0)]
  setnames(result, "mean0", sample_name)
  result
})

signal_matrix <- Reduce(
  function(x, y) merge(x, y, by = "region_id"),
  signal_tables
)

signal_matrix <- merge(metadata, signal_matrix, by = "region_id")

stopifnot(nrow(signal_matrix) == 103812)

treated_reps <- c(
  "MOLM13_Treated_12hr_Rep1",
  "MOLM13_Treated_12hr_Rep2",
  "MOLM13_Treated_12hr_Rep3"
)

vehicle_reps <- c(
  "MOLM13_Veh_12hr_Rep1",
  "MOLM13_Veh_12hr_Rep2",
  "MOLM13_Veh_12hr_Rep3"
)

cor(
  as.matrix(signal_matrix[, ..treated_reps]),
  method = "spearman",
  use = "pairwise.complete.obs"
)

cor(
  as.matrix(signal_matrix[, ..vehicle_reps]),
  method = "spearman",
  use = "pairwise.complete.obs"
)

TSS_table <- fread("metadata/TSS_window_metadata.tsv")
TSS_table[grepl("EXT1", gene_label)]$region_id

signal_matrix[region_id %chin% TSS_table[grepl("EXT1", gene_label)]$region_id]

# What am I really asking here? 
	# We want to know which HS genes have convincing PU.1 binding over the control 
	# And among these, where does EXT1 fall? 
	# The drug treatment is secondary, and can be used after answering these questions above with the vehicle controls 

# Lets start with looking at all genes 
	# Some number that summarizes the 3 reps over the control
	# See how many have "convincing PU.1 binding"

TSS_table <- fread("metadata/TSS_window_metadata.tsv")
TSS_table[, gene_symbol := sub("\\|.*$", "", gene_label)]

signal_matrix

veh_rep_cols <- paste0("MOLM13_Veh_12hr_Rep", 1:3)

veh_summary <- signal_matrix[, .(
  region_id,
  ctrl_signal = MOLM13_Veh_12hr_Ctrl,
  mean_signal = rowMeans(.SD, na.rm = TRUE),
  median_signal = apply(.SD, 1, median, na.rm = TRUE)
), .SDcols = veh_rep_cols]

veh_summary <- veh_summary[!region_id %chin% TSS_table[chr == "chrM"]$region_id]

veh_summary[mean_signal > 0][order(mean_signal)]
veh_summary[median_signal == 0] %>% nrow()

pseudocount <- 0.01

veh_summary[, `:=`(
  log2_enrichment_mean =
    log2((mean_signal + pseudocount) / (ctrl_signal + pseudocount)),
  log2_enrichment_median =
    log2((median_signal + pseudocount) / (ctrl_signal + pseudocount))
)]

pdf("out1.pdf")

ggplot(veh_summary, aes(x = log2_enrichment_mean, y = log2_enrichment_median)) + 
	geom_point() +
	geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") + 
	coord_equal()

dev.off()

veh_summary[
  TSS_table[, .(region_id, gene_symbol)],
  gene_symbol := i.gene_symbol,
  on = "region_id"
]

veh_summary <- veh_summary[,
  .SD[which.max(log2_enrichment_median)],
  by = gene_symbol
]


# On the above plot, show where the HS genes, pos ctrls, and neg ctrls fall 
hs_genes <- fread("~/steidl_lab/heparan-sulfates/useful_stuff/HS_genes_10182026.csv")
pos_ctrls <- c("CSF1R", "IRF8", "ITGAM", "CSF2RA", "FCGR2B")
neg_ctrls <- c("HBB", "EPCAM", "ALB", "MYOD1", "INS", "PAX6", "RHO", "KRT14")

label_dt <- veh_summary[gene_symbol %chin% hs_genes$Gene]

pdf("out2.pdf")

	ggplot(
	  veh_summary,
	  aes(x = log2_enrichment_mean, y = log2_enrichment_median)
	) +
	  geom_point(color = "grey40", alpha = 0.2, size = 0.6) +
	  geom_abline(
	    slope = 1, intercept = 0,
	    color = "red", linetype = "dashed"
	  ) +
	  geom_point(
	    data = label_dt,
	    color = "steelblue",
	    size = 1.5
	  ) +
	  geom_text_repel(
	    data = label_dt,
	    aes(label = gene_symbol),
	    size = 3,
	    max.overlaps = Inf,
	    show.legend = FALSE
	  ) +
	  coord_equal() +
	  labs(
	    x = "Mean log2 enrichment",
	    y = "Median log2 enrichment"
	  ) +
	  theme_bw()

dev.off()

plot_dt <- copy(veh_summary)

plot_dt[, plot_group := fcase(
  gene_symbol %chin% pos_ctrls, "Positive controls",
  gene_symbol %chin% neg_ctrls, "Negative controls",
  gene_symbol %chin% hs_genes[Category %in% c("Synthesis_1", "Synthesis_2")]$Gene, "HS pathway genes",
  default = "Other"
)]

genes_to_label <- unique(c(hs_genes[Category %in% c("Synthesis_1", "Synthesis_2")]$Gene, pos_ctrls, neg_ctrls))

label_dt <- plot_dt[gene_symbol %chin% genes_to_label]

pdf("out3.pdf")

	ggplot(
	  plot_dt,
	  aes(x = log2_enrichment_mean, y = mean_signal)
	) +
	  geom_point(color = "grey40", alpha = 0.2, size = 0.6) +
	  geom_point(
	    data = plot_dt[plot_group != "Other"],
	    aes(color = plot_group),
	    size = 1.5
	  ) +
	  geom_text_repel(
	    data = label_dt,
	    aes(label = gene_symbol, color = plot_group),
	    size = 3,
	    max.overlaps = Inf,
	    show.legend = FALSE
	  ) +
	  scale_color_manual(values = c(
	    "HS pathway genes" = "firebrick3",
	    "Positive controls" = "darkorange",
	    "Negative controls" = "purple"
	  )) +
	  # coord_equal() +
	  labs(
	    x = "log2 enrichment",
	    y = "Mean Signal",
	    color = NULL
	  ) +
	  theme_bw()

dev.off()

# Sort by enrichment -> rank -> x: rank, y: enrichment 
data.table::setorder(veh_summary, -log2_enrichment_mean)
veh_summary[, rank := .I]

plot_dt <- copy(veh_summary)

plot_dt[, plot_group := fcase(
  gene_symbol %chin% pos_ctrls, "Positive controls",
  gene_symbol %chin% neg_ctrls, "Negative controls",
  gene_symbol %chin%
    hs_genes[Category %in% c("Synthesis_1", "Synthesis_2")]$Gene,
  "HS pathway genes",
  default = "Other"
)]

genes_to_label <- unique(c(
  hs_genes[Category %in% c("Synthesis_1", "Synthesis_2")]$Gene,
  pos_ctrls,
  neg_ctrls
))

label_dt <- plot_dt[gene_symbol %chin% genes_to_label]

pdf("out4.pdf")

ggplot(
  plot_dt,
  aes(x = rank, y = log2_enrichment_mean)
) +
  geom_point(color = "grey60", alpha = 0.5, size = 0.7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_point(
    data = plot_dt[plot_group != "Other"],
    aes(color = plot_group),
    size = 1.5
  ) +
  geom_text_repel(
    data = label_dt,
    aes(label = gene_symbol, color = plot_group),
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "HS pathway genes" = "steelblue",
    "Positive controls" = "darkorange",
    "Negative controls" = "purple"
  )) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "Rank by mean log2 enrichment",
    y = "Mean log2 enrichment",
    color = NULL
  ) +
  theme_bw()

dev.off()



