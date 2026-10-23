library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Plan
	# Read in TSS metadata 
	# For each cell line
		# Load in and combine signal tables 
		# Check correlation for treated v. untreated 
		# Then just focus on the vehicle 
			# For each peak, calculate mean and median 
			# Calculate enrichment over control
			# Rank order genes and then lot
				# Label pos/neg	ctrls and HS genes 

# Read in metadata
TSS_table <- fread("metadata/TSS_window_metadata.tsv")
TSS_table[, gene_symbol := sub("\\|.*$", "", gene_label)]
metadata <- fread("metadata/TSS_window_metadata.tsv")

# Read in cell line files 
CELL_LINE <- "HL60"

summary_files <- list.files(
  paste0("summaries/", CELL_LINE, "_summaries"),
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

signal_matrix[grepl("^EXT1\\|", gene_label)]

# Check correlation

treated_reps <- grep("Treated_Rep", colnames(signal_matrix), value = TRUE)
vehicle_reps <- grep("Veh_Rep", colnames(signal_matrix), value = TRUE)

cols <- colnames(signal_matrix)

treated_reps <- cols[grepl("Treated", cols) & grepl("Rep", cols)]
vehicle_reps <- cols[grepl("Veh", cols) & grepl("Rep", cols)]

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

all_reps <- c(treated_reps, vehicle_reps)

cor(
  as.matrix(signal_matrix[, ..all_reps]),
  method = "spearman",
  use = "pairwise.complete.obs"
)

# Just focus on vehicle 
vehicle_reps <- cols[grepl("Veh", cols) & grepl("Rep", cols)]
vehicle_ctrl <- cols[grepl("Veh", cols) & grepl("Ctrl", cols)]

summary <- signal_matrix[, .(
  region_id,
  ctrl_signal = get(vehicle_ctrl),
  mean_signal = rowMeans(.SD, na.rm = TRUE),
  median_signal = apply(.SD, 1, median, na.rm = TRUE)
), .SDcols = vehicle_reps]

pseudocount <- 0.01

summary[, `:=`(
  log2_enrichment_mean =
    log2((mean_signal + pseudocount) / (ctrl_signal + pseudocount)),
  log2_enrichment_median =
    log2((median_signal + pseudocount) / (ctrl_signal + pseudocount))
)]

summary[
  TSS_table[, .(region_id, gene_symbol)],
  gene_symbol := i.gene_symbol,
  on = "region_id"
]

summary <- summary[,
  .SD[which.max(log2_enrichment_median)],
  by = gene_symbol
]

# Get ready to plot 
hs_genes <- fread("~/steidl_lab/heparan-sulfates/useful_stuff/HS_genes_10182026.csv")
pos_ctrls <- c("CSF1R", "IRF8", "ITGAM", "CSF2RA", "FCGR2B")
neg_ctrls <- c("HBB", "EPCAM", "ALB", "MYOD1", "INS", "PAX6", "RHO", "KRT14")

data.table::setorder(summary, -log2_enrichment_mean)
summary[, rank := .I]

plot_dt <- copy(summary)

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


pdf(paste0(CELL_LINE, ".pdf"))

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














