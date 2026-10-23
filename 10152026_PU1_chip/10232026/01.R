library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Read in metadata
TSS_table <- fread("metadata/TSS_window_metadata.tsv")
TSS_table[, gene_symbol := sub("\\|.*$", "", gene_label)]
metadata <- fread("metadata/TSS_window_metadata.tsv")

cell_lines <- c("HL60", "MOLM13", "MV411", "THP1")

cell_lines <- c("HL60", "MOLM13")

dt <- lapply(cell_lines, function(CELL_LINE) {
  print(paste0("Running for: ", CELL_LINE))

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

# Just focus on vehicle 
cols <- colnames(signal_matrix)

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

summary[, cell_line := CELL_LINE]


}) %>% rbindlist()

dt[gene_symbol == "EXT1"]


hs_genes <- fread("~/steidl_lab/heparan-sulfates/useful_stuff/HS_genes_10182026.csv")
pos_ctrls <- c("CSF1R", "IRF8", "ITGAM", "CSF2RA", "FCGR2B")
neg_ctrls <- c("HBB", "EPCAM", "ALB", "MYOD1", "INS", "PAX6", "RHO", "KRT14")

to_plot <- copy(dt)

to_plot <- to_plot[gene_symbol %in% hs_genes[Category %in% c("Synthesis_1", "Synthesis_2")]$Gene]



pdf("out1.pdf")

ggplot(to_plot, aes(
  x = gene_symbol,
  y = log2_enrichment_mean,
  color = cell_line
)) +
  geom_point(
    size = 3,
    position = position_dodge(width = 0.5)
  ) +
  labs(
    x = "Gene",
    y = "Log2 enrichment (mean)",
    color = "Cell line"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

dev.off()

plot_data <- to_plot %>%
  group_by(gene_symbol) %>%
  mutate(
    deviation = log2_enrichment_mean -
      median(log2_enrichment_mean, na.rm = TRUE)
  ) %>%
  ungroup()

pdf("out2.pdf", width = 10, height = 5)

ggplot(plot_data, aes(
  x = gene_symbol,
  y = deviation,
  color = cell_line,
  shape = cell_line
)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(
    size = 3,
    position = position_dodge(width = 0.5)
  ) +
  labs(
    x = "Gene",
    y = "Enrichment - gene median",
    color = "Cell line",
    shape = "Cell line"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

dev.off()

# Need to filter out ultra lowly expressed ones 

to_plot_filtered <- to_plot %>%
  group_by(gene_symbol) %>%
  filter(any(mean_signal >= 0.1)) %>%
  ungroup()

plot_data <- to_plot_filtered %>%
  group_by(gene_symbol) %>%
  mutate(
    deviation = log2_enrichment_mean -
      median(log2_enrichment_mean, na.rm = TRUE)
  ) %>%
  ungroup()

pdf("out3.pdf", width = 10, height = 5)

ggplot(plot_data, aes(
  x = gene_symbol,
  y = deviation,
  color = cell_line,
  shape = cell_line
)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(
    size = 3,
    position = position_dodge(width = 0.5)
  ) +
  labs(
    x = "Gene",
    y = "Enrichment - gene median",
    color = "Cell line",
    shape = "Cell line"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

dev.off()

