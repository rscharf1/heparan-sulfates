library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)
library(stringr)
library(tidyr)
library(parallel)

files <- list.files(path = "inputs/petti_2019", pattern = "\\.mtx.gz$", full.names = TRUE)
mut_files <- list.files(path = "inputs/petti_2019/scrna_mutations/cb_sniffer_results", full.names = TRUE)
aml <- parallel::mclapply(seq_along(files[1:5]), function(f) {
	id <- str_extract(files[f], "\\d{6}")

	barcodes <- files[f] %>% gsub("matrix", "barcodes", .) %>% gsub("mtx", "tsv", .)
	genes <- files[f] %>% gsub("matrix", "genes", .) %>% gsub("mtx", "tsv", .)
	matrix <- files[f]

	mat <- ReadMtx(
		mtx = matrix, 
		features = genes, 
		cells = barcodes
	)

	seu <- CreateSeuratObject(counts = mat)

	seu@meta.data$cell_barcode <- rownames(seu@meta.data)

	mut <- fread(mut_files[grepl(id, mut_files)][1])

	mut_dt <- mut[, .(mutant = any(AltReads > 0)), by = Cell][, mutation_label := ifelse(mutant, "mutant", "wildtype")]

	a <- left_join(
	  seu@meta.data,
	  mut_dt,
	  by = c("cell_barcode" = "Cell")
	)

	a <- a %>% replace_na(list(mutation_label = "Unknown"))

	seu@meta.data$mutation_label <- a$mutation_label

	seu@meta.data$sample_id <- id

	seu@meta.data$disease_status <- "AML"

	seu
}, mc.cores=8)

files <- list.files(path = "inputs/oetjen_2018", pattern = "\\.mtx.gz$", full.names = TRUE)

healthy <- parallel::mclapply(seq_along(files), function(f) {
	message(paste0("File ", as.character(f), " / ", length(files)))

	sample_id <- str_extract(files[f], "(?<=matrix_).+(?=\\.mtx\\.gz)")

	barcodes <- files[f] %>% gsub("matrix", "barcodes", .) %>% gsub("mtx", "tsv", .)
	genes <- files[f] %>% gsub("matrix", "genes", .) %>% gsub("mtx", "tsv", .)
	matrix <- files[f]

	if(file.exists(barcodes) & file.exists(genes) & file.exists(matrix)) {

		mat <- ReadMtx(
			mtx = matrix, 
			features = genes, 
			cells = barcodes
		)

		seu <- CreateSeuratObject(counts = mat, project = paste0("Sample_", sample_id))
		colnames(seu) <- paste0(sample_id, "_", colnames(seu))
		seu$sample <- sample_id  # Add sample metadata
		seu$disease_status <- "healthy"

		seu
	} 
	else (
		message("No bueno")
	)
}, mc.cores=8)

all_samples <- c(healthy, aml)

all_samples_norm <- parallel::mclapply(all_samples, function(seu) {
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = 30)
  seu
}, mc.cores=8) %>% unlist()

anchors <- FindIntegrationAnchors(
	object.list = all_samples_norm, 
	normalization.method = "LogNormalize",
	reduction = "rpca",
	dims = 1:20
)

seu_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

saveRDS(seu_integrated, "outputs/seu_integrated.rds")






