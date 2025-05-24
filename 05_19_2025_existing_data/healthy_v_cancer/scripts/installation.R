library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)
library(stringr)
library(tidyr)

################################
# LOAD AND COMBINE HEALTHY BONE MARROW scRNA DATA 
files <- list.files(path = "inputs/oetjen_2018", pattern = "\\.mtx.gz$", full.names = TRUE)

seurat_list <- lapply(seq_along(files), function(f) {
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

		seu
	} 
	else (
		message("No bueno")
	)

}) %>% unlist()

healthy <- merge(seurat_list[[1]], y = seurat_list[-1], project = "HealthyBM")

saveRDS(healthy, "inputs/oetjen_2018/all_healthy_seurat.rds")

################################
# BONE MARROW MAP
my_lib <- "/gs/gsfs0/users/rscharf/R/x86_64-pc-linux-gnu-library/4.4"
.libPaths(my_lib)

Sys.setenv(CXX17 = "g++ -std=c++17")
Sys.setenv(CXX17FLAGS = "-O2 -fPIC")

if(!require(BiocManager, quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AUCell", "doMC", "BiocNeighbors"))
if(!require(devtools, quietly = TRUE)) install.packages("devtools")
devtools::install_github("jaredhuling/jcolors") # dependency that is no longer on CRAN

devtools::install_github('andygxzeng/BoneMarrowMap', force = TRUE)

projection_path = '~/Tools/references/bonemarrowmap/'

# Download Bone Marrow Reference - 344 Mb
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_SymphonyReference.rds', 
                    destfile = paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))
# Download uwot model file - 221 Mb
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_uwot_model.uwot', 
                    destfile = paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot'))

################################
# AML DATA 

# barcodes <- "inputs/petti_2019/508084.barcodes.tsv.gz"
# genes <- "inputs/petti_2019/508084.genes.tsv.gz"
# matrix <- "inputs/petti_2019/508084.matrix.mtx.gz"

# mat <- ReadMtx(
# 	mtx = matrix, 
# 	features = genes, 
# 	cells = barcodes
# )

# # seu <- CreateSeuratObject(counts = mat, project = paste0("Sample_", sample_id))

# toupper(enzymes) %in% rownames(mat)

# seu <- readRDS("inputs/petti_2019/548327.seurat.rds")


# seu <- UpdateSeuratObject(seu)
# seu@meta.data %>% dim()
# rownames(seu@meta.data) %>% length()
# rownames(seu@meta.data) %>% unique() %>% length()

# mut <- fread("inputs/petti_2019/scrna_mutations/cb_sniffer_results/548327_CellCoverage_Coding_180813.txt")

# mut %>%
#   group_by(Cell) %>%
#   summarise(mutant = any(AltReads > 0)) %>%
#   mutate(mutation_label = ifelse(mutant, "mutant", "wildtype"))

# mut_dt <- mut[, .(mutant = any(AltReads > 0)), by = Cell][, mutation_label := ifelse(mutant, "mutant", "wildtype")]

# seu@meta.data$cell_barcode <- rownames(seu@meta.data)

# seu@meta.data <- left_join(
#   seu@meta.data,
#   mut_dt,
#   by = c("cell_barcode" = "Cell")
# )

####################
# AML Data 
# Using cellranger output to assemble seurat object 
mut <- fread("inputs/petti_2019/scrna_mutations/cb_sniffer_results/548327_CellCoverage_Coding_180813.txt")

mut_dt <- mut[, .(mutant = any(AltReads > 0)), by = Cell][, mutation_label := ifelse(mutant, "mutant", "wildtype")]

# MY SEURAT 
barcodes <- "inputs/petti_2019/548327.barcodes.tsv.gz"
genes <- "inputs/petti_2019/548327.genes.tsv.gz"
matrix <- "inputs/petti_2019/548327.matrix.mtx.gz"

mat <- ReadMtx(
	mtx = matrix, 
	features = genes, 
	cells = barcodes
)

seu <- CreateSeuratObject(counts = mat)

seu@meta.data$cell_barcode <- rownames(seu@meta.data)

# seu@meta.data <- left_join(
#   seu@meta.data,
#   mut_dt,
#   by = c("cell_barcode" = "Cell")
# )

a <- left_join(
  seu@meta.data,
  mut_dt,
  by = c("cell_barcode" = "Cell")
)

a <- a %>% replace_na(list(mutation_label = "Unknown"))

seu@meta.data$mutation_label <- a$mutation_label

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:10)

seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)

pdf("outputs/mine.pdf")
	DimPlot(seu, reduction = "umap", pt.size = 0.5)
	DimPlot(seu, reduction = 'umap', group.by = "mutation_label", pt.size = 0.5)
dev.off()

# Loop it -> create AML.rds 
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

aml <- lapply(seq_along(files[1:5]), function(f) {
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
})

saveRDS(aml, "inputs/petti_2019/all_AML_seurat.rds")

# files <- list.files(path = "inputs/oetjen_2018", pattern = "\\.mtx.gz$", full.names = TRUE)

# healthy <- parallel::mclapply(seq_along(files), function(f) {
# 	message(paste0("File ", as.character(f), " / ", length(files)))

# 	sample_id <- str_extract(files[f], "(?<=matrix_).+(?=\\.mtx\\.gz)")

# 	barcodes <- files[f] %>% gsub("matrix", "barcodes", .) %>% gsub("mtx", "tsv", .)
# 	genes <- files[f] %>% gsub("matrix", "genes", .) %>% gsub("mtx", "tsv", .)
# 	matrix <- files[f]

# 	if(file.exists(barcodes) & file.exists(genes) & file.exists(matrix)) {

# 		mat <- ReadMtx(
# 			mtx = matrix, 
# 			features = genes, 
# 			cells = barcodes
# 		)

# 		seu <- CreateSeuratObject(counts = mat, project = paste0("Sample_", sample_id))
# 		colnames(seu) <- paste0(sample_id, "_", colnames(seu))
# 		seu$sample <- sample_id  # Add sample metadata
# 		seu$disease_status <- "healthy"

# 		seu
# 	} 
# 	else (
# 		message("No bueno")
# 	)
# }, mc.cores=8) %>% unlist()

# all_samples <- c(healthy, aml)

# all_samples_norm <- parallel::mclapply(all_samples, function(seu) {
#   seu <- NormalizeData(seu)
#   seu <- FindVariableFeatures(seu)
#   seu <- ScaleData(seu)
#   seu <- RunPCA(seu, npcs = 30)
#   seu
# }, mc.cores=8) %>% unlist()

# anchors <- FindIntegrationAnchors(
# 	object.list = all_samples_norm, 
# 	normalization.method = "LogNormalize",
# 	reduction = "rpca",
# 	dims = 1:20
# )

# seu_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)



# # THEIR SEURAT
# seu2 <- readRDS("inputs/petti_2019/548327.seurat.rds")
# seu2 <- UpdateSeuratObject(seu2)
# saveRDS(seu2, "inputs/petti_2019/548327.seurat.updated.rds")

# seu2 <- readRDS("inputs/petti_2019/548327.seurat.updated.rds")

# intersect(rownames(seu2@meta.data), rownames(seu@meta.data)) %>% length()

# seu2@meta.data$cell_barcode <- rownames(seu2@meta.data)

# a <- left_join(
#   seu2@meta.data,
#   mut_dt,
#   by = c("cell_barcode" = "Cell")
# )

# a <- a %>% replace_na(list(mutation_label = "Unknown"))

# seu2@meta.data$mutation_label <- a$mutation_label

# seu2 <- RunUMAP(seu2, reduction = "pca", dims = 1:10)

# pdf("theirs.pdf")
# 	DimPlot(seu2, reduction = 'umap', group.by = "CellType", label = TRUE, pt.size = 0.5)
# 	DimPlot(seu2, reduction = 'umap', group.by = "mutation_label", pt.size = 0.5)
# dev.off()



out <- parallel::mclapply(files, function(f) {
	message("Ripping it")
	dt <- readRDS(f)
	dt[, id:=paste(aa_cdr1, aa_cdr2, aa_cdr3, sep="_")]
	dt <- dt[, c("id", "N"), with=FALSE]

	dt[, abundance:=N/sum(N)]

	dt <- hamming_collapse(dt)

	dt$sample <- basename(tools::file_path_sans_ext(f))

	dt

}, mc.cores = 5) %>% Reduce("rbind", .) 



healthy <- parallel::mclapply(seq_along(files[1:2]), function(f) {
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
}, mc.cores=8) %>% unlist()
