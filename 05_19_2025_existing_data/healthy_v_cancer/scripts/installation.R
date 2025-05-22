library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)
library(stringr)

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

projection_path = 'inputs/bonemarrowmap/'

# Download Bone Marrow Reference - 344 Mb
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_SymphonyReference.rds', 
                    destfile = paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))
# Download uwot model file - 221 Mb
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_uwot_model.uwot', 
                    destfile = paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot'))
