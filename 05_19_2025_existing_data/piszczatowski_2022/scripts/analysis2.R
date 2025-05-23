library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)
library(stringr)
library(BoneMarrowMap)

enzymes <- c("Hs2st1", "Hs3st1", "Hs3st3b1", "Hs6st1", "Hs6st2", "Ndst3", 
"Ndst2", "Ndst1", "Glce", "Ext2", "Extl2", "Extl1", "Extl3", 
"Ext1", "Gpc4", "Tgfbr3", "Agrn")

# Set up map 
projection_path = '~/Tools/references/bonemarrowmap/'

ref <- readRDS(paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))
ref$save_uwot_path <- paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot')

ReferenceSeuratObj <- create_ReferenceObject(ref)

ReferenceSeuratObj$CellType_Annotation_formatted <- ReferenceSeuratObj$CellType_Annotation_formatted %>% 
	gsub("\t\n|\t|\n", " ", .)

pdf("outputs/ref.pdf", width = 25, height = 15)
	DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CellType_Annotation_formatted', 
	        raster=FALSE, label=TRUE, label.size = 4) + NoAxes()
	DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CellType_Broad', 
	        raster=FALSE, label=TRUE, label.size = 4)
	FeaturePlot(ReferenceSeuratObj, reduction = 'umap', features = 'Pseudotime', raster=FALSE) + NoAxes()
dev.off()

# Prepare my data 
prep_mat <- function(path) {
	dt <- fread(path)

	mat_dense <- as.matrix(dt[, -1, with = FALSE])
	rownames(mat_dense) <- dt[[1]]

	rownames(mat_dense) <- toupper(rownames(mat_dense))

	seurat_obj <- CreateSeuratObject(counts = mat_dense)

	seurat_obj
}

high <- prep_mat("inputs/GSM6260301_hs3a8_high_counts.mat")
low <- prep_mat("inputs/GSM6260302_hs3a8_low_counts.mat")

high$group <- "high"
low$group <- "low"

seu_combined <- merge(high, y = low)

rna_counts_combined <- cbind(
	GetAssayData(seu_combined[["RNA"]], layer = "counts.1"),
	GetAssayData(seu_combined[["RNA"]], layer = "counts.2")
)

# rna_counts_combined <- as(rna_counts_combined, "dgCMatrix")

# seu_combined <- NormalizeData(seu_combined)
# seu_combined <- FindVariableFeatures(seu_combined)
# seu_combined <- ScaleData(seu_combined)
# seu_combined <- RunPCA(seu_combined, npcs = 30)

query <- map_Query(
    exp_query = rna_counts_combined, #query[['RNA']]@counts, 
    metadata_query = seu_combined@meta.data,
    ref_obj = ref,
    vars = 'group'
)

query <- query %>% calculate_MappingError(., reference = ref, MAD_threshold = 2.5) 

# Get QC Plots
QC_plots <- plot_MappingErrorQC(query)

pdf("outputs/qc.pdf", width = 15)
patchwork::wrap_plots(QC_plots, ncol = 4, widths = c(0.8, 0.3, 0.8, 0.3))
dev.off()

query <- subset(query, mapping_error_QC == 'Pass')

query <- predict_CellTypes(
  query_obj = query, 
  ref_obj = ref, 
  initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 

pdf("outputs/proj_umap.pdf", width = 15)
	DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType'), 
	        raster=FALSE, label=TRUE, label.size = 4)
	DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType_Broad'), 
	        raster=FALSE, label=TRUE, label.size = 4)
	DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType_Broad'), 
	        raster=FALSE)
	DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('group'), 
        raster=FALSE)
dev.off()

dt <- query@meta.data %>% data.table()

dt[grepl("Erythroid", dt$predicted_CellType_Broad)]$group %>% table

dt[grepl("MEP|Megakaryocyte", dt$predicted_CellType_Broad)]$group %>% table

pdf("outputs/v.pdf", width = 15)
VlnPlot(
  object = query, 
  features = toupper(enzymes[1]), 
  group.by = "predicted_CellType_Broad", 
  pt.size = 0
)
dev.off()

pdf("outputs/dot.pdf", width = 12)
	DotPlot(query, features = toupper(enzymes), group.by = "predicted_CellType_Broad") + RotatedAxis()
dev.off()


