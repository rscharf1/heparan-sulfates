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
projection_path = 'inputs/bonemarrowmap/'

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






