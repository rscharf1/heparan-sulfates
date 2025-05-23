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

# EXAMPLE
query <- readRDS('inputs/ExampleQuery_Roy2021.rds')

query <- map_Query(
    exp_query = query[['RNA']]@counts,
    metadata_query = query@meta.data,
    ref_obj = ref,
    vars = 'sampleID'
)

# Run QC based on mapping error score, flag cells with mapping error >= 2.5 MADs above median
query <- query %>% calculate_MappingError(., reference = ref, MAD_threshold = 2.5) 

# Get QC Plots
QC_plots <- plot_MappingErrorQC(query)

pdf("outputs/sample_qc.pdf", width = 15)
patchwork::wrap_plots(QC_plots, ncol = 4, widths = c(0.8, 0.3, 0.8, 0.3))
dev.off()

query <- subset(query, mapping_error_QC == 'Pass')

query <- predict_CellTypes(
  query_obj = query, 
  ref_obj = ref, 
  initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 

pdf("outputs/sample_umap.pdf", width = 15)
DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType'), 
        raster=FALSE, label=TRUE, label.size = 4)
DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType_Broad'), 
        raster=FALSE, label=TRUE, label.size = 4)
dev.off()

# MY DATA
healthy <- readRDS("inputs/oetjen_2018/all_healthy_seurat.rds")
combined_counts <- do.call(cbind, healthy[["RNA"]]@layers)
rownames(combined_counts) <- rownames(healthy)
colnames(combined_counts) <- colnames(healthy)

my_query <- map_Query(
    exp_query = combined_counts, 
    metadata_query = healthy@meta.data,
    ref_obj = ref,
    vars = 'sample'
)

my_query <- my_query %>% calculate_MappingError(., reference = ref, MAD_threshold = 2.5) 

# Get QC Plots
QC_plots <- plot_MappingErrorQC(my_query)

pdf("outputs/healthy_qc.pdf", width = 15)
patchwork::wrap_plots(QC_plots, ncol = 4, widths = c(0.8, 0.3, 0.8, 0.3))
dev.off()

my_query <- subset(my_query, mapping_error_QC == 'Pass')

my_query <- predict_CellTypes(
  query_obj = my_query, 
  ref_obj = ref, 
  initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 

saveRDS(my_query, "inputs/healthy_mapped.rds")

healthy_mapped <- readRDS("inputs/healthy_mapped.rds")

pdf("outputs/healthy_umap.pdf", width = 15)
DimPlot(subset(healthy_mapped, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType'), 
        raster=FALSE, label=TRUE, label.size = 4)
DimPlot(subset(healthy_mapped, mapping_error_QC == 'Pass'), reduction = 'umap_projected', group.by = c('predicted_CellType_Broad'), 
        raster=FALSE, label=TRUE, label.size = 4)
dev.off()

pdf("outputs/dot.pdf")
DotPlot(healthy_mapped, features = toupper(enzymes), group.by = "predicted_CellType_Broad") + RotatedAxis()
dev.off()

pdf("outputs/test.pdf")
	FeaturePlot(healthy_mapped, features = "EXT2", reduction = "umap_projected")
dev.off()

hspc_labels <- c(
  "HSC MPP", "LMPP", "MEP", "Early GMP", "Cycling Progenitor",
  "Late GMP", "Early Erythroid", "Early Lymphoid", "EoBasoMast Precursor"
)

pdf("outputs/test2.pdf")
	FeaturePlot(subset(healthy_mapped, predicted_CellType_Broad %in% hspc_labels), features = "EXT2", reduction = "umap_projected")
dev.off()

# AML Data 






















