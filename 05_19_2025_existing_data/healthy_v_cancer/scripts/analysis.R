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
    DimPlot(healthy_mapped, reduction = 'umap_projected', group.by = c('predicted_CellType'), 
            raster=FALSE, label=TRUE, label.size = 4)
    DimPlot(healthy_mapped, reduction = 'umap_projected', group.by = c('predicted_CellType_Broad'), 
            raster=FALSE, label=TRUE, label.size = 4)
dev.off()

hspc_labels <- c(
  "HSC", "LMPP", "CLP", "MLP", "MLP-II",
  "MPP–MkEry", "MPP–MyLy", "Cycling Progenitor",
  "Early GMP", "GMP–Cycle", "GMP–Mono", "GMP–Neut",
  "Early ProMono", "Late ProMono",
  "MEP", "CFU–E", "BFU–E", "EoBasoMast Precursor",
  "Megakaryocyte Precursor", "Pro–Erythroblast"
)

healthy_hspc <- subset(healthy_mapped, subset = predicted_CellType %in% hspc_labels)

hsc_markers <- c("AVP", "HOPX", "MECOM", "HES1", "GATA2", "CD34", "PROCR", "MEG3", "MALAT1")

pdf("outputs/hsc.pdf", width = 12)
    DimPlot(healthy_hspc, group.by = "predicted_CellType", reduction = "umap_projected", label = TRUE)

    VlnPlot(healthy_hspc, features = hsc_markers,
            group.by = "predicted_CellType", pt.size = 0)

    DotPlot(healthy_hspc, features = hsc_markers,
            group.by = "predicted_CellType") + RotatedAxis()

    FeaturePlot(healthy_hspc, features = hsc_markers)
dev.off()


pdf("outputs/hspc_HS_enymes_dot.pdf", width = 12)
    DotPlot(healthy_hspc, features = toupper(enzymes), group.by = "predicted_CellType") + RotatedAxis()
    # FeaturePlot(healthy_hspc, features = "HS2ST1", reduction = "umap_projected")
dev.off()

# AML
aml <- readRDS("inputs/petti_2019/all_aml_seurat.rds")
aml <- merge(aml[[1]], y = aml[-1], project = "AML_BM")
combined_counts <- do.call(cbind, aml[["RNA"]]@layers)
rownames(combined_counts) <- rownames(aml)
colnames(combined_counts) <- colnames(aml)

my_query <- map_Query(
    exp_query = combined_counts, 
    metadata_query = aml@meta.data,
    ref_obj = ref,
    vars = 'sample_id'
)

my_query <- my_query %>% calculate_MappingError(., reference = ref, MAD_threshold = 2.5) 

# Get QC Plots
QC_plots <- plot_MappingErrorQC(my_query)

pdf("outputs/aml_qc.pdf", width = 15)
patchwork::wrap_plots(QC_plots, ncol = 4, widths = c(0.8, 0.3, 0.8, 0.3))
dev.off()

my_query <- subset(my_query, mapping_error_QC == 'Pass')

my_query <- predict_CellTypes(
  query_obj = my_query, 
  ref_obj = ref, 
  initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 

saveRDS(my_query, "inputs/aml_mapped.rds")

aml_mapped <- readRDS("inputs/aml_mapped.rds")

subset(aml_mapped, subset = predicted_CellType %in% hspc_labels)

pdf("outputs/hspc_HS_enymes_AML_dot.pdf", width = 12)
    DotPlot(
        subset(aml_mapped, subset = predicted_CellType %in% hspc_labels & mutation_label == "wildtype"), 
        features = toupper(enzymes), 
        group.by = "predicted_CellType") + 
    RotatedAxis() + 
    ggtitle("HSCs")

    DotPlot(
        subset(aml_mapped, subset = predicted_CellType %in% hspc_labels & mutation_label == "mutant"), 
        features = toupper(enzymes), 
        group.by = "predicted_CellType") + 
    RotatedAxis() + 
    ggtitle("LSCs")
dev.off()

sub <- subset(aml_mapped, subset = predicted_CellType == "HSC")

pdf("outputs/test.pdf")
VlnPlot(subset(aml_mapped, subset = predicted_CellType == "HSC"), 
    features = "HS2ST1", group.by = "mutation_label", pt.size = 0)
dev.off()


sub <- NormalizeData(sub)
sub <- ScaleData(sub)
table(FetchData(sub, "HS2ST1") > 0)

summary(FetchData(sub, "HS2ST1"))












