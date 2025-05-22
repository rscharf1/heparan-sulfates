library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)

# Enzymes
	# Hs2st1, Hs3st1-6, Hs6st1-3
	# Ndst1-4
	# Glce
		# Ext1-2

# try reading in the unzipped version

prep_mat <- function(path) {
	dt <- fread(path)

	mat_dense <- as.matrix(dt[, -1, with = FALSE])
	rownames(mat_dense) <- dt[[1]]

	seurat_obj <- CreateSeuratObject(counts = mat_dense)

	seurat_obj
}

genes <- c("Hs2st1", "Hs3st1", "Hs3st3b1", "Hs6st1", "Hs6st2", "Ndst3", 
"Ndst2", "Ndst1", "Glce", "Ext2", "Extl2", "Extl1", "Extl3", 
"Ext1", "Gpc4", "Tgfbr3", "Agrn")

high <- prep_mat("inputs/GSM6260301_hs3a8_high_counts.mat")
low <- prep_mat("inputs/GSM6260302_hs3a8_low_counts.mat")

high$group <- "high"
low$group <- "low"

combined <- merge(high, y = low)

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)

combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:10)

combined <- FindNeighbors(combined, dims = 1:10)

combined <- FindClusters(combined, resolution = 0.5)

pdf("outputs/umap.pdf")
DimPlot(combined, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.5) + 
	ggtitle("All Cells")
DimPlot(subset(combined, subset = group == "high"), reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.5) + 
	ggtitle("HS3A8 High Cells")
DimPlot(subset(combined, subset = group == "low"), reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.5) + 
	ggtitle("HS3A8 Low Cells")

# subset(combined, subset = group == "high")
dev.off()

#############
markers <- read.csv(text = "
	Cluster,Cell_Type,Marker1,Marker2,Marker3
	1,CEP,Hba-a1,Hbb-bs,Klf1
	2,EEP,Epor,Klf1,Tfrc
	3,MEP,Gata1,Zfpm1,Mpl
	4,Mk,Itga2b,Pf4,Gp9
	5,Baso/Mast,Cpa3,Ms4a2,Fcer1a
	6,Mono,Csf1r,Ly6c2,Fcgr3
	7,Gran,Elane,Mpo,S100a9
	8,Dendritic,Flt3,Cd209a,Zbtb46
	9,Ly-T/NK,Cd3d,Nkg7,Gzmb
	10,Ly-B,Ms4a1,Cd79a,Bcl11a
	11,GMP,Mpo,Elane,Cebpe
	12,CLP,Dntt,Rag1,Il7r
	13,CMP,Cebpa,Csf1r,Runx1
	14,MPP,Mecom,Hlf,Procr"
)

genes_to_plot <- c("Kit", "Gata1", "Gfi1b", "Mpl", "Itga2b", "Flt3", "Klf1", "Tfrc", "Hbb-bs", "Hba-a1", "Ermap", "Cd36", "Gypa", "Epor", "Slc4a1", "Alas2")

pdf("outputs/dot.pdf")
DotPlot(combined, features = genes_to_plot, group.by = "seurat_clusters") + RotatedAxis()
dev.off()

new_labels <- c(
  "3" = "MEP",
  "5" = "CEP",
  "8" = "EEP",
  "9" = "EEP"
)

vec <- combined$seurat_clusters %>% as.vector()

relabelled_vector <- ifelse(vec %in% names(new_labels),
                            new_labels[vec],
                            vec)

combined$cluster_labels <- relabelled_vector

Idents(combined) <- combined$cluster_labels

# Now plot with those labels
pdf("outputs/erythroid.pdf")
DimPlot(combined, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.5) +
  ggtitle("Defined Erythroid Subsets")
dev.off()

pdf("outputs/HS_enzyme_expr.pdf")

lapply(genes, function(gene) {

	p1 <- FeaturePlot(subset(combined, subset = group == "high"), features = gene, 
                  reduction = "umap") + ggtitle(paste(gene, "in HS3A8 High"))
	p2 <- FeaturePlot(subset(combined, subset = group == "low"), features = gene, 
	                  reduction = "umap") + ggtitle(paste(gene, "in HS3A8 Low"))

	p1 + p2 
})

dev.off()

mep_cluster <- subset(combined, subset = cluster_labels == "MEP")
mep_cluster %>% .[[4]] %>% table

mep_cluster <- NormalizeData(mep_cluster)
mep_cluster <- FindVariableFeatures(mep_cluster)
mep_cluster <- ScaleData(mep_cluster)
mep_cluster <- RunPCA(mep_cluster)
mep_cluster <- RunUMAP(mep_cluster, dims = 1:10)
mep_cluster <- FindNeighbors(mep_cluster, dims = 1:10)
mep_cluster <- FindClusters(mep_cluster, resolution = 0.5)

pdf("outputs/MEP.pdf")
DimPlot(mep_cluster, reduction = "umap", group.by = "group") + 
	ggtitle("MEP cluster")
dev.off()


pdf("outputs/MEP_enzymes.pdf")
lapply(genes, function(gene) {

	p1 <- FeaturePlot(mep_cluster, features = gene, 
                  reduction = "umap") + ggtitle(paste(gene, "in MEP cluster"))

	p1
})
dev.off()

clusters <- c("MEP", "EEP", "CEP")

dt <- combined@meta.data %>% data.table
dt[cluster_labels %in% clusters]$group %>% table
dt[cluster_labels %in% clusters]$cluster_labels %>% table

dt[cluster_labels %in% clusters][, .N, by = .(cluster_labels, group)]

plots <- lapply(clusters, function(cluster) {
	tmp_cluster <- subset(combined, subset = cluster_labels == cluster)
	tmp_cluster %>% .[[4]] %>% table

	tmp_cluster <- NormalizeData(tmp_cluster)
	tmp_cluster <- FindVariableFeatures(tmp_cluster)
	tmp_cluster <- ScaleData(tmp_cluster)
	tmp_cluster <- RunPCA(tmp_cluster)
	tmp_cluster <- RunUMAP(tmp_cluster, dims = 1:10)
	tmp_cluster <- FindNeighbors(tmp_cluster, dims = 1:10)
	tmp_cluster <- FindClusters(tmp_cluster, resolution = 0.5)

	
	# DimPlot(tmp_cluster, reduction = "umap", group.by = "group") + 
	# 	ggtitle(paste(cluster, "cluster"))

	for(gene in genes) {

	}

})

pdf("outputs/a.pdf")
print(plots)
dev.off()

####

pdf("outputs/test.pdf")
VlnPlot(subset(combined, subset = cluster_labels %in% clusters), features = "Hs6st1", pt.size = 0.1) + 
  ggtitle("Hs6st1 Expression Across Clusters")
dev.off()


pdf("outputs/test.pdf")
DotPlot(combined, features = genes, group.by = "cluster_labels") + RotatedAxis()
dev.off()












