library(data.table)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)

# QC
counts <- readRDS("inputs/H1/counts_sparse.rds")

seu <- CreateSeuratObject(
  counts = counts,
  project = "scRNA",
  min.cells = 3,
  min.features = 200
)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-|^mt-")

pdf("outputs/01_H1-QC.pdf")

  VlnPlot(
    seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  )

  FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")

  seu_filt <- subset(
    seu,
    subset =
      nFeature_RNA > 500 &
      nFeature_RNA < 7000 &
      percent.mt < 20
  )

  VlnPlot(
    seu_filt,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  )

  FeatureScatter(seu_filt, feature1 = "nCount_RNA", feature2 = "percent.mt")

dev.off()

# saveRDS(seu, "intermediates/H1_seu_filtered.rds")

counts <- readRDS("inputs/H2/counts_sparse.rds")

seu <- CreateSeuratObject(
  counts = counts,
  project = "scRNA",
  min.cells = 3,
  min.features = 200
)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-|^mt-")

pdf("outputs/02_H2-QC.pdf")

  VlnPlot(
    seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  )

  FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")

  seu <- subset(
    seu,
    subset =
      nFeature_RNA > 500 &
      nFeature_RNA < 7000 &
      percent.mt < 20
  )

  VlnPlot(
    seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  )

  FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")

dev.off()

# saveRDS(seu, "intermediates/H2_seu_filtered.rds")

rm(list = c("counts", "seu", "seu_filt"))

#############
seu_H1 <- readRDS("intermediates/H1_seu_filtered.rds")
seu_H2 <- readRDS("intermediates/H2_seu_filtered.rds")

seu_H1$condition <- "preleukemic"
seu_H2$condition <- "leukemic"

seu_H1 <- NormalizeData(seu_H1)
seu_H2 <- NormalizeData(seu_H2)

objs <- list(seu_H1, seu_H2)

features <- SelectIntegrationFeatures(object.list = objs)

objs <- lapply(objs, function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

anchors <- FindIntegrationAnchors(object.list = objs, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)

DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:20)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

# saveRDS(combined, "intermediates/combined.rds")

combined <- readRDS("intermediates/combined.rds")

pdf("outputs/03_markers.pdf")

  DimPlot(combined)
  DimPlot(combined, group.by = "condition")
  FeaturePlot(combined, features = c("Kit","Ly6a","Procr","Hlf","Gata2"))
  VlnPlot(combined, features = c("Kit","Ly6a","Hlf"), group.by="seurat_clusters")

dev.off()

# Call clusters 
DefaultAssay(combined) <- "RNA"
combined <- JoinLayers(combined)
combined <- NormalizeData(combined)

cluster_markers <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

cluster_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>% print(n = 30)

lineage_markers <- list(
  HSC = c("Hlf","Procr","Kit","Ly6a","Gata2"),
  Myeloid = c("Mpo","Elane","S100a8","S100a9","Lyz2"),
  Macrophage = c("Csf1r","Trem2","Clec7a"),
  Erythroid = c("Gata1","Hbb-bt","Hba-a1"),
  Megakaryocyte = c("Pf4","Itga2b","Gp9"),
  Bcell = c("Cd79a","Ms4a1"),
  Tcell = c("Cd3d","Cd3e")
)

for (name in names(lineage_markers)) {
  combined <- AddModuleScore(
    combined,
    features = list(lineage_markers[[name]]),
    name = paste0(name,"_score")
  )
}

pdf("outputs/04_clusters.pdf")

  DimPlot(combined, label = TRUE)

  VlnPlot(
    combined,
    features = grep("_score1", colnames(combined@meta.data), value = TRUE),
    group.by = "seurat_clusters",
    pt.size = 0
  )

dev.off()

# Expression of HS Genes
HS_genes <- fread("inputs/HS_relevant_genes.csv") %>% .[[1]]

Idents(combined) <- "seurat_clusters"

genes_of_interest <- c("EXT1","Ext2","Ndst1","Hs3st1","Hlf") %>% stringr::str_to_title()
genes_of_interest <- HS_genes %>% stringr::str_to_title()

genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(combined)] %>% unique()

clusters <- levels(combined)

de_results <- list()

for (cl in clusters) {
  
  subset_obj <- subset(combined, idents = cl)
  Idents(subset_obj) <- "condition"
  
  de <- FindMarkers(
    subset_obj,
    ident.1 = "leukemic",
    ident.2 = "preleukemic",
    features = genes_of_interest,
    logfc.threshold = 0,
    min.pct = 0
  )
  
  de$cluster <- cl
  de$gene <- rownames(de)
  
  de_results[[cl]] <- de
}

de_table <- do.call(rbind, de_results)

de_table$cluster <- factor(
  de_table$cluster,
  levels = sort(unique(as.numeric(de_table$cluster)))
)

combined <- AddModuleScore(
  combined,
  features = list(genes_of_interest),
  name = "HS_score"
)

pdf("outputs/05_HS_genes.pdf")

  ggplot(de_table, aes(x = gene, y = cluster, fill = avg_log2FC)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0
    ) +
    labs(
      x = "Gene",
      y = "Cluster",
      fill = "log2FC\n(Leuk - Pre)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  DotPlot(
    combined,
    features = genes_of_interest,
    group.by = "seurat_clusters",
    split.by = "condition"
  ) + RotatedAxis()

  DotPlot(
    combined,
    features = genes_of_interest,
    group.by = "seurat_clusters",
    split.by = "condition"
  ) + RotatedAxis()

  DimPlot(combined, label = TRUE)
  FeaturePlot(combined, features = c("Hlf","Kit","Ly6a","Mpo","Csf1r"))

  VlnPlot(
    subset(combined, idents = 1),
    features = "HS_score1",
    group.by = "condition"
  )

  VlnPlot(combined,
        features = c("nFeature_RNA","nCount_RNA","percent.mt"),
        group.by = "seurat_clusters")

dev.off()

# Cluster 1 is a problem, get rid of it 

Idents(combined) <- "seurat_clusters"

combined_clean <- subset(combined, idents = setdiff(levels(combined), "1"))

table(combined_clean$seurat_clusters)

DefaultAssay(combined_clean) <- "RNA"
combined_clean <- JoinLayers(combined_clean)

clusters <- levels(combined_clean)
de_results <- list()

for (cl in clusters) {
  
  subset_obj <- subset(combined_clean, idents = cl)
  
  if (length(unique(subset_obj$condition)) < 2) next
  
  Idents(subset_obj) <- "condition"
  
  de <- FindMarkers(
    subset_obj,
    ident.1 = "leukemic",
    ident.2 = "preleukemic",
    features = genes_of_interest,
    logfc.threshold = 0,
    min.pct = 0
  )
  
  de$cluster <- cl
  de$gene <- rownames(de)
  
  de_results[[cl]] <- de
}

de_table_clean <- do.call(rbind, de_results)

de_table_clean$cluster <- factor(
  de_table_clean$cluster,
  levels = sort(unique(as.numeric(de_table_clean$cluster)))
)

combined_clean <- AddModuleScore(
  combined_clean,
  features = list(unique(genes_of_interest)),
  name = "HS_score"
)

pdf("outputs/06_HS_genes_clean.pdf")

  ggplot(de_table_clean, aes(x = gene, y = cluster, fill = avg_log2FC)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(
      x = "Gene",
      y = "Cluster",
      fill = "log2FC\n(Leuk - Pre)"
    )

  VlnPlot(
    combined_clean,
    features = "HS_score1",
    group.by = "condition"
  )

  VlnPlot(
    combined_clean,
    features = "HS_score1",
    group.by = "condition",
    split.by = "seurat_clusters"
  )

dev.off()

cluster9 <- subset(combined_clean, idents = 9)

pdf("outputs/07_clusters.pdf")

  VlnPlot(cluster9, features = "Extl1", group.by = "condition")

dev.off()

######
# percent expressing per gene
DefaultAssay(combined_clean) <- "RNA"

expr_matrix <- GetAssayData(combined_clean, layer = "data")

pct_expr <- Matrix::rowMeans(expr_matrix > 0)

df_expr <- data.frame(
  gene = names(pct_expr),
  pct_expressing = pct_expr
)

df_expr$category <- ifelse(df_expr$gene %in% genes_of_interest,
                           "HS_genes",
                           "Other_genes")

pdf("outputs/08_low_expression.pdf")

ggplot(df_expr, aes(x = pct_expressing, fill = category)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(
    x = "Fraction of cells expressing gene",
    y = "Density"
  )

expr_matrix <- GetAssayData(combined_clean, layer = "data")
mean_expr <- Matrix::rowMeans(expr_matrix)

df_expr$mean_expression <- mean_expr

ggplot(df_expr, aes(x = mean_expression, fill = category)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(
    x = "Mean normalized expression",
    y = "Density"
  )

dev.off()













##############
df_counts <- as.data.frame(
  table(combined_clean$seurat_clusters, combined_clean$condition)
)

colnames(df_counts) <- c("cluster", "condition", "n_cells")

pdf("outputs/04_clusters2.pdf", width = 12)

  ggplot(df_counts, aes(x = factor(cluster), y = n_cells, fill = condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(
      x = "Cluster",
      y = "Number of Cells",
      fill = "Condition"
    )

  VlnPlot(
    combined_clean,
    features = "HS_score1",
    group.by = "condition",
    split.by = "seurat_clusters"
  )

dev.off()










