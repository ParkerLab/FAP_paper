#Now clustering the Seurat object of FUSION data.

library(Seurat)
library(miloR)
library(dplyr)
library(tidyr)

setwd("/scratch/scjp_root/scjp1/christav/milo/get_seurat_obj")

rna_filtered <- readRDS("FUSION_seurat.Rds") #loading in the Seurat object with metadata added

rna_filtered <- NormalizeData(rna_filtered)
rna_filtered <- FindVariableFeatures(rna_filtered, selection.method = "vst", nfeatures = 2000)
rna_filtered <- ScaleData(rna_filtered)
rna_filtered <- RunPCA(rna_filtered, features = VariableFeatures(object = rna_filtered))
ElbowPlot(rna_filtered) #pick where levels off, 15
PCS <- 15
rna_filtered <- FindNeighbors(rna_filtered, dims = 1:PCS)
rna_filtered <- FindClusters(rna_filtered, resolution = 0.5)
rna_filtered <- RunUMAP(rna_filtered, dims = 1:10)
DimPlot(rna_filtered, reduction = "umap", group.by = "coarse_cluster_name", label = TRUE)
DimPlot(rna_filtered, reduction = "umap", group.by = "T2D", alpha = 0.5)

saveRDS(rna_filtered, "FUSION_seurat.Rds") #save updated Seurat object
