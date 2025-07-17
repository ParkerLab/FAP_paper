library(sceasy)
library(reticulate)
setwd("~/links/christav/fap_village_multiome/results/RNA_clustering")
library(Seurat)
library(glue)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sctransform)
library(Matrix)
library(SingleCellExperiment)

rna <- readRDS("RNA_final.rds")

library(zellkonverter)
DefaultAssay(rna) <- "RNA"
rna <- JoinLayers(rna)
DefaultAssay(rna) <- "SCT"
rna_sce <- as.SingleCellExperiment(rna)
reducedDims(rna_sce)
reducedDimNames(rna_sce) <- gsub("UMAP", "X_umap", reducedDimNames(rna_sce))
writeH5AD(rna_sce, 'rna_final.h5ad')
