setwd("/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/RNA_clustering/")

library(Seurat)
library(miloR)
library(dplyr)
library(tidyr)
library(SingleCellExperiment)
library(scater)
library(patchwork)
library(Matrix)
library(BiocParallel)
library(scran)
set.seed(42)

rna <- readRDS("RNA_final.rds")

#modifying donor names to separate out insulin and basal
rna$donor <- gsub(".*\\(([^)]+)\\).*", "\\1", rna$donor)
rna$donor_enviro <- paste(rna$donor, rna$enviro, sep="_")

DefaultAssay(rna) <- "RNA"
rna <- JoinLayers(rna)
DefaultAssay(rna) <- "SCT"
rna_sce <- as.SingleCellExperiment(rna)
rna_milo <- Milo(rna_sce)

rna_milo <- buildGraph(rna_milo, k = 40, d = 25) #got values from Arushi's attempt at this
rna_milo <- makeNhoods(rna_milo, prop = 0.05, k = 40, d=25, refined = TRUE)
plotNhoodSizeHist(rna_milo) #want between 50-100
rna_milo <- countCells(rna_milo, meta.data = data.frame(colData(rna_milo)), sample="donor_enviro")

saveRDS(rna_milo, "ipscfap_milo.Rds")
rna_milo <- readRDS("ipscfap_milo.Rds")

traj_design <- data.frame(colData(rna_milo))[,c("donor_enviro", "age", "bmi", "ogtt", "sex", "enviro")] #change based on what testing for
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$donor_enviro

rna_milo <- calcNhoodDistance(rna_milo, d=25)
da_results <- testNhoods(rna_milo, design = ~ age + bmi + ogtt + sex + enviro, design.df = traj_design)
rna_milo <- buildNhoodGraph(rna_milo)

plotUMAP(rna_milo) + plotNhoodGraphDA(rna_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")

da_results <- annotateNhoods(rna_milo, da_results, coldata_col = "subcluster")
plotDAbeeswarm(da_results, group.by = "subcluster") 

plotDAbeeswarm(da_results, group.by = "subcluster") + #ylim(-2, 2) +
  scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, na.value = "lightgrey")
