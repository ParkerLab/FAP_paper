#After making the Seurat object, adding metadata to the FUSION data.

library(Seurat)
library(miloR)
library(dplyr)
library(tidyr)

setwd("/scratch/scjp_root/scjp1/christav/milo/get_seurat_obj")

#reading in the Seurat object formed in 1_make_seurat.R
rna <- readRDS("FUSION_seurat.Rds")

indiv_metadata <- read.table("/scratch/scjp_root/scjp1/christav/milo/metadata/fusion_sample_info.tsv", header=T) #metadata on donors
nuc_metadata <- read.table("/scratch/scjp_root/scjp1/christav/milo/metadata/cluster_info_qc.tsv", header=T) # metadata on nuclei

metadata <- nuc_metadata[,c("index", "modality", "barcode", "SNG.1ST", "coarse_cluster_name")]
rna_metadata <- subset(metadata, modality == "rna")
rna_metadata2 <- subset(rna_metadata, SNG.1ST != "multiome")

rna_metadata3 <- rna_metadata2 %>% separate(index, c("rna", "batch", "NM", "batch2", "barcode2"), sep = "\\.")
rna_metadata3$seurat <- paste(rna_metadata3$barcode, "_", rna_metadata3$batch, "-NM-", rna_metadata3$batch2, ".decont", sep="")

rna_metadata4 <- merge(rna_metadata3, indiv_metadata, by="SNG.1ST")
small_metadata <- rna_metadata4[,c("SNG.1ST", "batch.x", "batch2", "barcode", "coarse_cluster_name", "seurat", "age", "sex", "bmi")]

metadata_sorted <- small_metadata[match(colnames(rna), small_metadata$seurat), ]
metadata_sorted <- metadata_sorted[, -which(names(metadata_sorted) == "seurat")]
rna <- AddMetaData(rna, metadata = metadata_sorted)
cells_to_keep <- !is.na(rna@meta.data$bmi)
rna_filtered <- rna[, cells_to_keep]

#adding ogtt, fasting glucose, and fasting insulin
big_meta <- read.csv("/scratch/scjp_root/scjp1/christav/milo/metadata/tissue.csv", header = T, sep = ",")
small_meta <- big_meta[,c("labelcode", "ogtt_status_biopsy", "glu_fast_biopsy")]
small_meta <- small_meta[complete.cases(small_meta),]
colnames(small_meta) <- c("SNG.1ST", "OGTT", "fasting_glucose")

fast_ins <- read.table("/scratch/scjp_root/scjp1/christav/milo/metadata/fasting_insulin_fusion_muscle_mRNA_301.txt", header = T, sep = "\t")
fast_ins <- fast_ins[,c("labelcode", "S_Insu")]
colnames(fast_ins) <- c("SNG.1ST", "fasting_insulin")

small_meta <- merge(small_meta, fast_ins)
small_meta <- small_meta %>%
  mutate(T2D = if_else(OGTT == "T2D", 1, 0))

seurat <- rna@meta.data$SNG.1ST
seurat <- as.data.frame(seurat)
colnames(seurat) <- "SNG.1ST"

meta_added <- merge(seurat, small_meta)
metadata_sorted <- meta_added[match(rna@meta.data$SNG.1ST, meta_added$SNG.1ST), ]
metadata_sorted <- metadata_sorted[, -which(names(metadata_sorted) == "SNG.1ST")]
rna <- AddMetaData(rna, metadata = metadata_sorted)
cells_to_keep <- !is.na(rna@meta.data$T2D)
rna_filtered <- rna[, cells_to_keep]

saveRDS(rna_filtered, "FUSION_seurat.Rds") #saving updated seurat object
