#This script was for examining the percentage of FAPs in a muscle biopsy from the FUSION study, figure 1E and D.

library(dplyr)
library(tidyr)
library(ggplot2)
library("DESeq2")

setwd("/scratch/scjp_root/scjp1/christav/milo/percent_faps")


#Step 1: wrangling metadata
indiv_metadata <- read.table("/scratch/scjp_root/scjp1/christav/milo/metadata/fusion_sample_info.tsv", header=T) #metadata on each donor
nuc_metadata <- read.table("/scratch/scjp_root/scjp1/christav/milo/metadata/cluster_info_qc.tsv", header=T) #metadata on each nuclei 

metadata <- nuc_metadata[,c("index", "modality", "barcode", "SNG.1ST", "coarse_cluster_name")]
#rna_metadata <- subset(metadata, modality == "rna")
#rna_metadata2 <- subset(rna_metadata, SNG.1ST != "multiome")
#rna_metadata3 <- rna_metadata2 %>% separate(index, c("rna", "batch", "NM", "batch2", "barcode2"), sep = "\\.")

#if want just RNA, comment out the 2 lines below and run the 3 lines above
rna_metadata2 <- subset(metadata, SNG.1ST != "multiome")
rna_metadata3 <- rna_metadata2 %>% separate(index, c("rna", "batch", "NM", "batch2", "barcode2"), sep = "\\.")
rna_metadata3$seurat <- paste(rna_metadata3$barcode, "_", rna_metadata3$batch, "-NM-", rna_metadata3$batch2, ".decont", sep="")

rna_metadata4 <- merge(rna_metadata3, indiv_metadata, by="SNG.1ST")
small_metadata <- rna_metadata4[,c("SNG.1ST", "batch.x", "coarse_cluster_name", "age", "sex", "bmi")]

big_meta <- read.csv("/scratch/scjp_root/scjp1/christav/milo/metadata/tissue.csv", header = T, sep = ",")
small_meta <- big_meta[,c("labelcode", "ogtt_status_biopsy", "glu_fast_biopsy")]
small_meta <- small_meta[complete.cases(small_meta),]
colnames(small_meta) <- c("SNG.1ST", "OGTT", "fasting_glucose")
small_meta <- small_meta %>%
  mutate(T2D = if_else(OGTT == "T2D", 1, 0))

fast_ins <- read.table("/scratch/scjp_root/scjp1/christav/milo/metadata/fasting_insulin_fusion_muscle_mRNA_301.txt", header = T, sep = "\t") #fasting insulin metadata was separate
fast_ins <- fast_ins[,c("labelcode", "S_Insu")]
colnames(fast_ins) <- c("SNG.1ST", "fasting_insulin")

small_meta <- merge(small_meta, fast_ins)

indiv_metadata <- merge(indiv_metadata, small_meta)

cell_type_counts_individual <- small_metadata %>%
  group_by(SNG.1ST, coarse_cluster_name) %>%
  summarise(count = n(), .groups = 'drop')

wide_format <- cell_type_counts_individual %>%
  pivot_wider(names_from = SNG.1ST, values_from = count, values_fill = list(count = 0))
rows <- wide_format$coarse_cluster_name

rownames(indiv_metadata) <- indiv_metadata$SNG.1ST
indiv_metadata$SNG.1ST <- NULL

#checking format
all(rownames(indiv_metadata) %in% colnames(wide_format))
wide_format <- wide_format[, rownames(indiv_metadata)]
rownames(wide_format) <- rows
all(rownames(indiv_metadata) == colnames(wide_format))


#Step 2: using DESeq2 to test if there are differences in FAP composition
#centering and scaling all continuous variables, setting T2D as factor
indiv_metadata$bmi_scaled <- scale(indiv_metadata$bmi, center = TRUE)
indiv_metadata$age_scaled <- scale(indiv_metadata$age, center = TRUE)
indiv_metadata$fastglu_scaled <- scale(indiv_metadata$fasting_glucose, center = TRUE)
indiv_metadata$fastins_scaled <- scale(indiv_metadata$fasting_insulin, center = TRUE)
indiv_metadata$T2D <- as.factor(indiv_metadata$T2D)

#DESeq2 model
dds <- DESeqDataSetFromMatrix(countData=wide_format, 
                              colData=indiv_metadata, 
                              design=~batch + sex + bmi_scaled + age_scaled + T2D + fastglu_scaled + fastins_scaled)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
resOrdered
write.csv(as.data.frame(resOrdered), file="no_ogtt_all_sexM_bmi_age_t2d_fastglu_fastins.csv") #significant

#individual factors from model
resultsNames(dds)
res.sex <- results(dds, name = "sex_M_vs_F")
res.sex
write.csv(as.data.frame(res.sex), file="noogttallsexMbmiaget2dfastglufastins_sex.csv") #significant

res.bmi <- results(dds, name = "bmi_scaled")
res.bmi
write.csv(as.data.frame(res.bmi), file="noogttallsexMbmiaget2dfastglufastins_bmi.csv") #not significant

res.age <- results(dds, name = "age_scaled")
res.age
write.csv(as.data.frame(res.age), file="noogttallsexMbmiaget2dfastglufastins_age.csv") #not significant

res.t2d <- results(dds, name = "T2D_1_vs_0")
res.t2d
write.csv(as.data.frame(res.t2d), file="noogttallsexMbmiaget2dfastglufastins_t2d.csv") #significant

res.glu <- results(dds, name = "fastglu_scaled")
res.glu
write.csv(as.data.frame(res.glu), file="noogttallsexMbmiaget2dfastglufastins_fastglu.csv") #significant

res.ins <- results(dds, name = "fastins_scaled")
res.ins
write.csv(as.data.frame(res.ins), file="noogttallsexMbmiaget2dfastglufastins_fastins.csv") #significant


#violin plot of percent FAP by T2D status
cell_type_counts_individual_wide <- cell_type_counts_individual %>%
  pivot_wider(
    names_from = coarse_cluster_name,
    values_from = count,
    values_fill = 0  # fill missing combinations with 0
  )
indiv_metadata$SNG.1ST <- row.names(indiv_metadata)
cell_type_counts_individual_wide_meta <- merge(cell_type_counts_individual_wide, indiv_metadata, by = "SNG.1ST")
cell_type_counts_individual_wide_meta$T2D <- as.factor(cell_type_counts_individual_wide_meta$T2D)

ggplot(cell_type_counts_individual_wide_meta, aes(x=T2D, y=Mesenchymal_Stem_Cell, fill = T2D)) + geom_violin() +
  ggtitle("Number of FAPs in a Skeletal Muscle Biopsy by T2D Status")

cell_type_counts_individual_wide_meta$Total_Nuclei <- cell_type_counts_individual_wide_meta$Adipocyte+cell_type_counts_individual_wide_meta$Endothelial+cell_type_counts_individual_wide_meta$Macrophage+cell_type_counts_individual_wide_meta$Mesenchymal_Stem_Cell+cell_type_counts_individual_wide_meta$Muscle_Fiber_Mixed+cell_type_counts_individual_wide_meta$Neuromuscular_junction+cell_type_counts_individual_wide_meta$Neuronal+cell_type_counts_individual_wide_meta$Satellite_Cell+cell_type_counts_individual_wide_meta$Smooth_Muscle+cell_type_counts_individual_wide_meta$T_cell+cell_type_counts_individual_wide_meta$Type_1+cell_type_counts_individual_wide_meta$Type_2a+cell_type_counts_individual_wide_meta$Type_2x
cell_type_counts_individual_wide_meta$Percent_FAP <- (cell_type_counts_individual_wide_meta$Mesenchymal_Stem_Cell/cell_type_counts_individual_wide_meta$Total_Nuclei)*100

ggplot(cell_type_counts_individual_wide_meta, aes(x=T2D, y=Percent_FAP, fill = T2D)) + geom_violin() +
  ggtitle("Percent of FAPs in a Skeletal Muscle Biopsy by T2D Status")
