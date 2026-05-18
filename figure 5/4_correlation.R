setwd("~/links/christav/fap_village_multiome/results/RNA_clustering/comp_to_fusion")
library(Seurat)
library(glue)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sctransform)
library(Matrix.utils)
library(GenomicFeatures)

#getting gene expression matrix from differentiation data
rna <- readRDS("../RNA_final.rds")
#rna <- readRDS("../RNA_basal.rds")
#rna <- readRDS("../RNA_insulin.rds")
rna@meta.data$donor_clean <- gsub(".*\\(([^)]+)\\).*", "\\1", rna@meta.data$donor)
Idents(rna) <- "donor_clean"
raw_matrix <- GetAssayData(rna, slot = "counts")
pseudobulk_counts <- aggregate.Matrix(t(raw_matrix), groupings = Idents(rna), fun = "sum")
ipsc_fap <- t(pseudobulk_counts)
ipsc_fap <- as.data.frame(ipsc_fap)

#getting tpm values for ipsc-faps
## make TxDb from GTF file 
txdb <- makeTxDbFromGFF("/scratch/scjp_root/scjp0/shared_data/reference/human/hg38/topmed/gencode.v30.annotation.ERCC92.gtf")

## get gene information
all.genes <- genes(txdb)

#filtering for genes with at least 10 counts in 25% of donors
filter_by_expression <- function(df, min_count = 10, min_fraction = 0.25) {
  keep <- rowSums(df >= min_count) >= (ncol(df) * min_fraction)
  df[keep, ]
}
ipsc_fap_filtered <- filter_by_expression(ipsc_fap)

#now getting gene lengths for tpm calculation
ipsc_fap_filtered$genes <- rownames(ipsc_fap_filtered)
ipsc_fap_filtered <- ipsc_fap_filtered %>% separate(genes, c("ENS_ID_ext","Name"), sep = " ")
my.genes <- ipsc_fap_filtered$ENS_ID_ext
ipsc_fap_filtered$Name <- gsub(".*\\(([^)]+)\\).*", "\\1", ipsc_fap_filtered$Name)

## get the length of each of those genes
my.genes.lengths <- width(all.genes[my.genes])
## put the names back on the lengths
names(my.genes.lengths) <- my.genes
ipsc_fap_filtered$gene_length <- my.genes.lengths

ipsc_fap_filtered$gene_lengths_kb <- ipsc_fap_filtered$gene_length / 1000

gene_length_kb <- ipsc_fap_filtered$gene_lengths_kb
ipsc_fap_filtered$ENS_ID_ext <- NULL
ipsc_fap_filtered$Name <- NULL
ipsc_fap_filtered$gene_length <- NULL
ipsc_fap_filtered$gene_lengths_kb <- NULL

rpk <- sweep(ipsc_fap_filtered, 1, gene_length_kb, FUN = "/")
scaling_factors <- colSums(rpk)
tpm <- sweep(rpk, 2, scaling_factors, FUN = "/") * 1e6
tpm <- as.data.frame(tpm)

#average tpm
tpm$mean <- rowSums(tpm)/29
#add gene names back
tpm$genes <- rownames(tpm)
tpm <- tpm %>% separate(genes, c("ENS_ID_ext","Name"), sep = " ")
tpm$Name <- gsub(".*\\(([^)]+)\\).*", "\\1", tpm$Name)


#getting gene expression data from FUSION
fusion <- read.table("/nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/eqtl/eqtl_11/results/input-features/fusion.Adipocyte.tsv", header = F, sep = "\t")
#filtering for genes with at least 10 counts in 25% of donors
rownames(fusion) <- fusion$V4
fusion$gene_lengths_kb <- fusion$V7 / 1000
fusion_genes_lengths <- cbind(fusion$V4, fusion$gene_lengths_kb)
fusion$V1 <- NULL
fusion$V2 <- NULL
fusion$V3 <- NULL
fusion$V4 <- NULL
fusion$V5 <- NULL
fusion$V6 <- NULL
fusion$V7 <- NULL
fusion$gene_lengths_kb <- NULL
fusion_filtered <- filter_by_expression(fusion)
fusion_filtered$V1 <- rownames(fusion_filtered)
fusion_filtered <- merge(fusion_filtered, fusion_genes_lengths)
fusion_gene_length_kb <- fusion_filtered$V2
fusion_gene_length_kb <- as.numeric(fusion_gene_length_kb)
rownames(fusion_filtered) <- fusion_filtered$V1
fusion_filtered$V1 <- NULL
fusion_filtered$V2 <- NULL
#tpm per donor for fusion
fusion_rpk <- sweep(fusion_filtered, 1, fusion_gene_length_kb, FUN = "/")
fusion_scaling_factors <- colSums(fusion_rpk)
fusion_tpm <- sweep(fusion_rpk, 2, fusion_scaling_factors, FUN = "/") * 1e6
fusion_tpm <- as.data.frame(fusion_tpm)

fusion_tpm$mean <- rowSums(fusion_tpm)/ncol(fusion_tpm)
fusion_tpm$genes <- rownames(fusion_tpm)

#subset to same genes
common_genes <- intersect(tpm$Name, fusion_tpm$genes)
ipsc_fap_small <- tpm[c("Name", "mean")]
colnames(ipsc_fap_small) <- c("Gene", "ipsc_fap_tpm")
fusion_small <- fusion_tpm[c("genes", "mean")]
colnames(fusion_small) <- c("Gene", "fusion_tpm")

ipsc_fap_common <- ipsc_fap_small[ipsc_fap_small$Gene %in% common_genes, ]
fusion_fap_common <- fusion_small[fusion_small$Gene %in% common_genes, ]
ipsc_fap_common <- ipsc_fap_common[order(ipsc_fap_common$Gene), ]
fusion_fap_common <- fusion_fap_common[order(fusion_fap_common$Gene), ]
ipsc_fap_common <- ipsc_fap_common[!duplicated(ipsc_fap_common$Gene), ]
fusion_fap_common <- fusion_fap_common[!duplicated(fusion_fap_common$Gene), ]

ipsc_fap_common$fusion_tpm <- fusion_fap_common$fusion_tpm

#calculate correlation
correlation_result <- cor.test(ipsc_fap_common$ipsc_fap_tpm, ipsc_fap_common$fusion_tpm, method = "spearman")
print(correlation_result)
print(length(common_genes))

#insulin
#FAP: 0.59391 p-value < 2.2e-16, number of common genes = 3005
#Smooth muscle: 0.5888085 p-value < 2.2e-16, number of common genes = 2124
#Endothelial: 0.5445251 p-value < 2.2e-16, number of common genes = 2831
#Satellite cells: 0.5371937 p-value < 2.2e-16, number of common genes = 645
#NMJ: 0.5291458 p-value < 2.2e-16, number of common genes = 1263
#MFM: 0.4869289 p-value < 2.2e-16, number of common genes = 218
#Type 2x: 0.4663332 p-value < 2.2e-16, number of common genes = 5339
#Type 2a: 0.4624044 p-value < 2.2e-16, number of common genes = 6190
#Type 1: 0.4275677 p-value < 2.2e-16, number of common genes = 7454
#Neuronal: 0.407066 p-value = 7.624e-08, number of common genes = 165
#Macrophage: 0.384056 p-value = 2.572e-10, number of common genes = 257
#Adipocyte: 0.3300548 p-value = 2.085e-06, number of common genes = 200

#Tcells: 0.8928571 p-value = 0.0123, number of common genes = 7


#basal
#FAP: 0.6443935 p-value < 2.2e-16, number of common genes = 2931
#Smooth muscle: 0.6236363 p-value < 2.2e-16, number of common genes = 2064
#Endothelial: 0.5786001 p-value < 2.2e-16, number of common genes = 2743
#Satellite_cells: 0.5709351 p-value < 2.2e-16, number of common genes = 628
#NMJ: 0.5420168 p-value < 2.2e-16, number of common genes = 1221
#Type 2x: 0.5035689 p-value < 2.2e-16, number of common genes = 5023
#MFM: 0.4886154 p-value < 2.2e-16, number of common genes = 211
#Type 2a: 0.483594 p-value < 2.2e-16, number of common genes = 5729
#Neuronal: 0.4543908 p-value = 4.282e-09, number of common genes = 155
#Type 1: 0.4263726 p-value < 2.2e-16, number of common genes = 6709
#Adipocyte: 0.416834 p-value = 2.566e-09, number of common genes = 192
#Macrophage: 0.3932289 p-value = 1.739e-10, number of common genes = 249

#Tcells: 0.8928571 p-value = 0.0123, number of common genes = 7


#basal and insulin
#FAP: 0.5911684 p-value < 2.2e-16, number of common genes = 3084
#Smooth_Muscle: 0.5782736 p-value < 2.2e-16, number of common genes = 2170
#Endothelial: 0.5268661 p-value < 2.2e-16, number of common genes = 2908
#Satellite cell: 0.5080873 p-value < 2.2e-16, number of common genes = 663
#NMJ: 0.5060416 p-value < 2.2e-16, number of common genes = 1298
#Type 2a: 0.466226 p-value < 2.2e-16, number of common genes = 6550
#Type 2x: 0.462895 p-value < 2.2e-16, number of common genes = 5602
#MFM: 0.4502332 p-value = 6.451e-13, number of common genes = 228
#Type 1: 0.4466631 p-value < 2.2e-16, number of common genes = 8053
#Neuronal: 0.3744287 p-value = 3.496e-07, number of common genes = 177
#Macrophage: 0.3423811 p-value = 1.593e-08, number of common genes = 262
#Adipocyte: 0.2461619 p-value = 0.0003184, number of common genes = 211

#T cell: 0.8928571 p-value = 0.0123, number of common genes = 7 (dropped due to common genes being less than 10)



cell_types <- c("Adipocyte", "Macrophage", "Neuronal", "Type_1", "MFM", "Type_2x", "Type_2a", "NMJ", "Satellite_cell", "Endothelial", "Smooth_muscle", "FAP")
correlation <- c(0.2461619, 0.3423811, 0.3744287, 0.4466631, 0.4502332, 0.462895, 0.466226, 0.5060416, 0.5080873, 0.5268661, 0.5782736, 0.5911684)

df <- rbind(cell_types, correlation)
df <- as.data.frame(df)
df <- t(df)
df <- as.data.frame(df)
df$correlation <- as.numeric(df$correlation)
df$norm_coeff <- df$correlation/max(df$correlation)
df$cell_types <- factor(df$cell_types, levels=c("Adipocyte", "Macrophage", "Neuronal", "Type_1", "MFM", "Type_2x", "Type_2a", "NMJ", "Satellite_cell", "Endothelial", "Smooth_muscle", "FAP"))

library(ggplot2)
ggplot(df, aes(x=cell_types, y=1, fill=norm_coeff)) + geom_tile(color = NA) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Reference Cell Type") + ylab("iPSC-FAPs") + labs(fill = "Normalized Spearman Correlation Coefficient") +
  scale_fill_gradientn(values = scales::rescale(x=c(min(df$norm_coeff),0.95, max(df$norm_coeff)), to = c(0,1), from = c(min(df$norm_coeff), max(df$norm_coeff))),  colors = c('white', "white", "red"))

