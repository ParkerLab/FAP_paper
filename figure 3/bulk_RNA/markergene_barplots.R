#use this to merge all of the bulk RNA together into one matrix, calculate TPM, and plot

setwd("/scratch/scjp_root/scjp1/christav/fap_bulk_rna/FAP_RNA/results/forRanalysis/")
library("DESeq2")
library("matrixStats")
library("dplyr")
library("biomaRt")
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)

cols <- c("Geneid", "id")
id_dds1 <- read.table("id_Rmatrix.txt", header = T, col.names = cols)

#merge all samples together
#need to make sure genes are in same order, but this should work
dds <- cbind(id_dds1, id_dds2)
rownames(dds) <- dds[,1]
dds[,1] <- NULL
dds[, "Geneid" == names(dds)] <- NULL

write.table(dds, "fusion_matrix_corrected.txt", row.names = TRUE, col.names = TRUE, quote = FALSE)

dds <- read.table("fusion_matrix_corrected.txt", header = TRUE)
dds <- dds[,c("id1", "id2")] 


#TPM normalization 
#need gene lengths
library(biomaRt)
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

genes <- rownames(dds)

gene_info <- getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id",
    "start_position",
    "end_position"
  ),
  filters = "hgnc_symbol",
  values = genes,
  mart = mart
)
gene_info$gene_length <- gene_info$end_position - gene_info$start_position
gene_length <- gene_info$gene_length
names(gene_length) <- gene_info$hgnc_symbol
gene_length <- gene_length[rownames(dds)]

#calculate tpm
rpk <- dds / (gene_length / 1000)
scaling_factor <- colSums(rpk, na.rm = TRUE)
tpm <- t( t(rpk) / scaling_factor ) * 1e6
log_tpm <- log2(tpm + 1)

#prep the matrix to plot
genes <- c("FBN1", "COL1A2", "LUM", "COL3A1", "POU5F1", "SOX2", "NANOG", "CD34")

genes_found <- genes[genes %in% rownames(tpm)]

df <- tpm[genes_found, ] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "donor", values_to = "expression")

df$gene <- factor(df$gene, levels = genes)

#plotting
ggplot(df, aes(x = gene, y = expression, fill = donor)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(
    y = "Log TPM",
    x = "Gene",
    title = "Gene Expression Across Donors"
  )
