nuc_metadata <- read.table("/scratch/scjp_root/scjp1/christav/milo/metadata/cluster_info_qc.tsv", header=T)
barcodes <- nuc_metadata[,c("barcode", "SNG.1ST")]
write.table(barcodes, "barcode_assignments.tsv", col.names = F, row.names = F, quote = F, sep = "\t")

setwd("/gpfs/accounts/scjp_root/scjp1/christav/sig_peak_test")
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load counts
counts <- read.delim("peak_counts_indiv.txt", comment.char="#")
rownames(counts) <- counts$Geneid
count_data <- counts[, -(1:6)]  # remove annotation columns

# Sample metadata
metadata <- colnames(count_data)
metadata <- data.frame(metadata)
metadata <- metadata %>% separate(metadata, c("folder", "meta", "bam"), sep = "\\.")
metadata <- metadata %>% separate(meta, c("subcluster", "donor"), sep = "_")
metadata <- metadata[,c("donor", "subcluster")]

# filter genes: keep if they have >=20 counts in >=60 (little over 20%) samples
keep <- rowSums(count_data >= 20) >= 60
count_data_filtered <- count_data[keep, ]

peaks_per_sample <- colSums(count_data_filtered > 0)
peaks_per_sample <- data.frame(peaks_per_sample)
ggplot(peaks_per_sample, aes(x=peaks_per_sample)) + geom_histogram()

#filter samples, only keep individuals with more than 20,000 peaks
peaks_per_sample <- colSums(count_data_filtered > 0)
threshold <- 20000
good_samples <- peaks_per_sample >= threshold
count_data_filtered <- count_data_filtered[, good_samples]

#add filter for individuals with all 3
mtadata2 <- metadata[good_samples, ]
all_3_samples <- table(mtadata2$donor)
all_3_samples <- names(all_3_samples[all_3_samples == 3])

keep <- mtadata2$donor %in% all_3_samples 
count_data_filtered <- count_data_filtered[,keep]
metadata_filtered <- mtadata2[keep,]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_data_filtered,
  colData = metadata_filtered,
  design = ~ donor + subcluster
)

# Run DESeq2
dds <- DESeq(dds)
#peak_154762     chr3    155078965       155080492
resultsNames(dds)
res.adipose.fibrous <- results(dds, name = "subcluster_preadipogenic_vs_fibrous")
res.prog.fibrous <- results(dds, name = "subcluster_progenitors_vs_fibrous")
res.adipose.prog <- results(dds, contrast=c("subcluster","preadipogenic","progenitors"))

res.adipose.fibrous["peak_154762", "pvalue"] #2.711978e-19
res.prog.fibrous["peak_154762", "pvalue"] #0.2085152
res.adipose.prog["peak_154762", "pvalue"] #2.278215e-31

res.adipose.fibrous["peak_154762", "padj"] #2.950771e-25
res.prog.fibrous["peak_154762", "padj"] #0.3166537
res.adipose.prog["peak_154762", "padj"] #2.178749e-29
