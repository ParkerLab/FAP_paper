#This was used to generate the line plots of differential abundance by subcluster for different metadata traits in Supp Fig 1.

setwd("/scratch/scjp_root/scjp1/christav/milo/get_seurat_obj")
library(ggplot2)
library(dplyr)


#Step 1: Read in data for the subclusters
subclusters <- read.table("/nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/subclustering/subclustering_FAP_02_2025/results/factorize/fusion.bgstringent.Mesenchymal_Stem_Cell.numgenes_2000__k_5__lambda_5__epoch_5__iter_4__batchsize_5000.umap_cluster_res_0.025.tsv", header=T)
subclusters$barcode_match <- paste(subclusters$barcode, "_", subclusters$library, ".decont", sep="")
colnames(subclusters) <- c("UMAP_1", "UMAP_2", "modality", "library", "barcode_short", "batch", "extra", "lib1", "cluster", "barcode")


#Step 2: Read in barcodes and differential abundance values from milo.
barcodes <- read.table("barcodes_all_bmi.tsv", header=T)
barcodes_msc <- subset(barcodes, coarse_cluster_name == "Mesenchymal_Stem_Cell")


#Step 3: Combine them together, then separate out by subcluster.
heatmap <- merge(barcodes_msc, subclusters)
heatmap$cluster <- as.factor(heatmap$cluster)
range(heatmap$logFC)
pread_cluster <- subset(heatmap, cluster == "0")
fib_cluster <- subset(heatmap, cluster == "1")
prog_cluster <- subset(heatmap, cluster == "2")

pread_cluster <- pread_cluster[order(pread_cluster$logFC),]
pread_cluster$rank <- seq.int(nrow(pread_cluster))
pread_summary <- pread_cluster %>%
  mutate(percentile = ntile(rank, 10)) %>%
  group_by(percentile) %>%
  summarise(avg_logFC = mean(logFC, na.rm = TRUE))
pread_summary$cluster <- "pre-adipogenic"

fib_cluster <- fib_cluster[order(fib_cluster$logFC),]
fib_cluster$rank <- seq.int(nrow(fib_cluster))
fib_summary <- fib_cluster %>%
  mutate(percentile = ntile(rank, 10)) %>%
  group_by(percentile) %>%
  summarise(avg_logFC = mean(logFC, na.rm = TRUE))
fib_summary$cluster <- "pre-fibrogenic"

prog_cluster <- prog_cluster[order(prog_cluster$logFC),]
prog_cluster$rank <- seq.int(nrow(prog_cluster))
prog_summary <- prog_cluster %>%
  mutate(percentile = ntile(rank, 10)) %>%
  group_by(percentile) %>%
  summarise(avg_logFC = mean(logFC, na.rm = TRUE))
prog_summary$cluster <- "progenitors"


#Step 4: Recombine the subclusters and plot!
all_clusters <- rbind(pread_summary, fib_summary, prog_summary)
all_clusters$percentile.x <- all_clusters$percentile*10

ggplot(all_clusters, aes(x=percentile.x, y=avg_logFC, group=cluster, color=cluster)) + geom_hline(yintercept = 0) +
  geom_line() + geom_point() + ggtitle("Change in Differential Abundance by BMI for FAP Subclusters")
