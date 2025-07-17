setwd("/gpfs/accounts/scjp_root/scjp1/christav/fap_village_multiome/results/scDRS/enrichment_per_cluster")

rna <- readRDS("../../RNA_clustering/RNA_final.rds")
DimPlot(rna, reduction = "umap", group.by = "subclusters", label = TRUE)

scores <- read.table("../compute_score_out/UKB_460K.biochemistry_HbA1c.score.gz", header = T, sep = "\t")

scores_small <- scores[,c("X", "norm_score")]
colnames(scores_small) <- c("barcode", "score") 
#range(scores_small$score)
scores_small$score <- scores_small$score + 5

rna$subcluster <- as.character(rna$subcluster)
clusters <- cbind(colnames(rna), rna$subcluster)
clusters <- as.data.frame(clusters)
colnames(clusters) <- c("barcode", "cluster")

clusters_scores <- merge(clusters, scores_small, by = "barcode")
clusters_prog <- subset(clusters_scores, cluster == "Progenitors")
clusters_adipo <- subset(clusters_scores, cluster == "Adipogenic")
clusters_fibro <- subset(clusters_scores, cluster == "Fibrogenic")

library(dplyr)

avg_scores_subcluster <- clusters_scores %>%
  group_by(cluster) %>%
  summarise(mean_score = mean(score, na.rm = TRUE))

avg_scores_subcluster$mean_score <- avg_scores_subcluster$mean_score - 5

library(ggplot2)
ggplot(avg_scores_subcluster, aes(x = reorder(cluster, -mean_score), y = mean_score, fill = cluster)) + geom_bar(stat = "identity") +
  ggtitle("Basal Metabolic Rate Enrichment by Subcluster") + xlab("FAP Subcluster") + ylab("Mean Enrichment Score per Nuclei")

t.test(clusters_fibro$score, clusters_adipo$score)
t.test(clusters_prog$score, clusters_adipo$score)
t.test(clusters_fibro$score, clusters_prog$score)

#T2D: fibro/adipo pval = 0.09218; adipo/prog pval = 0.0001109
#HbA1c: fibro/prog pval = 5.173e-06; prog/adipo pval = 0.008152
#triglycerides: fibro-adipo pval = < 2.2e-16; adipo/prog pval = 7.627e-07
#cardio: prog/adipo pval = 1.362e-09; adipo/fibro pval = < 2.2e-16
#basal metabolic: prog/adipo pval = < 2.2e-16; adipo/fibro pval = < 2.2e-16



#code below generated a big table with average enrichment per subcluster for all of the traits that I had, took some time to run but was handy to look at many traits at once

rna$subcluster <- as.character(rna$subcluster)
clusters <- cbind(colnames(rna), rna$subcluster)
clusters <- as.data.frame(clusters)
colnames(clusters) <- c("barcode", "cluster")

# Read the list of score files
score_files <- read.table("/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/scDRS/compute_score_out/files.txt", header = TRUE)
score_files <- as.character(score_files$files.txt)

# Initialize an empty list to hold average score data frames
scores_list <- list()

for (filename in score_files) {
  # Read the score file
  full_path <- file.path("..", "compute_score_out", filename)
  score_df <- read.table(full_path, header = TRUE, sep = "\t")
  
  # Prepare the score data
  scores_small <- score_df[, c("X", "norm_score")]
  colnames(scores_small) <- c("barcode", "score") 
  scores_small$score <- scores_small$score + 5
  
  # Merge with clusters data
  clusters_scores <- merge(clusters, scores_small, by = "barcode")
  
  # Compute average score per subcluster
  avg_scores_subcluster <- clusters_scores %>%
    group_by(cluster) %>%
    summarise(!!filename := mean(score, na.rm = TRUE) - 5, .groups = 'drop')
  
  # Add to list
  scores_list[[filename]] <- avg_scores_subcluster
}

# Merge all trait-wise average scores by cluster
all_scores <- Reduce(function(x, y) merge(x, y, by = "cluster", all = TRUE), scores_list)

write.table(all_scores, "all_scores_final.tsv", row.names = F, col.names = T, quote = F, sep = "\t")
