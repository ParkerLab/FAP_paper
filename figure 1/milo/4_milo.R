#Generalizable code, I tweaked this depending on what metadata variables I was interested in testing for.
#This gives how to get Figure 1E-G and Supp Fig 1A-J

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

setwd("/scratch/scjp_root/scjp1/christav/milo/get_seurat_obj")


#Step 1: convert the Seurat object to a milo object. Can save this afterwards and then comment this part out.
rna_filtered <- readRDS("FUSION_seurat.rds") #load FUSION Seurat object that has been clustered.
rna_filtered_sce <- as.SingleCellExperiment(rna_filtered) #convert to SCE
rna_milo <- Milo(rna_filtered_sce) #convert to milo
rna_milo <- buildGraph(rna_milo, k = 40, d = 25) 
rna_milo <- makeNhoods(rna_milo, prop = 0.05, k = 40, d=25, refined = TRUE)
plotNhoodSizeHist(rna_milo) 
rna_milo <- countCells(rna_milo, meta.data = data.frame(colData(rna_milo)), sample="SNG.1ST")
saveRDS(rna_milo, "FUSION_milo.Rds")


#Step 2: running milo!
rna_milo <- readRDS("FUSION_milo.Rds")

traj_design <- data.frame(colData(rna_milo))[,c("SNG.1ST", "age", "batch.x", "bmi", "fasting_glucose", "fasting_insulin", "T2D", "sex")] #change based on what testing for, last item is the one being tested (in this case, sex)
traj_design$T2D <- as.factor(traj_design$T2D)
colData(rna_milo)$T2D <- as.factor(colData(rna_milo)$T2D)
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$SNG.1ST

rna_milo <- calcNhoodDistance(rna_milo, d=25)
da_results <- testNhoods(rna_milo, design = ~ age + batch.x + bmi + fasting_glucose + fasting_insulin + T2D + sex, design.df = traj_design) #change based on what testing for, last item is the one being tested (in this case, sex)
rna_milo <- buildNhoodGraph(rna_milo)


#Step 3: saving plots
jpeg(file="sex_all_milo.jpeg", res=600, width=9600, height=4800)
plotUMAP(rna_milo) + plotNhoodGraphDA(rna_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")
dev.off()
da_results <- annotateNhoods(rna_milo, da_results, coldata_col = "coarse_cluster_name")
jpeg(file="sex_all_beeswarm.jpeg", res=600, width=4800, height=9600)
plotDAbeeswarm(da_results, group.by = "coarse_cluster_name")
dev.off()
jpeg(file="sex_all_beeswarm_fig_color.jpeg", res=600, width=4800, height=9600)
plotDAbeeswarm(da_results, group.by = "coarse_cluster_name") + ylim(-2, 2) +
  scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, na.value = "lightgrey") #adjusting axis and colors
dev.off()
da_results <- groupNhoods(rna_milo, da_results, max.lfc.delta = 10, overlap=1) 


#Step 4: getting barcodes of nuclei in different neighborhoods
nhoods_to_nhood_index_df <- data.frame(
  Nhood = seq_along(as.character(nhoodIndex(rna_milo))),  
  NhoodIndex = as.character(nhoodIndex(rna_milo))
)

# Let's convert the barcode x cell index nHood sparse matrix to a regular dataframe that maps
# barcode to cell index, which we can then map to nHood name, which we can then match to
# the da_results info. 
long_df <- as(nhoods(rna_milo), "TsparseMatrix") %>%
  {data.frame(
    barcode = rownames(.)[.@i + 1],  
    NhoodIndex = colnames(.)[.@j + 1],
    value = .@x
  )} %>%
  filter(value > 0) %>% # Only keep non-zero entries if needed
    select(-value)

barcode_index_mapping <- left_join(long_df, nhoods_to_nhood_index_df, by = join_by("NhoodIndex" == "NhoodIndex"))
final_df <- left_join(barcode_index_mapping, da_results, by = join_by("Nhood" == "Nhood"))
write.table(final_df, "barcodes_all_sex.tsv", quote = FALSE, row.names = FALSE, sep = "\t")


#Step 5: getting marker genes for neighborhoods
#remove 0 count genes
logcounts(rna_milo) <- log1p(counts(rna_milo))
keep.rows <- rowSums(logcounts(rna_milo)) != 0
rna_milo <- rna_milo[keep.rows, ]
## Find highly variable genes
set.seed(42)
dec <- modelGeneVar(rna_milo)
hvgs <- getTopHVGs(dec, n=2000)

dge_smp <- findNhoodGroupMarkers(rna_milo, da_results, subset.row = hvgs, subset.groups = c("10"))
write.table(dge_smp, "de_all_sex.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

#testing within FAPs, change this as well for whatever variable interested in
dge_12 <- testDiffExp(rna_milo, da_results, design = ~ sex, meta.data = data.frame(colData(rna_milo)),
                     subset.row = hvgs, subset.nhoods=da_results$NhoodGroup=="10")
dge_12 <- as.data.frame(dge_12)
dge_12$GeneID <- row.names(dge_12)
write.table(dge_12, "de_mscs_sex.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

#plot heatmap of DE genes
dge_smp <- dge_smp[order(dge_smp$adj.P.Val_10),]

markers <- dge_smp$GeneID[1:10]

genes <- intersect(rownames(rna_milo), markers)
rownames(rna_milo) <- gsub(" ", "_", rownames(rna_milo))
rownames(rna_milo) <- gsub("\\(", "", rownames(rna_milo))
rownames(rna_milo) <- gsub("\\)", "", rownames(rna_milo))

genes <- gsub(" ", "_", genes)
genes <- gsub("\\(", "", genes)
genes <- gsub("\\)", "", genes)

jpeg(file="sex_marker_heatmap.jpeg", res=600, width=4800, height=4800)
plotNhoodExpressionGroups(rna_milo, da_results, genes,
                          subset.nhoods = da_results$NhoodGroup %in% c('10'), 
                          scale=TRUE,
                          grid.space = "fixed")
dev.off()

#testing within FAPs
dge_12 <- dge_12[order(dge_12$X10.adj.P.Val),]

markers <- dge_12$GeneID[1:10]
markers <- gsub(" ", "_", markers)
markers <- gsub("\\(", "", markers)
markers <- gsub("\\)", "", markers)

jpeg(file="sex_marker_mscs_heatmap.jpeg", res=600, width=4800, height=4800)
plotNhoodExpressionGroups(rna_milo, da_results, markers,
                          subset.nhoods = da_results$coarse_cluster_name %in% c("Mesenchymal_Stem_Cell"), 
                          scale=TRUE,
                          grid.space = "fixed")
dev.off()
