#This code would plot marker genes identified though milo as being DE between FAP neighborhoods on the FAP subclusters, I tested a lot of genes so there's a long list commented out towards the end.

setwd("/scratch/scjp_root/scjp1/christav/milo/get_seurat_obj/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)
library(cowplot)
ggplot2::theme_set(theme_cowplot())
library(data.table)
library(viridis)
library(gridExtra)

#functions
savep = function(p, f, h=8, w=4){
  png(glue("{prefix}.{f}.png"), height=h, width=w, units="in", res = 150)
  print(p)
  dev.off()
}

expumap = function(d, g, title){
  p = ggplot(d, aes(x=UMAP_1, y=UMAP_2)) +
    geom_point(aes_string(color=g),shape=16, size=.4) +
    scale_color_viridis(direction=-1) +
    theme(#legend.position = "none", 
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()) +
    labs(title=title)
  return(p)
}    

makeplot = function(d, g){
  
  p = expumap(d, g, g)
  savep(p, glue("{g}.umap1"), 5,6)
  
  p = ggplot(d, aes_string(x="cluster", y=g)) +
    geom_violin(aes_string(color="cluster")) +
    geom_boxplot(aes_string(color="cluster")) +
    
    theme(legend.position = "none") +
    labs(title=g)
  savep(p, glue("{g}.violin1"), 3.5,7)
  
  
}

#loading FAP expression data
d = read.csv("/gpfs/accounts/scjp_root/scjp99/arushiv/muscle-sn/analyses_hg38/subclustering/subclustering_FAP_02_2025/results/normcounts/fusion.bgstringent.Mesenchymal_Stem_Cell.genes.tsv", sep='\t', header=T, as.is=T, row.names = 1)
d = log10(d +1)
d = as.data.frame(d)
d = t(d)
d = as.data.frame(d)
d$index = rownames(d)

#loading UMAP coordinates for the FAPs
c = read.csv("/gpfs/accounts/scjp_root/scjp99/arushiv/muscle-sn/analyses_hg38/subclustering/subclustering_FAP_02_2025/results/factorize/fusion.bgstringent.Mesenchymal_Stem_Cell.numgenes_2000__k_5__lambda_5__epoch_5__iter_4__batchsize_5000.umap_cluster_res_0.025.tsv", sep='\t', header=T, as.is=T, row.names = 1)
print(head(c))

#merging the expression data with the UMAP coordinates
c$index = rownames(c)
c$modality = gsub(".*_", "", c$modality)
c = c[c$modality == "rna",]
c = c[, c("UMAP_1","UMAP_2", "cluster", "index")]

d2 = merge(d, c, by = "index")
#print(head(d2))
genes = colnames(d2)
genes = genes[!genes %in% c("UMAP_1", "UMAP_2", "cluster", "filtering", "index")]

#lists of genes that I tested, and the metadata trait they were associated with (other were marker genes for other cell types as a sanity check)
#chart_genes = list("TNNT1" = "other", "MYH1" = "other", "ATP2A1" = "other", 
#                   "ATP2A2" = "other", "LGR5" = "other", "TNNT3" = "other", "TTN" = "other", "TPM3" = "other")

#chart_genes = list("CBLB" = "age", "CDK8" = "age", "SLC2A13" = "age", "PLCB4" = "age",
#                   "PDE5A" = "age", "ITGA8" = "age", "AFF3" = "age", "PLXDC1" = "age")

#chart_genes = list("ZBTB7C" = "BMI", "KCNQ5" = "BMI", "ITGA6" = "BMI", "TTN" = "BMI",
#                   "BMPR1B" = "BMI", "SLC22A3" = "BMI", "TRDN" = "BMI", "PDZRN4" = "BMI", "NEB" = "BMI")

#chart_genes = list("ARHGAP24" = "Fasting glucose", "NEB" = "Fasting glucose",
#                   "TRDN" = "Fasting glucose")

#chart_genes = list("ARHGEF28" = "Fasting insulin", "TENM2" = "Fasting insulin", 
#                   "LMO7" = "Fasting insulin",
#                   "PDZRN4" = "Fasting insulin", "TTN" = "Fasting insulin", 
#                   "BMPR1B" = "Fasting insulin", "NEB" = "Fasting insulin", "SLC22A3" = "Fasting insulin")

#chart_genes = list("ZNF385D" = "OGTT", "LAMA2" = "OGTT", "TP63" = "OGTT",
#                   "FAM129A" = "OGTT", "COL15A1" = "OGTT", "EGFR" = "OGTT", "PSD3" = "OGTT", "COL6A3" = "OGTT", "LRRTM4" = "OGTT", "RSPO3" = "OGTT")

#chart_genes = list("GREB1L" = "sex")

#chart_genes = list("HSPG2" = "T2D", "COL6A3" = "T2D", "RSPO3" = "T2D", 
#                   "LAMA2" = "T2D", "MME" = "T2D", "ARHGAP24" = "T2D", 
#                   "COL15A1" = "T2D", "SCN7A" = "T2D", "LRRTM4" = "T2D")

chart_genes = list("COL6A3" = "T2D", 
                   "LAMA2" = "T2D", "MME" = "T2D", 
                   "COL15A1" = "T2D", "SCN7A" = "T2D", "LRRTM4" = "T2D")


#CHECK GENES FIRST OR IT GETS STUCK
"SCN7A" %in% genes

centers <- d2 %>% group_by(.data[['cluster']]) %>% summarize(
  UMAP_1 = median(x = .data[['UMAP_1']]),
  UMAP_2 = median(x = .data[['UMAP_2']])
)

plist <- list()

## plotting on umap
plist[['clust']] = ggplot(d2, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(aes(color=cluster), shape=16, size=.4) +
  geom_text(data = centers, mapping = aes(label = cluster), colour = "black", size = 3) +
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

for (g in names(chart_genes)){
  plist[[g]] = expumap(d2, g, glue("{g}\n{chart_genes[[g]]}")) + theme(plot.title = element_text(size=7), axis.title = element_blank())
}

#saving
prefix = "t2d_fig"
png(glue("{prefix}.composite.png"), height=8, width=13, units="in", res=150)
do.call(grid.arrange, c(plist, nrow=2))
dev.off()
