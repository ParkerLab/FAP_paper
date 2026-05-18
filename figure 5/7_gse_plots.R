setwd("/gpfs/accounts/scjp_root/scjp1/christav/fap_village_multiome/results/mast/gse")

#go <- read.delim("/scratch/scjp_root/scjp99/vthihong/gsea_go_HV.txt", header = T, sep = "\t")
#kegg <- read.delim("/scratch/scjp_root/scjp99/vthihong/gsea_kegg_HV.txt", header = T, sep = "\t")

go <- read.table("gsea_go.0.35.txt", sep = "\t", header = T)
kegg <- read.table("gsea_kegg.0.35.txt", sep = "\t", header = T)

go <- go[order(go$p.adjust),]
kegg <- kegg[order(kegg$p.adjust),]

#insulin <- subset(go, ID == "GO:0048009")
#glucose
#lipids
wnt <- subset(go, ID == "GO:0016055" | ID == "GO:0060070" | ID == "GO:0030111")
atp <- subset(go, ID == "GO:0006754" | ID == "GO:0046034" | ID == "GO:0042773")
morph <- subset(go, ID == "GO:0000902" | ID == "GO:0048729")
muscle <- subset(go, ID == "GO:0061061" | ID == "GO:0042692")
go_plot <- rbind(wnt, atp, morph, muscle)
#go_plot <- go_plot[order(go_plot$p.adjust),]

#collagen <- subset(go, ID == "GO:0030199" | ID == "GO:0032963" | ID == "GO:0038065")
#ecm <- subset(go, ID == "GO:0030198" | ID == "GO:0085029")
#go_plot <- rbind(collagen, ecm)
#go_plot <- go_plot[order(go_plot$p.adjust),]

#cell_diff <- subset(go, ID == "GO:0045597")
#ecm <- subset(go, ID == "GO:0030198")
#mes_cell_diff <- subset(go, ID == "GO:0048762")
#go_plot <- rbind(cell_diff, ecm, mes_cell_diff)
#go_plot <- go_plot[order(go_plot$p.adjust),]

kegg_plot <- subset(kegg, ID == "hsa04931" | ID == "hsa00190" | ID == "hsa03040" | ID == "hsa04820" | ID == "hsa04512" | ID == "hsa04151")

go_plot$Description <- factor(go_plot$Description, levels=c("muscle cell differentiation", "muscle structure development", "cell morphogenesis", "tissue morphogenesis", "ATP metabolic process", "ATP biosynthetic process", "ATP synthesis coupled electron transport", "regulation of Wnt signaling pathway", "canonical Wnt signaling pathway", "Wnt signaling pathway"))
kegg_plot$Description <- factor(kegg_plot$Description, levels=c("Cytoskeleton in muscle cells", "ECM-receptor interaction", "Oxidative phosphorylation", "PI3K-Akt signaling pathway", "Spliceosome", "Insulin resistance"))

#plotting
library(ggplot2)
ggplot(go_plot, aes(x=enrichmentScore, y=Description, color = pvalue)) + geom_point(size=5) 
ggplot(kegg_plot, aes(x=enrichmentScore, y=Description, color = pvalue)) + geom_point(size = 5)


#KEGG genes
library(KEGGREST)
library(biomaRt)
library(dplyr)

# KEGG numeric gene IDs from core_enrichment
nums <- c("2932","4790","5728","9882","5599","4792","10000","6197","2673","10724")

# Add prefix
kegg_ids <- paste0("hsa:", nums)

# Step 1: Extract Entrez IDs from KEGG
entrez_ids <- sapply(kegg_ids, function(x) {
  entry <- tryCatch(keggGet(x)[[1]], error = function(e) NULL)
  if (is.null(entry)) return(NA)
  sub("hsa:", "", entry$ENTRY) # returns Entrez ID
})

# Build small data.frame
df <- data.frame(
  KEGG_ID = kegg_ids,
  Entrez = as.integer(entrez_ids),  
  stringsAsFactors = FALSE
)

# Step 2: Map Entrez → Ensembl using biomaRt
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(
  attributes = c("entrezgene_id", "ensembl_gene_id", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = df$Entrez,
  mart = mart
)

# Step 3: Join (now works fine)
final <- left_join(df, mapping, by = c("Entrez" = "entrezgene_id"))

final

ens <- final$ensembl_gene_id


