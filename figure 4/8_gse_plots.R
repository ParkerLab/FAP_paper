setwd("/gpfs/accounts/scjp_root/scjp1/christav/fap_village_multiome/results/mast/gse")

#go <- read.delim("/scratch/scjp_root/scjp99/vthihong/gsea_go_HV.txt", header = T, sep = "\t")
#kegg <- read.delim("/scratch/scjp_root/scjp99/vthihong/gsea_kegg_HV.txt", header = T, sep = "\t")

go <- read.table("gsea_go.final.txt", sep = "\t", header = T)
kegg <- read.table("gsea_kegg.final.txt", sep = "\t", header = T)

go <- go[order(go$p.adjust),]
kegg <- kegg[order(kegg$p.adjust),]

insulin <- subset(go, ID == "GO:0048009")
wnt <- subset(go, ID == "GO:0060070" | ID == "GO:0060828" | ID == "GO:0030177")
atp <- subset(go, ID == "GO:0042773" | ID == "GO:0042775" | ID == "GO:0006754")
morph <- subset(go, ID == "GO:0000902" | ID == "GO:0048729")
muscle <- subset(go, ID == "GO:0060537" | ID == "GO:0042692")
go_plot <- rbind(insulin, wnt, atp, morph, muscle)
go_plot <- go_plot[order(go_plot$p.adjust),]

#plotting
library(ggplot2)
ggplot(go_plot, aes(x=enrichmentScore, y=Description, color = p.adjust)) + geom_point(size=5)
