setwd("/gpfs/accounts/scjp_root/scjp1/christav/fall22_multiome/demuxlet")

demuxlet <- read.table("/scratch/scjp_root/scjp1/christav/fall22_multiome/demuxlet/E-hg38-ATAC.best.txt", header = T, sep = "\t")

demuxlet.sng <- subset(demuxlet, DROPLET.TYPE == "SNG")
barcode_indiv <- demuxlet.sng[c("BARCODE", "SNG.BEST.GUESS")]

write.table(barcode_indiv, "E-ATAC-barcode-indiv.tsv", quote = F, row.names = F, col.names = F, sep = "\t")


barcodes <- read.table("/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/RNA_clustering/insulin_barcodes.txt", header = T, sep = "\t")
barcode_list <- barcodes$barcode
barcode_list <- as.data.frame(barcode_list)
#convert to atac
rna_barcodes <- read.table("/scratch/scjp_root/scjp1/christav/fap_village_multiome/Github/snRNAseq-NextFlow/737-arc-v1.txt")
colnames(rna_barcodes) <- "RNA"
atac_barcodes <- read.table("/scratch/scjp_root/scjp1/christav/fap_village_multiome/Github/snATACseq-NextFlow/737K-arc-v1.txt")
colnames(atac_barcodes) <- "ATAC"
convert <- cbind(atac_barcodes, rna_barcodes)
colnames(barcode_list) <- "RNA"
barcode_list_convert <- merge(barcode_list, convert)
barcode_list_atac <- barcode_list_convert$ATAC
barcode_list_atac <- as.data.frame(barcode_list_atac)
write.table(barcode_list_atac, "/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/RNA_clustering/insulin_atac_barcode_list.txt", col.names = F, row.names = F, quote = F, sep = "\t")
