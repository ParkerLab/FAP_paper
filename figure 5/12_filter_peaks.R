setwd("~/links/christav/fap_village_multiome/results/differential_peaks/glmnb/scripts")

library(tidyr)

#applying FDR correction and subsetting peaks from there 
all_peaks <- read.table("../output/FAP_all.tsv", header = T, row.names = NULL)
converged_peaks <- subset(all_peaks, p_LRT != "NA")
#FDR correction using Storey qvalue
library(qvalue)
pvalues <- converged_peaks$p_LRT
qobj <- qvalue(p = pvalues, fdr.level = 0.05)
converged_peaks$qvalues <- qobj$qvalues
converged_peaks$pi0 <- qobj$pi0
converged_peaks$lfdr <- qobj$lfdr
converged_peaks$significant <- qobj$significant

library(tidyr)
sig_peaks <- subset(converged_peaks, significant == TRUE)
sig_peaks_bed <- sig_peaks %>% separate(index, c("chr", "coord"), sep = ":")
sig_peaks_bed <- sig_peaks_bed %>% separate(coord, c("start", "end"), sep = "-")
sig_peaks_bed <- sig_peaks_bed[,c("chr", "start", "end")]
non_sig_peaks <- subset(converged_peaks, significant == FALSE)
non_sig_peaks_bed <- non_sig_peaks %>% separate(index, c("chr", "coord"), sep = ":")
non_sig_peaks_bed <- non_sig_peaks_bed %>% separate(coord, c("start", "end"), sep = "-")
non_sig_peaks_bed <- non_sig_peaks_bed[,c("chr", "start", "end")]

converged_peaks_bed <- converged_peaks %>% separate(index, c("chr", "coord"), sep = ":")
converged_peaks_bed <- converged_peaks_bed %>% separate(coord, c("start", "end"), sep = "-")
converged_peaks_bed <- converged_peaks_bed[,c("chr", "start", "end")]

write.table(sig_peaks, file = "../output/sig_peaks/FAP_sig_0.05.tsv", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(sig_peaks_bed, file = "../output/sig_peaks/FAP_sig_0.05.bed", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(non_sig_peaks, file = "../output/not_sig_peaks/FAP_sig_0.05.tsv", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(non_sig_peaks_bed, file = "../output/not_sig_peaks/FAP_sig_0.05.bed", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(converged_peaks, file = "../output/FAP_converged.tsv", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE) 
write.table(converged_peaks_bed, file = "../output/FAP_converged.bed", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE) 
