setwd("~/links/christav/fall22_multiome/scripts")

peaks <- read.table("../macs2/22042.merged.broadPeak")
peaks$V4 <- paste("peak_", row.names(peaks))
write.table(peaks, "../macs2/22042.label.broadPeak", col.names = F, row.names = F, quote = F, sep = "\t")
