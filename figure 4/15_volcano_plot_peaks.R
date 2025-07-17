#volcano plot of differential peaks

setwd("~/links/christav/fap_village_multiome/results/differential_peaks/glmnb/scripts")

peaks <- read.table("../output/FAP_converged.tsv", header = F, row.names = NULL)
colnames(peaks) <- c("row", "index", "sample", "estimate", "se", "z", "p_GLMMNB", "p_LRT", "qvalues", "pi0", "lfdr", "significant")
#converged_peaks <- subset(peaks, p_GLMMNB != "NA")

#making volcano plot!
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

#add diffexpression
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
peaks$diffexpressed <- "NO"
peaks$diffexpressed[peaks$estimate > 0 & peaks$significant == "TRUE"] <- "UP"
peaks$diffexpressed[peaks$estimate < 0 & peaks$significant == "TRUE"] <- "DOWN"
table(peaks$diffexpressed)

#plot!
range(peaks$estimate)

#work of art plot
ggplot(data = peaks, aes(x = estimate, y = -log10(p_LRT), col = diffexpressed)) +
  #geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 20), xlim = c(-3, 4)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Significance', #legend_title, 
       x = expression("Estimate"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('Differential Peaks in iPSC-FAPs in Insulin Conditions Compared to Basal') 
