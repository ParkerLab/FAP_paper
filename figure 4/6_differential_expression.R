setwd("~/links/christav/fap_village_multiome/results/mast")
library(Seurat)
library(glue)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sctransform)
library(MAST)
library(scater)

rna <- readRDS("../RNA_clustering/RNA_final.rds")

#get counts matrix for QC volcano plot
counts <- rna[["RNA"]]$counts
counts <- as.data.frame(counts)
zeronuc <- rowSums(counts==0)
totnuc <- ncol(counts)
zerop <- zeronuc/totnuc
zerop <- data.frame(zerop) %>% arrange(rownames(.))

#mast
DefaultAssay(rna) <- "RNA"
rna <- JoinLayers(rna)
DefaultAssay(rna) <- "SCT"
rna <- as.SingleCellExperiment(rna)
#cdr
rna$CDR <- as.numeric(colSums(counts(rna)>0)/nrow(counts(rna)))

model <- "~ enviro + sex + bmi+ ogtt + age + scaled_gene_num + CDR + (1|donor)"

counts(rna) <- calculateTPM(rna)
logcounts(rna) <- log2(counts(rna) + 1)
rna <- SceToSingleCellAssay(rna)

#finding genes
rna <- filterLowExpressedGenes(rna, threshold = 0.2)
rna$enviro <- factor(rna$enviro, levels = unique(rna$enviro))
rna$donor <- factor(rna$donor, levels = unique(rna$donor))
rna$sex <- factor(rna$sex, levels = unique(rna$sex))
rna$bmi <- as.numeric(rna$bmi)
rna$ogtt <- factor(rna$ogtt, levels = unique(rna$ogtt))
rna$age <- as.numeric(rna$age)

cdr <- colSums(assay(rna)>0)
colData(rna)$scaled_gene_num <- scale(cdr)

cond <- factor(colData(rna)$enviro)
cond <- relevel(cond, ref = "basal")
colData(rna)$enviro <- cond

zlmCond.glmer <- zlm(as.formula(model), rna, method = 'glmer', parallel = TRUE,
                     ebayes=FALSE, exprs_value = 'logcounts',
                     strictConvergence = FALSE, fitArgsD = list(nAGQ = 0))

contrasts <- unique(summary(zlmCond.glmer)$datatable$contrast)
test.contrast_abbrev <- as.character(contrasts[grepl("enviro", contrasts)])
summaryCond.glmer <- summary(zlmCond.glmer, doLRT = test.contrast_abbrev)
summaryDt.glmer <- summaryCond.glmer$datatable
summaryDt <- summaryDt.glmer

fcHurdle <- merge(summaryDt[contrast==test.contrast_abbrev & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast==test.contrast_abbrev & component=='logFC', .(primerid, coef, ci.hi, ci.lo, z)], #logFC coefficients
                    by = 'primerid')

#also looking at discrete and continuous coefficients as a QC measure
fcHurdle.d <- merge(summaryDt[contrast==test.contrast_abbrev & component=='D',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast==test.contrast_abbrev & component=='logFC', .(primerid, coef, ci.hi, ci.lo, z)], #logFC coefficients
                  by = 'primerid')
fcHurdle.c <- merge(summaryDt[contrast==test.contrast_abbrev & component=='C',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast==test.contrast_abbrev & component=='logFC', .(primerid, coef, ci.hi, ci.lo, z)], #logFC coefficients
                    by = 'primerid')

all.fcHurdle <- fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
all.fcHurdle$lfcSE <- abs(all.fcHurdle$ci.hi)
all.fcHurdle <- rename(all.fcHurdle, "primerid" = "Gene", `Pr(>Chisq)` = "pvalue", "coef" = "log2FC")

#discrete and continuous
all.fcHurdle.d <- fcHurdle.d[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
all.fcHurdle.d$lfcSE <- abs(all.fcHurdle.d$ci.hi)
all.fcHurdle.d <- rename(all.fcHurdle.d, "primerid" = "Gene", `Pr(>Chisq)` = "d.pvalue", "coef" = "d.log2FC")
all.fcHurdle.c <- fcHurdle.c[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
all.fcHurdle.c$lfcSE <- abs(all.fcHurdle.c$ci.hi)
all.fcHurdle.c <- rename(all.fcHurdle.c, "primerid" = "Gene", `Pr(>Chisq)` = "c.pvalue", "coef" = "c.log2FC")

#merging all 3 together
all.fcHurdle.e <- merge(all.fcHurdle, all.fcHurdle.c, by = "Gene")
all.fcHurdle.f <- merge(all.fcHurdle.e, all.fcHurdle.d, by = "Gene")

#making volcano plot!
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

#add diffexpression
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
all.fcHurdle$diffexpressed <- "NO"
all.fcHurdle$diffexpressed[all.fcHurdle$log2FC > 0.6 & all.fcHurdle$pvalue < 0.05] <- "UP"
all.fcHurdle$diffexpressed[all.fcHurdle$log2FC < -0.6 & all.fcHurdle$pvalue < 0.05] <- "DOWN"

table(all.fcHurdle$diffexpressed)
#DOWN   NO   UP 
#613 7189  329 

#plot!

# Create a basic volcano plot to get xlim/ylim
#ggplot(data = all.fcHurdle, aes(x = log2FC, y = -log10(pvalue))) +
#  geom_point()

#work of art plot, Figure 4H
ggplot(data = all.fcHurdle, aes(x = log2FC, y = -log10(pvalue), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 350), xlim = c(-4, 4)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Significance', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('iPSC-FAPs in Insulin Conditions Compared to Basal') 

# volcano of continuous coefficient
all.fcHurdle.c$diffexpressed <- "NO"
all.fcHurdle.c$diffexpressed[all.fcHurdle.c$c.log2FC > 0.6 & all.fcHurdle.c$c.pvalue < 0.05] <- "UP"
all.fcHurdle.c$diffexpressed[all.fcHurdle.c$c.log2FC < -0.6 & all.fcHurdle.c$c.pvalue < 0.05] <- "DOWN"
ggplot(data = all.fcHurdle.c, aes(x = c.log2FC, y = -log10(c.pvalue), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 350), xlim = c(-4, 4)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Cont Significance', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('iPSC-FAPs in Insulin Conditions Compared to Basal') 

table(all.fcHurdle.c$diffexpressed)
#DOWN   NO   UP 
#613 7239  279 

#volcano of discrete coefficient
all.fcHurdle.d$diffexpressed <- "NO"
all.fcHurdle.d$diffexpressed[all.fcHurdle.d$d.log2FC > 0.6 & all.fcHurdle.d$d.pvalue < 0.05] <- "UP"
all.fcHurdle.d$diffexpressed[all.fcHurdle.d$d.log2FC < -0.6 & all.fcHurdle.d$d.pvalue < 0.05] <- "DOWN"
ggplot(data = all.fcHurdle.d, aes(x = d.log2FC, y = -log10(d.pvalue), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 350), xlim = c(-4, 4)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Discrete Significance', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('iPSC-FAPs in Insulin Conditions Compared to Basal') 

table(all.fcHurdle.d$diffexpressed)
#DOWN   NO   UP 
#581 7222  328 

#colored by D/C logFC-values
ggplot(data = all.fcHurdle.f, aes(x = log2FC, y = -log10(pvalue), col = c.log2FC)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2, alpha = 0.5) + 
  coord_cartesian(ylim = c(0, 350), xlim = c(-5, 5)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Cont log2FC', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('iPSC-FAPs in Insulin Conditions Compared to Basal') 

ggplot(data = all.fcHurdle.f, aes(x = log2FC, y = -log10(pvalue), col = d.log2FC)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2, alpha = 0.5) + 
  coord_cartesian(ylim = c(0, 350), xlim = c(-5, 5)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Discrete log2FC', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('iPSC-FAPs in Insulin Conditions Compared to Basal') 

#volcano plot colored by zero nuclei
all.fcHurdle <- read.table("all.fcHurdle.txt", header = T, sep = "\t")
zerop$Gene <- rownames(zerop)

fchurdle2 <- merge(all.fcHurdle, zerop, by = "Gene")

ggplot(data = fchurdle2, aes(x = log2FC, y = -log10(pvalue), col = zerop)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2, alpha = 0.5) + 
  coord_cartesian(ylim = c(0, 350), xlim = c(-5, 5)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Ratio 0 Counts Nuclei to Total', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('iPSC-FAPs in Insulin Conditions Compared to Basal') 
