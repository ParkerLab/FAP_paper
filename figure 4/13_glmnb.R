setwd("/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/differential_peaks/glmnb/")
suppressPackageStartupMessages({library(SingleCellExperiment)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(glue)
  library(cowplot)
  library(pbmcapply)
  library(DHARMa)
  library(magrittr)
  library(glmmTMB)
  library(pscl)
  library(optparse)
  library(lmtest)
})

#from snapATAC2, also see peak_mat_R_obj.R for the second one
merged_peaks <- read.table("../merged_peaks.tsv", sep = "\t", header = T)
load("peak_mat.Rda")

#prepping metadata
#metadata: individual
#insulin <- read.table("../../doublet_detection_results5/demuxlet/processed/insulin.assignments.txt", sep="\t", header=F)
#insulin <- subset(insulin, V2 == "SNG")
#insulin$V1 <- paste(insulin$V1,"insulin",sep="_")
#basal <- read.table("../../doublet_detection_results5/demuxlet/processed/basal.assignments.txt", sep="\t", header=F)
#basal <- subset(basal, V2 == "SNG")
#basal$V1 <- paste(basal$V1,"basal",sep="_")
#metadata <- rbind(insulin, basal)
#metadata <- metadata[,c("V1", "V3")]
#colnames(metadata) <- c("barcode", "donor")
#metadata: sex, ogtt, bmi, age
#indiv_metadata <- read.table("/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/RNA_clustering/village_metadata.txt", header=T, sep = ",")
#colnames(indiv_metadata) <- c("donor", "sex", "bmi", "ogtt", "age")
#indiv_all <- merge(metadata, indiv_metadata, by = "donor")
#indiv_all <- indiv_all %>% separate(barcode, c("rna", "enviro"), sep = "_")
#indiv_all$barcode_sample <- paste(indiv_all$rna, "_", indiv_all$enviro, sep="")
#ATAC_barcodes <- read.table("/scratch/scjp_root/scjp1/christav/fap_village_multiome/Github/snATACseq-NextFlow/737K-arc-v1.txt", col.names = "ATAC")
#RNA_barcodes <- read.table("/scratch/scjp_root/scjp1/christav/fap_village_multiome/Github/snRNAseq-NextFlow/737-arc-v1.txt", col.names = "RNA")
#ATAC_barcodes$rna <- RNA_barcodes$RNA
#indiv_all <- merge(indiv_all, ATAC_barcodes)
#indiv_all$atac_barcode_sample <- paste(indiv_all$ATAC, "_", indiv_all$enviro, sep="")
#add TSS enrichment
#ataqv <- read.table("/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/ATAC_results/ataqv/single-nucleus/basal-hg38.txt", sep = "\t", header = T)
#ataqv$atac_barcode_sample <- paste(ataqv$name, "_basal", sep = "")
#ataqv_basal <- ataqv[,c("atac_barcode_sample", "tss_enrichment")]
#ataqv <- read.table("/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/ATAC_results/ataqv/single-nucleus/insulin-hg38.txt", sep = "\t", header = T)
#ataqv$atac_barcode_sample <- paste(ataqv$name, "_insulin", sep = "")
#ataqv_insulin <- ataqv[,c("atac_barcode_sample", "tss_enrichment")]
#ataqv_both <- rbind(ataqv_basal, ataqv_insulin)
#metadata <- merge(indiv_all, ataqv_both)
#subset to just final nuclei
#basal_passqc <- read.table("../../RNA_clustering/basal_atac_barcode_list.txt", header = F, sep = "\t")
#basal_passqc$atac_barcode_sample <- paste(basal_passqc$V1, "_basal", sep = "")
#insulin_passqc <- read.table("../../RNA_clustering/insulin_atac_barcode_list.txt", header = F, sep = "\t")
#insulin_passqc$atac_barcode_sample <- paste(insulin_passqc$V1, "_insulin", sep = "")
#passqc <- rbind(basal_passqc, insulin_passqc)
#metadata <- metadata %>%
#  mutate(passqc = if_else(metadata$atac_barcode_sample %in% passqc$atac_barcode_sample,"Y","N"))
#metadata <- subset(metadata, passqc == "Y")
#write.table(metadata, "metadata.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep=",")

#load metadata
metadata <- read.table("metadata.txt", sep = ",", header = T)
metadata$atac_barcode_sample <- paste(metadata$enviro, ":", metadata$ATAC, sep = "")

#following Arushi's script from here: "/gpfs/accounts/scjp_root/scjp99/arushiv/pfizer_hfpef/atac/glmmnb/scripts/bin/glmnb.R"
expr = assay(pads)
meta =  droplevels(colData(pads))
meta$atac_barcode_sample <- rownames(meta)
meta <- merge(meta, metadata)
meta$sex = as.factor(meta$sex)
meta$ogtt = as.factor(meta$ogtt)
meta$donor = as.factor(meta$donor)
meta$age = as.numeric(meta$age)
meta$tss_enrichment = as.numeric(meta$tss_enrichment)

fmla <- as.formula("GENE ~ sample + sex + age + ogtt + bmi + tss_enrichment + offset(log(total_counts)) + (1|donor)")
null_fmla <- as.formula("GENE ~ sex + age + ogtt + bmi + tss_enrichment + offset(log(total_counts)) + (1|donor)")

glmnb = function(index, fmla, null_fmla, fmla_random){
  test_dat = data.frame(
    GENE = expr[index,],
    sample = meta$sample,
    sex = as.factor(meta$sex),
    age = meta$age,
    ogtt = as.factor(meta$ogtt),
    bmi = meta$bmi,
    donor = as.factor(meta$donor),
    tss_enrichment = meta$tss_enrichment
  )
  test_dat %<>% mutate(total_counts = as.numeric(colSums(expr)+1))
  
  ## initialize with NA in case of error late
  hfpf = estimate = se = z = p_mod1 = p_val = NA
  
  tryCatch({
    mod1 = glmmTMB(fmla, test_dat, family = nbinom2, REML = FALSE)
    mod_null = glmmTMB(null_fmla, test_dat, family = nbinom2, REML = FALSE)
    hfpf_index = grep("sample", rownames(summary(mod1)$coefficients$cond))
    hfpf = rownames(summary(mod1)$coefficients$cond)[hfpf_index]
    estimate = summary(mod1)$coefficients$cond[hfpf,'Estimate']
    se = summary(mod1)$coefficients$cond[hfpf,'Std. Error']
    z = summary(mod1)$coefficients$cond[hfpf,'z value']
    p_mod1 = summary(mod1)$coefficients$cond[hfpf,'Pr(>|z|)']
    
    print(summary(mod1))
    
    lr_result = lrtest(mod1, mod_null)
    p_val <- lr_result$`Pr(>Chisq)`[2] 
    
  }, error = function(e) {print(e)})
  result = list(index, hfpf, estimate, se, z, p_mod1, p_val) 
  print(result)
  return(result)
  
}

n_chunks = 50
chunk_number = 1 #change this number to run a different batch, or change above line to change number of batches
npeaks = length(rownames(expr))
len_chunk = round(npeaks/n_chunks)
start = (chunk_number - 1) * (len_chunk) + 1
end = start + len_chunk
if (chunk_number == n_chunks){
  end = npeaks
}

set.seed(1234)
peaks_to_run = rownames(expr)[start:end]
#print(glue("selecting {n_chunks} out of {npeaks}")) 
out = pbmclapply(peaks_to_run, mc.cores = 1, function(i){glmnb(i, fmla, null_fmla, fmla_random)})

#out.trim <- out$value #only need this for interactive session

res = do.call(rbind, lapply(out, function(i){as.data.frame(i, col.names = c("index", "sample", "estimate", "se", "z", "p_GLMMNB", "p_LRT"))}))

write.table(res, sep='\t', quote = FALSE, file = "output/FAP_1.tsv")
