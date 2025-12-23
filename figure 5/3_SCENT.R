#SCENT

#getting the rna/atac matrices
setwd("~/links/christav/fap_village_multiome/results/RNA_clustering")
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(glue)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sctransform)
library(harmony)
library(SCENT)
library(Matrix)

rna <- readRDS("ATAC_RNA_seuratobj.rds")

#getting rna input
rna_insulin <- GetAssayData(
  object = rna,
  assay = "RNA",
  layer = "counts.insulin"
)

rna_basal <- GetAssayData(
  object = rna,
  assay = "RNA",
  layer = "counts.basal"
)
#check that gene names same and cells different
#stopifnot(
#  identical(rownames(rna_insulin), rownames(rna_basal)),
#  length(intersect(colnames(rna_insulin), colnames(rna_basal))) == 0
#)

rna_data <- cbind(rna_insulin, rna_basal)

#getting atac input
atac_data <- GetAssayData(
  object = rna,
  assay = "ATAC",
  layer = "counts"
)
#check atac/rna overlap and check format
#common_cells <- intersect(colnames(rna_data), colnames(atac_data))

#rna_data  <- rna_data[, common_cells]
#atac_data <- atac_data[, common_cells]

#stopifnot(
#  inherits(rna_data, "dgCMatrix"),
#  inherits(atac_data, "dgCMatrix"),
#  identical(colnames(rna_data), colnames(atac_data))
#)

#metadata
meta.data <- rna@meta.data
meta.data$logUMI <- log(meta.data$umis)
meta.data$celltype <- "iPSC-FAP"
meta.data_scent <- meta.data[, c("barcode", "fraction_mitochondrial", "celltype", "logUMI", "orig.ident", "seurat_clusters", "donor", "sex", "bmi", "ogtt", "age")]
colnames(meta.data_scent) <- c("cell", "fraction_mitochondrial", "celltype", "logUMI", "orig.ident", "seurat_clusters", "donor", "sex", "bmi", "ogtt", "age")

#stopifnot(
#  all(meta.data_scent$cell == colnames(rna_data)),
#  nrow(meta.data_scent) == ncol(rna_data)
#)

#peak info (testing near my interesting peak)
#MARKERS <- c("SIPA1") 
#unlist(lapply(glue('\\({MARKERS}\\)'), function(x){grep(x, rownames(rna), value=T, ignore.case=T)}))

genes <- c("ENSG00000142186.17 (SCYL1)", "ENSG00000168056.15 (LTBP3)", "ENSG00000176973.8 (FAM89B)", "ENSG00000173442.13 (EHBP1L1)", "ENSG00000173327.8 (MAP3K11)")
genes <- data.frame(genes)
genes$V2 <- "chr11-65559734-65561243"
#peaks <- rownames(atac_data)
#peaks_chr11 <- peaks[grep("^chr11-", peaks)]
#peaks_chr11 <- data.frame(peaks_chr11)
#peak_to_check <- "chr11-65559734-65561243"
#peak_to_check %in% rownames(atac_data)
colnames(genes) <- c("V1", "V2")

#SCENT
SCENT_obj <- CreateSCENTObj(rna = rna_data, atac = atac_data, meta.data = meta.data_scent,
                            peak.info = genes,
                            covariates = c("logUMI","fraction_mitochondrial","orig.ident", "donor", "sex", "bmi", "ogtt", "age"), 
                            celltypes = "celltype")

SCENT_obj <- SCENT_algorithm(object = SCENT_obj, celltype = "iPSC-FAP", ncores = 1, regr = 'poisson', bin = TRUE)

result <- SCENT_obj@SCENT.result
write.table(result, "scent_result.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
#LTBP3 is significant! 
