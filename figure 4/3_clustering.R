setwd("~/links/christav/fap_village_multiome/results/RNA_clustering")
library(Seurat)
library(glue)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sctransform)

rna <- readRDS("RNA_final.rds")
rna_basal <- readRDS("RNA_basal_final.rds")
rna_insulin <- readRDS("RNA_insulin_final.rds")

load_mm <- function(matrix_file, features_file, barcodes_file) {
  tmp <- as(Matrix::readMM(matrix_file), 'dgCMatrix')
  features <- read.table(features_file, as.is=T, sep='\t', head=F)
  features <- paste0(features$V1, ' (', features$V2, ')')
  barcodes <- read.table(barcodes_file, as.is=T, head=F)[,1]
  dimnames(tmp) <- list(features, barcodes)
  return(tmp)
}

RNA_MTX <- file.path('/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/doublet_detection_results5/pass-qc-nuclei-counts-with-doublets/basal.matrix.mtx')
RNA_FEATURES <- file.path('/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/doublet_detection_results5/pass-qc-nuclei-counts-with-doublets/basal.features.tsv')
RNA_BARCODES <- file.path('/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/doublet_detection_results5/pass-qc-nuclei-counts-with-doublets/basal.barcodes.tsv')
mtx_basal <- load_mm(RNA_MTX, RNA_FEATURES, RNA_BARCODES)
#getting pass QC barcodes as singlets
pass_qc_barcodes <- read.table("basal_barcodes.txt", sep="\t", header=T) #file generated below after final QC filtering
#next two lines are for if using demuxlet's output
#pass_qc_singlets <- subset(pass_qc_barcodes, V2 == "SNG")
#pass_qc_barcodes <- pass_qc_singlets["V1"]
pass_qc_barcodes <- pass_qc_barcodes["barcode"]
# subset to our pass QC barcodes
mtx_basal <- mtx_basal[,eval(parse(text=pass_qc_barcodes))]
#add library to barcode
colnames(mtx_basal)<-paste(colnames(mtx_basal),"basal",sep="_")

#for generating just basal sample clustering
rna_basal <- CreateSeuratObject(counts = mtx_basal, min.cells=5)

RNA_MTX <- file.path('/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/doublet_detection_results5/pass-qc-nuclei-counts-with-doublets/insulin.matrix.mtx')
RNA_FEATURES <- file.path('/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/doublet_detection_results5/pass-qc-nuclei-counts-with-doublets/insulin.features.tsv')
RNA_BARCODES <- file.path('/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/doublet_detection_results5/pass-qc-nuclei-counts-with-doublets/insulin.barcodes.tsv')
mtx_insulin <- load_mm(RNA_MTX, RNA_FEATURES, RNA_BARCODES)
#getting pass QC barcodes as singlets
pass_qc_barcodes <- read.table("insulin_barcodes.txt", sep="\t", header=T) #file generated below after final QC filtering
#next two lines are for if using demuxlet's output
#pass_qc_singlets <- subset(pass_qc_barcodes, V2 == "SNG")
#pass_qc_barcodes <- pass_qc_singlets["V1"]
pass_qc_barcodes <- pass_qc_barcodes["barcode"]
# subset to our pass QC barcodes
mtx_insulin <- mtx_insulin[,eval(parse(text=pass_qc_barcodes))]
#add library to barcode
colnames(mtx_insulin)<-paste(colnames(mtx_insulin),"insulin",sep="_")

#for generating just insulin sample clustering
rna_insulin <- CreateSeuratObject(counts = mtx_insulin, min.cells=5)

#for joint clustering, concatenate samples together
mtx.big <- cbind(mtx_insulin, mtx_basal)

rna <- CreateSeuratObject(counts = mtx.big, min.cells=5)

#metadata: environment
enviro <- rownames(rna@meta.data)
enviro <- as.list(enviro)
enviro <- data.frame(enviro)
enviro <- t(enviro)
enviro <- data.frame(enviro)
colnames(enviro) <- c("V1")
rownames(enviro) <- 1:nrow(enviro)
enviro <- enviro %>% separate(V1, c("barcode","enviro"), sep = "_")
rna@meta.data$enviro <- enviro$enviro

#metadata: individual 
insulin <- read.table("../doublet_detection_results5/demuxlet/processed/insulin.assignments.txt", sep="\t", header=F)
insulin <- subset(insulin, V2 == "SNG")
insulin$V1 <- paste(insulin$V1,"insulin",sep="_")

basal <- read.table("../doublet_detection_results5/demuxlet/processed/basal.assignments.txt", sep="\t", header=F)
basal <- subset(basal, V2 == "SNG")
basal$V1 <- paste(basal$V1,"basal",sep="_")

metadata <- rbind(insulin, basal)
metadata <- metadata[,c("V1", "V3")]
colnames(metadata) <- c("barcode", "donor")
metadata_sorted <- metadata[match(colnames(rna), metadata$barcode), ]
rna <- AddMetaData(rna, metadata = metadata_sorted)

#adding individual metadata to basal/insulin clustering
colnames(basal) <- c("barcode", "nuclei", "donor")
basal_sorted <- basal[match(colnames(rna_basal), basal$barcode), ]
rna_basal <- AddMetaData(rna_basal, metadata = basal_sorted)
colnames(insulin) <- c("barcode", "nuclei", "donor")
insulin_sorted <- insulin[match(colnames(rna_insulin), insulin$barcode), ]
rna_insulin <- AddMetaData(rna_insulin, metadata = insulin_sorted)

#metadata: sex, ogtt, bmi, age
indiv_metadata <- read.table("village_metadata.txt", header=T, sep = ",")
colnames(indiv_metadata) <- c("donor", "sex", "bmi", "ogtt", "age")
indiv_all <- merge(metadata_sorted, indiv_metadata, by = "donor")
indiv_all <- indiv_all[match(colnames(rna), indiv_all$barcode), ]
rna <- AddMetaData(rna, metadata = indiv_all)
#adding to basal/insulin clustering
indiv_basal <- merge(basal_sorted, indiv_metadata, by = "donor")
indiv_basal <- indiv_basal[match(colnames(rna_basal), indiv_basal$barcode), ]
rna_basal <- AddMetaData(rna_basal, metadata = indiv_basal)
indiv_insulin <- merge(insulin_sorted, indiv_metadata, by = "donor")
indiv_insulin <- indiv_insulin[match(colnames(rna_insulin), indiv_insulin$barcode), ]
rna_insulin <- AddMetaData(rna_insulin, metadata = indiv_insulin)

#adding qc variables from RNA pipeline
basal_qc <- read.table("/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/RNA_results/qc/basal-hg38.qc.txt", sep="\t", header=T)
basal_qc$barcode <- paste(basal_qc$barcode, "basal", sep = "_")
insulin_qc <- read.table("/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/RNA_results/qc/insulin-hg38.qc.txt", sep="\t", header=T)
insulin_qc$barcode <- paste(insulin_qc$barcode, "insulin", sep = "_")
qc_meta <- rbind(insulin_qc, basal_qc)
qc_meta <- qc_meta[match(colnames(rna), qc_meta$barcode), ]
rna <- AddMetaData(rna, metadata = qc_meta)

#subset to just donors we have information for (shouldn't be necessary but as a sanity check)
rna <- rna[, !is.na(rna$donor)]

rna_basal <- rna_basal[, !is.na(rna_basal$donor)]
rna_insulin <- rna_insulin[, !is.na(rna_insulin$donor)]

#scTransform and Harmony
rna[["RNA"]] <- split(rna[["RNA"]], f = rna$enviro)
rna <- SCTransform(rna, vars.to.regress = c("donor", "age", "sex", "bmi", "ogtt"), verbose = FALSE)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))
options(future.globals.maxSize = 3e+09)
rna <- IntegrateLayers( 
  object = rna,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  verbose = F
)
ElbowPlot(rna) #pick where levels off, range of like 8-20 maybe, 10 for final clustering
PCS <- 10
rna <- FindNeighbors(rna, dims = 1:PCS)
rna <- FindClusters(rna, resolution = 0.5)
rna <- RunUMAP(rna, dims = 1:PCS)
#plots
DimPlot(rna, reduction = "umap", label = TRUE)
DimPlot(rna, reduction = "umap", group.by = "donor")
#removing small side cluster of less than 100 nuclei
#rna <- rna[, rna$seurat_clusters != 7]
#then rerun scTransform, PCA, Harmony, UMAP
saveRDS(rna, file = "RNA_final.rds")

#saving barcodes/donors of finalized clustering, getting plot of donors by sample
barcodes <- colnames(rna)
barcodes <- as.data.frame(barcodes)
donors <- rna$donor
donors <- as.data.frame(donors)
barcodes <- cbind(barcodes, donors)
barcodes$donors <- gsub(".*\\(([^)]+)\\).*", "\\1", barcodes$donors)
barcodes <- barcodes %>% separate(barcodes, c("barcode","enviro"), sep = "_")
basal_barcodes <- subset(barcodes, enviro == "basal")
insulin_barcodes <- subset(barcodes, enviro == "insulin")
write.table(basal_barcodes, "basal_barcodes.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(insulin_barcodes, "insulin_barcodes.txt", row.names = F, col.names = T, quote = F, sep = "\t")
basal_donors <- table(basal_barcodes$donors)
basal_donors <- as.data.frame(basal_donors)
insulin_donors <- table(insulin_barcodes$donors)
insulin_donors <- as.data.frame(insulin_donors)
colnames(basal_donors) <- c("donor", "freq")
basal_donors$enviro <- "basal"
colnames(insulin_donors) <- c("donor", "freq")
insulin_donors$enviro <- "insulin"
both_donors <- rbind(basal_donors, insulin_donors)
ggplot(both_donors, aes(x=reorder(donor, -freq), y=freq, fill = enviro)) +
  ggtitle("Multiome Donors") + xlab("Donors") + ylab("Frequency") +
  geom_bar(stat = "identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#getting all of the UMAPs and QC metrics for clustering
DimPlot(rna, reduction = "umap", group.by = "enviro")
DimPlot(rna, reduction = "umap", group.by = "donor") 
DimPlot(rna, reduction = "umap", group.by = "sex") 
DimPlot(rna, reduction = "umap", group.by = "ogtt") 
FeaturePlot(rna, reduction = "umap", features="age", cols = c("red", "blue"))
FeaturePlot(rna, reduction = "umap", features="bmi", cols = c("red", "blue"))
#plotting qc metrics
FeaturePlot(rna, reduction = "umap", features="umis", cols = c("lightblue", "blue"))
FeaturePlot(rna, reduction = "umap", features="fraction_mitochondrial", cols = c("lightblue", "blue"))
FeaturePlot(rna, reduction = "umap", features="nFeature_RNA", cols = c("lightblue", "blue"))
FeaturePlot(rna, reduction = "umap", features="total_reads", cols = c("lightblue", "blue"))
FeaturePlot(rna, reduction = "umap", features="uniquely_mapped_reads", cols = c("lightblue", "blue"))
VlnPlot(rna, features = c("umis"))
VlnPlot(rna, features = c("fraction_mitochondrial"))
VlnPlot(rna, features = c("nFeature_RNA"))
VlnPlot(rna, features = c("total_reads"))
VlnPlot(rna, features = c("uniquely_mapped_reads"))
#marker genes 
MARKERS <- c("GLI1") 
unlist(lapply(glue('\\({MARKERS}\\)'), function(x){grep(x, rownames(mtx.big), value=T, ignore.case=T)}))
#fap markers
VlnPlot(rna, features = c("ENSG00000135318.12 (NT5E)"))
VlnPlot(rna, features = c("ENSG00000118523.6 (CCN2)"))
VlnPlot(rna, features = c("ENSG00000134363.12 (FST)"))
FeaturePlot(rna, reduction = "umap", features = c("ENSG00000135318.12 (NT5E)", "ENSG00000118523.6 (CCN2)", "ENSG00000134363.12 (FST)"))
#fap subtype markers: fibrogenic
VlnPlot(rna, features = c("ENSG00000196352.15 (CD55)"))
VlnPlot(rna, features = c("ENSG00000075223.14 (SEMA3C)"))
FeaturePlot(rna_insulin, reduction = "umap", features = c("ENSG00000196352.15 (CD55)", "ENSG00000075223.14 (SEMA3C)"))
FeaturePlot(rna, reduction = "umap", features = c("ENSG00000075223.14 (SEMA3C)"))
#fap subtype markers: adipogenic
VlnPlot(rna, features = c("ENSG00000090339.9 (ICAM1)"))
VlnPlot(rna, features = c("ENSG00000117525.14 (F3)"))
FeaturePlot(rna, reduction = "umap", features = c("ENSG00000090339.9 (ICAM1)", "ENSG00000117525.14 (F3)"))
FeaturePlot(rna, reduction = "umap", features = c("ENSG00000117525.14 (F3)"))
#fap subtype markers: tenogenic
VlnPlot(rna, features = c("ENSG00000154864.12 (PIEZO2)")) #not much there
FeaturePlot(rna, reduction = "umap", features = c("ENSG00000154864.12 (PIEZO2)"))
#fap subtype markers: progenitors
VlnPlot(rna, features = c("ENSG00000154096.13 (THY1)"))
FeaturePlot(rna, reduction = "umap", features = c("ENSG00000154096.13 (THY1)"))
#other markers
VlnPlot(rna, features = c("ENSG00000132170.21 (PPARG)"))
VlnPlot(rna, features = c("ENSG00000133392.18 (MYH11)"))
VlnPlot(rna, features = c("ENSG00000081237.20 (PTPRC)"))
FeaturePlot(rna, reduction = "umap", features = c("ENSG00000132170.21 (PPARG)", "ENSG00000133392.18 (MYH11)", "ENSG00000081237.20 (PTPRC)"))
#insulin stimulated genes
VlnPlot(rna, features = c("ENSG00000177606.6 (JUN)"))
VlnPlot(rna, features = c("ENSG00000170345.10 (FOS)"))
VlnPlot(rna, features = c("ENSG00000159069.14 (FBXW5)"))
VlnPlot(rna, features = c("ENSG00000124299.14 (PEPD)"))
VlnPlot(rna, features = c("ENSG00000145431.11 (PDGFC)"))
FeaturePlot(rna, reduction = "umap", features = c("ENSG00000177606.6 (JUN)", "ENSG00000170345.10 (FOS)", "ENSG00000159069.14 (FBXW5)", "ENSG00000124299.14 (PEPD)", "ENSG00000145431.11 (PDGFC)"))

#individual donor
#unique(rna@meta.data$donor)
#[1] "SNG (32102)" "SNG (32032)" "SNG (32057)" "SNG (32012)" "SNG (32120)" "SNG (12181)"
#[7] "SNG (32028)" "SNG (12171)" "SNG (12039)" "SNG (12066)" "SNG (12180)" "SNG (32015)"
#[13] "SNG (32100)" "SNG (12040)" "SNG (12145)" "SNG (12034)" "SNG (12182)" "SNG (12150)"
#[19] "SNG (32019)" "SNG (32049)" "SNG (12031)" "SNG (32086)" "SNG (32044)" "SNG (32051)"
#[25] "SNG (32031)" "SNG (12006)" "SNG (32021)" "SNG (12134)" "SNG (22025)"
cells_to_keep <- rna@meta.data$donor == "SNG (22025)"
rna_filtered <- rna[, cells_to_keep]
DimPlot(rna_filtered, reduction = "umap", group.by = "donor")

#adding how long lines were thawed for to the metadata/FACS expression of FAP marker and plotting
days_thawed <- read.table("days_thawed.txt", sep = "\t", header = T)
barcodes <- colnames(rna)
barcodes <- as.data.frame(barcodes)
donors <- rna$donor
donors <- as.data.frame(donors)
barcode_donor <- cbind(barcodes, donors)
colnames(barcode_donor) <- c("barcode", "donor")
barcode_donor_days <- merge(barcode_donor, days_thawed, by = "donor")
barcode_donor_days$donor <- NULL
barcode_donor_days <- barcode_donor_days[match(colnames(rna), barcode_donor_days$barcode), ]
rna <- AddMetaData(rna, metadata = barcode_donor_days)
FeaturePlot(rna, reduction = "umap", features="Days_Thawed", cols = c("red", "blue"))
FeaturePlot(rna, reduction = "umap", features="FACS", cols = c("red", "blue"))

#clustering just basal or insulin
rna_basal <- SCTransform(rna_basal, vars.to.regress = c("donor", "age", "sex", "bmi", "ogtt"), verbose = FALSE)
rna_basal <- RunPCA(rna_basal, features = VariableFeatures(object = rna_basal))
ElbowPlot(rna_basal) #pick where levels off, 9 before donor specific
PCS <- 9
rna_basal <- FindNeighbors(rna_basal, dims = 1:PCS)
rna_basal <- FindClusters(rna_basal, resolution = 0.5)
rna_basal <- RunUMAP(rna_basal, dims = 1:PCS)
#plots
DimPlot(rna_basal, reduction = "umap", label = TRUE)
DimPlot(rna_basal, reduction = "umap", group.by = "donor")
saveRDS(rna_basal, file = "RNA_basal_final.rds")

rna_insulin <- SCTransform(rna_insulin, vars.to.regress = c("donor", "age", "sex", "bmi", "ogtt"), verbose = FALSE)
rna_insulin <- RunPCA(rna_insulin, features = VariableFeatures(object = rna_insulin))
ElbowPlot(rna_insulin) #pick where levels off, 11 before donor specific
PCS <- 11
rna_insulin <- FindNeighbors(rna_insulin, dims = 1:PCS)
rna_insulin <- FindClusters(rna_insulin, resolution = 0.5)
rna_insulin <- RunUMAP(rna_insulin, dims = 1:PCS)
#plots
DimPlot(rna_insulin, reduction = "umap", label = TRUE)
DimPlot(rna_insulin, reduction = "umap", group.by = "donor")
saveRDS(rna_insulin, file = "RNA_insulin_final.rds")



#using clustering of individual samples to mark where those are on joint clustering, since the subclusters are clearer on the single samples
rna <- readRDS("RNA_final.rds")
rna_basal <- readRDS("RNA_basal_final.rds")
rna_insulin <- readRDS("RNA_insulin_final.rds")

DimPlot(rna, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(rna_basal, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(rna_insulin, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#getting barcodes of basal fibrogenic/progenitors
Idents(rna_basal) <- 'seurat_clusters'
barcodes <- CellsByIdentities(rna_basal, idents = 3)
barcodes2 <- CellsByIdentities(rna_basal, idents = 2)
barcodes <- unlist(barcodes)
barcodes2 <- unlist(barcodes2)
barcodes <- c(barcodes, barcodes2)
rna$basal_fibrogenic <- colnames(rna) %in% barcodes
table(rna$basal_fibrogenic)
DimPlot(rna, reduction = "umap", group.by = "basal_fibrogenic")

#Idents(rna_basal) <- 'seurat_clusters'
barcodes <- CellsByIdentities(rna_basal, idents = 0)
barcodes2 <- CellsByIdentities(rna_basal, idents = 1)
barcodes3 <- CellsByIdentities(rna_basal, idents = 6)
barcodes4 <- CellsByIdentities(rna_basal, idents = 4)
barcodes5 <- CellsByIdentities(rna_basal, idents = 5)
barcodes6 <- CellsByIdentities(rna_basal, idents = 7)
barcodes <- unlist(barcodes)
barcodes2 <- unlist(barcodes2)
barcodes3 <- unlist(barcodes3)
barcodes4 <- unlist(barcodes4)
barcodes5 <- unlist(barcodes5)
barcodes6 <- unlist(barcodes6)
barcodes <- c(barcodes, barcodes2, barcodes3, barcodes4, barcodes5, barcodes6)
rna$basal_progenitors <- colnames(rna) %in% barcodes
table(rna$basal_progenitors)
DimPlot(rna, reduction = "umap", group.by = "basal_progenitors")

#barcodes of insulin adipogenic/fibrogenic/progenitors
rna_insulin <- FindClusters(rna_insulin, resolution = 1) #increasing the clusters to get the nuclei I want from trajectory analysis
PCS = 11
rna_insulin <- RunUMAP(rna_insulin, dims = 1:PCS)
DimPlot(rna_insulin, reduction = "umap", label = TRUE) + NoLegend()

Idents(rna_insulin) <- 'seurat_clusters'
barcodes <- CellsByIdentities(rna_insulin, idents = 0)
barcodes <- unlist(barcodes)
rna$insulin_adipogenic <- colnames(rna) %in% barcodes
table(rna$insulin_adipogenic)
DimPlot(rna, reduction = "umap", group.by = "insulin_adipogenic")

#Idents(rna_insulin) <- 'seurat_clusters'
barcodes <- CellsByIdentities(rna_insulin, idents = 1)
barcodes2 <- CellsByIdentities(rna_insulin, idents = 8)
barcodes <- unlist(barcodes)
barcodes2 <- unlist(barcodes2)
barcodes <- c(barcodes, barcodes2)
rna$insulin_fibrogenic <- colnames(rna) %in% barcodes
table(rna$insulin_fibrogenic)
DimPlot(rna, reduction = "umap", group.by = "insulin_fibrogenic")

#Idents(rna_insulin) <- 'seurat_clusters'
barcodes <- CellsByIdentities(rna_insulin, idents = 7)
barcodes2 <- CellsByIdentities(rna_insulin, idents = 6)
barcodes3 <- CellsByIdentities(rna_insulin, idents = 2)
barcodes <- unlist(barcodes)
barcodes2 <- unlist(barcodes2)
barcodes3 <- unlist(barcodes3)
barcodes <- c(barcodes, barcodes2, barcodes3)
rna$insulin_progenitors <- colnames(rna) %in% barcodes
table(rna$insulin_progenitors)
DimPlot(rna, reduction = "umap", group.by = "insulin_progenitors")

#label based on this
rna <- FindClusters(rna, resolution = 1.8) #increasing the clusters to the nuclei I want from trajectory analysis
PCS = 10
rna <- RunUMAP(rna, dims = 1:PCS)
DimPlot(rna, reduction = "umap", label = TRUE) + NoLegend()

new.cluster.ids <- c("Adipogenic", "Fibrogenic", "Adipogenic", "Progenitors", "Fibrogenic", "Fibrogenic", "Progenitors", "Progenitors", "Adipogenic", "Progenitors", "Fibrogenic", "Progenitors", "Fibrogenic", "Progenitors", "Progenitors", "Progenitors", "Progenitors", "Fibrogenic")
names(new.cluster.ids) <- levels(rna)
rna <- RenameIdents(rna, new.cluster.ids)
rna$subcluster <- Idents(rna)
DimPlot(rna, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(rna, reduction = "umap", label = FALSE, pt.size = 0.5) 

subclusters <- table(rna$seurat_clusters, rna$enviro)
#basal insulin
#0   577     590
#1   345     811
#2   261     719
#3   278     189
#4   110     355
#5   222     175
#6    96     122
#7    49      49
#8    67      26
#9    34      51
subclusters <- as.data.frame(subclusters)
subclusters$new.cluster.ids <- c("Adipogenic", "Fibrogenic", "Adipogenic", "Progenitors", "Fibrogenic", "Fibrogenic", "Progenitors", "Progenitors", "Adipogenic", "Progenitors", "Fibrogenic", "Progenitors", "Fibrogenic", "Progenitors", "Progenitors", "Progenitors", "Progenitors", "Fibrogenic")
subclusters
adipo <- subset(subclusters, new.cluster.ids == "Adipogenic")
fibro <- subset(subclusters, new.cluster.ids == "Fibrogenic")
prog <- subset(subclusters, new.cluster.ids == "Progenitors")

#fisher's exact test for subcluster composition
adipogenic <- subclusters[subclusters$new.cluster.ids == "Adipogenic", ]
not_adipogenic <- subclusters[subclusters$new.cluster.ids != "Adipogenic", ]

A <- sum(adipogenic$Freq[adipogenic$Var2 == "basal"])    # basal in Adipogenic
C <- sum(adipogenic$Freq[adipogenic$Var2 == "insulin"])  # insulin in Adipogenic
B <- sum(not_adipogenic$Freq[not_adipogenic$Var2 == "basal"])    # basal in NOT Adipogenic
D <- sum(not_adipogenic$Freq[not_adipogenic$Var2 == "insulin"])  # insulin in NOT Adipogenic

contingency_table <- matrix(c(A, C, B, D), nrow = 2, byrow = TRUE,
                            dimnames = list(c("Adipogenic", "Not_Adipogenic"),
                                            c("Basal", "Insulin")))

fisher.test(contingency_table)
#	Fisher's Exact Test for Count Data
#data:  contingency_table
#p-value = 1.024e-11
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  0.5440624 0.7189159
#sample estimates:
#  odds ratio 
#0.6257545 

fibrogenic <- subclusters[subclusters$new.cluster.ids == "Fibrogenic", ]
not_fibrogenic <- subclusters[subclusters$new.cluster.ids != "Fibrogenic", ]

A <- sum(fibrogenic$Freq[fibrogenic$Var2 == "basal"])    # basal in fibrogenic
C <- sum(fibrogenic$Freq[fibrogenic$Var2 == "insulin"])  # insulin in fibrogenic
B <- sum(not_fibrogenic$Freq[not_fibrogenic$Var2 == "basal"])    # basal in NOT fibrogenic
D <- sum(not_fibrogenic$Freq[not_fibrogenic$Var2 == "insulin"])  # insulin in NOT fibrogenic

contingency_table <- matrix(c(A, C, B, D), nrow = 2, byrow = TRUE,
                            dimnames = list(c("fibrogenic", "Not_fibrogenic"),
                                            c("Basal", "Insulin")))

fisher.test(contingency_table)
#		Fisher's Exact Test for Count Data
#data:  contingency_table
#p-value = 0.02823
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.013754 1.290821
#sample estimates:
#  odds ratio 
#1.144008 

progenitors <- subclusters[subclusters$new.cluster.ids == "Progenitors", ]
not_progenitors <- subclusters[subclusters$new.cluster.ids != "Progenitors", ]

A <- sum(progenitors$Freq[progenitors$Var2 == "basal"])    # basal in progenitors
C <- sum(progenitors$Freq[progenitors$Var2 == "insulin"])  # insulin in progenitors
B <- sum(not_progenitors$Freq[not_progenitors$Var2 == "basal"])    # basal in NOT progenitors
D <- sum(not_progenitors$Freq[not_progenitors$Var2 == "insulin"])  # insulin in NOT progenitors

contingency_table <- matrix(c(A, C, B, D), nrow = 2, byrow = TRUE,
                            dimnames = list(c("progenitors", "Not_progenitors"),
                                            c("Basal", "Insulin")))

fisher.test(contingency_table)
#	Fisher's Exact Test for Count Data
#data:  contingency_table
#p-value = 0.000205
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.106986 1.396084
#sample estimates:
#  odds ratio 
#1.243177 
