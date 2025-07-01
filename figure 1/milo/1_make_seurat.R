#This describes the first step for the milo figures in figure 1 and supp fig 1, forming a Seurat object with the FUSION biopsy data and saving it.

setwd("/scratch/scjp_root/scjp1/christav/milo/Arushi_data")

library(Seurat)
library(miloR)

load_mm <- function(matrix_file, features_file, barcodes_file) {
  tmp <- as(Matrix::readMM(matrix_file), 'dgCMatrix')
  features <- read.table(features_file, as.is=T, sep='\t', head=F)
  features <- paste0(features$V1, ' (', features$V2, ')')
  barcodes <- read.table(barcodes_file, as.is=T, head=F)[,1]
  dimnames(tmp) <- list(features, barcodes)
  return(tmp)
}

folders <- read.table("/scratch/scjp_root/scjp1/christav/milo/get_seurat_obj/folders.txt")
folders <- as.list(folders$V1)

# Initialize an empty list to hold matrices
matrix_list <- list()

# Loop through each folder to load the data
for (folder in folders) {
  RNA_MTX <- file.path(folder, 'matrix.mtx')
  RNA_FEATURES <- file.path(folder, 'features.tsv')
  RNA_BARCODES <- file.path(folder, 'barcodes.tsv')
  
  # Load the matrix and add to the list
  folder_matrix <- load_mm(RNA_MTX, RNA_FEATURES, RNA_BARCODES)
  
  # Append folder name to column names for uniqueness
  colnames(folder_matrix) <- paste(colnames(folder_matrix), folder, sep = "_")
  
  # Add this matrix to the list
  matrix_list[[folder]] <- folder_matrix
}

# Combine all matrices by column
mtx.big <- do.call(cbind, matrix_list)

# Create Seurat object
rna <- CreateSeuratObject(counts = mtx.big, min.cells = 5)

saveRDS(rna, file = "/scratch/scjp_root/scjp1/christav/milo/get_seurat_obj/FUSION_seurat.Rds")
