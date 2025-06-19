# FAP_paper
This repo contains the code used to develop the manuscript "Discovering stimulatory state specific type 2 diabetes GWAS mechanisms with single-cell multi-omics on iPSC-derived fibro-adipogenic progenitors."

Rough outline of code to be included:

FUSION analyses (figure 1):
1. FAP abundance analyses
2. Milo analyses

Differentiation analyses (figure 2):
1. FACS differentiation analyses
2. Endpoint FACS histogram
3. Logistic regression to FUSION

Insulin stimulated glucose uptake assay analyses (figure 3):
1. Luminesence analyses
2. LMM analyses

FAP Village Multiome analyses (figure 4-6):
1. Census-seq analysis
2. Initial RNA/ATAC analyses
3. Multiome QC/doublet detection (initial, redo with sctricter parameters)
4. Setting up UCSC browser
5. RNA clustering, RNA/ATAC clustering, clustering with FUSION
6. Trajectory analyses
7. Differential gene expression (MAST, GSEA/gProfiler)
8. Allelic bias (WASP)
9. Differential peak identification, GWAS enrichment of differential peaks
10. Gene expression GWAS enrichment
