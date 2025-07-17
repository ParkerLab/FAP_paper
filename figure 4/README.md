This folder contains the bulk of the snRNA+ATAC analyses, it contains code used to generate Figure 4 and Fig S4.

Before multiome analyses was started, census-seq was ran on light whole genome sequencing of the cell village as outlined in: https://github.com/chventresca/census-seq

1_census_seq.R - Contains code to generate Figure 4B plot of census-seq data.

For multiome processing, snRNA and snATAC processing was performed as outlined in these pipelines: 
https://github.com/ParkerLab/snRNAseq-NextFlow 
https://github.com/porchard/snATACseq-NextFlow

2_qc.ipynb - A notebook for selecting pass QC nuclei based on outlined parameters.

Prior to doublet detection the vcf had to be sorted using this code:

```
samtools view -H /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/ATAC_results/prune/lipids-hg38.pruned.bam | grep SQ | perl -pe 's/.*SN:(.*)\tLN.*/$1/' > header.txt

python /scratch/scjp_root/scjp1/christav/fap_village_multiome/scripts/sort-vcf-specify-order.py --sort-order-file /scratch/scjp_root/scjp1/christav/fap_village_multiome/input/vcf/header.txt --vcf-in /scratch/scjp_root/scjp1/christav/fap_village_multiome/input/vcf/fusion_iPSC_hg38_norm.variant.sorted2.rsid.qual.vcf --vcf-out /scratch/scjp_root/scjp1/christav/fap_village_multiome/input/vcf/fusion_iPSC_hg38_norm.variant.sorted3.rsid.qual.vcf

gzip /scratch/scjp_root/scjp1/christav/fap_village_multiome/input/vcf/fusion_iPSC_hg38_norm.variant.sorted3.rsid.qual.vcf
```

Doublet detection was then performed as outlined here: https://github.com/ParkerLab/Multiome-Doublet-Detection-NextFlow 
Nuclei that were not labelled as doublets according to demuxlet were selected as final pass QC nuclei.

The provided pipelines generate bigwig files for ATAC, to generate for RNA can do something like this:

```
bedtools genomecov -ibam ../prune/basal-hg38.before-dedup.bam -bg > basal.bedgraph

LC_COLLATE=C sort -k1,1 -k2,2n basal.bedgraph > basal.sorted.bedgraph

bedGraphToBigWig basal.sorted.bedgraph /scratch/scjp_root/scjp0/shared_data/reference/human/hg38/hg38.chrom_sizes basal_rna.bw
```

The resulting UCSC browser session is available at: https://genome.ucsc.edu/s/chventresca/FAP%20Village%20Multiome

snRNA based analyses:

3_clustering.R - R code to generate Figure 4C-G and Figure S4A.

4_correlation.R - R code to generate Figure SB.

5_trajectory.ipynb - Code to generate Figure 4D and Figure S4B and C.

6_differential_expression.R - Code to run MAST and identify differentially expressed genes between insulin and basal samples. Contains instructions for Figure 4H.

7_gene_set_enrich_analysis.ipynb - Code to run clusterProfiler and identify enriched gene sets in the differentially expressed genes.

8_gse_plots.R - Code to generate Figure 4J based on gene set enrichment. (Can use both 7_gene_set_enrich_analysis.ipynb and/or 8_gse_plots.R to get Table S1.)

9_milo.R - Code for Figure S4E and F, milo on the iPSC-FAPs.

snATAC based analyses:

First need to filter and sort fragment files from initial ATAC analyses to just final pass QC nuclei. 

10_filter.py - Script used here for filtering.

```
/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/differential_peaks/scripts/filter.py /scratch/scjp_root/scjp1/christav/fap_village_multiome/work/ATAC_work/bc/3a5e4e271d1d1412976b96aebd5e0f/basal-hg38.frags.bed /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/RNA_clustering/basal_atac_barcode_list.txt 4 > /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/ATAC_results/fragment-file/basal.frags.filtered.final.bed

sort -k4,4 -k1,1 -k2n,2 -k3n,3 /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/ATAC_results/fragment-file/basal.frags.filtered.final.bed > /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/ATAC_results/fragment-file/basal.frags.filtered.final.sorted.bed
```

11_snapatac_joint.ipynb - Then can use snapATAC2 to get joint peaks between insulin and basal samples.

12_peak_R_obj.R - R script to extract peak information from snapATAC2 object, also needed.

13_glmnb.R - R script to perform differential peak identification using a Negative Binomial Generalized Linear Mixed Model (NBGLMM), this was performed in 50 batches (there's a line at the end to select which batch to run). Commented out section at the beginning gives how to generate the metadata file used in the linear mixed model.

After that you'll have 50 files of FAP peaks, can use the line below to concat them all into the same file.

```
awk 'NR==1 {header=$_} FNR==1 && NR!=1 { $_ ~ $header getline; } {print}' *.tsv > FAP_all.tsv
```

14_filter_peaks.R - Gives code to apply FDR correction to the peaks, identify significantly different peaks, and find the values for Figure 4L.

15_volcano_plot_peaks.R - Code to generate Figure 4M.
