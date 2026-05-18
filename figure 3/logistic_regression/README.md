First, ran snATAC pipeline (expanded on under Figure 4).
Output of this includes an aligned BAM file and demuxlet assignment of nuclei to donors.

From here:

1_get_barcode_indiv.R - This is R code to get nuclei barcodes for each individual.

2_get_indiv_bam.ipynb - Code to separate the aligned BAM file by donor based on demuxlet assignments.

3_macs2 - MACS2 commands to call peaks on BAM files specific to each individual and merge together across replicates. Also includes conversion to bigWig files for the UCSC browser screenshot in Fig SD.

4_add_peak_number.R - Adding a column for peak number to the MACS2 results to do the logistic regression.

5_overall_log_reg.py - Gives how to perform the logistic regression comparing chromatin accessibility peaks to the FUSION cell types. This can be run with:

```
python /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/scripts/1_overall_log_reg.py --peaks /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/fusion_data/day-2-12345.label.broadPeak /scratch/scjp_root/scjp1/christav/fall22_multiome/macs2/day21-12345.label.broadPeak --zhang-matrix /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/fusion_data/fusion.matrix.mtx --zhang-cell-types /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/fusion_data/fusion.celltypes.tsv --zhang-peaks /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/fusion_data/fusion.peaks.tsv --chrom-sizes /scratch/scjp_root/scjp0/shared_data/reference/human/hg38/hg38.chrom_sizes --tss /scratch/scjp_root/scjp1/christav/reference/hg38.gencode.tss.bed.gz --prefix "12345-"
```

6_comp_celltypes.R - R code to generate the heatmap used in the final figure for Fig 2G.
