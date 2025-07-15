First, ran snATAC pipeline (expanded on under Figure 4).
Output of this includes an aligned BAM file, demuxlet assignment of nuclei to donors, and a broadPeak file from MACS2.

From here:

1_overall_log_reg.py - Gives how to get Figure 2F, the logistic regression comparing the mix of all 10 individuals to the FUSION cell types. This can be run with:

```
python /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/scripts/1_overall_log_reg.py --peaks /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/fusion_data/day-2.broadPeak /scratch/scjp_root/scjp1/christav/fall22_multiome/macs2/day21.broadPeak --zhang-matrix /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/fusion_data/fusion.matrix.mtx --zhang-cell-types /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/fusion_data/fusion.celltypes.tsv --zhang-peaks /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/fusion_data/fusion.peaks.tsv --chrom-sizes /scratch/scjp_root/scjp0/shared_data/reference/human/hg38/hg38.chrom_sizes --tss /scratch/scjp_root/scjp1/christav/reference/hg38.gencode.tss.bed.gz --prefix "ipscfap-fusion-"
```

2_get_barcode_indiv.R - The rest is for Fig S2C-J. This is R code to get nuclei barcodes for each individual.

3_get_indiv_bam.ipynb - Code to separate the aligned BAM file by donor based on demuxlet assignments.

4_macs2 - MACS2 commands to call peaks on BAM files specific to each individual and merge together across replicates. Also includes conversion to bigWig files for the UCSC browser screenshot in Fig SK.

5_add_peak_number.R - Adding a column for peak number to the MACS2 results to do the logistic regression.

Now can run the individual logistic regression using something like:

```
python /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/scripts/logistic-regression-comparison-to-zhang-atlas.py --peaks /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/fusion_data/day-2.broadPeak /scratch/scjp_root/scjp1/christav/fall22_multiome/macs2/12181.label.broadPeak --zhang-matrix /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/fusion_data/fusion.matrix.mtx --zhang-cell-types /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/fusion_data/fusion.celltypes.tsv --zhang-peaks /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/log_reg_celltype/fusion_data/fusion.peaks.tsv --chrom-sizes /scratch/scjp_root/scjp0/shared_data/reference/human/hg38/hg38.chrom_sizes --tss /scratch/scjp_root/scjp1/christav/reference/hg38.gencode.tss.bed.gz --prefix "12181"
```

6_comp_celltypes.R - R code to generate the heatmap used in the final figures for Fig 2F and Fig S2C-J.
