This folder contains the analysis of the iPSC-FAP sub-types, it contains code used to generate Figure 5 and Fig S5.

1_trajectory.ipynb - Pseudotime analyses used in both Fig 5B and Fig S5D-E.

2_milo.R - Neighborhood analyses on the iPSC-FAP sub-types used for Fig 5C and D.

For the binomial test, can use R's function binom.test:
```
binom.test(838, 1221, 0.6)
```
Where here it tests the number of adipogenic nuclei in the insulin sample, the total adipogenic nuclei, and the expected proportion of insulin nuclei (0.6, as insulin nuclei make up 0.6 of the full sample).

A similar method used in Figure 2 to separate by individual was used to split the ATAC bams by subcluster and call peaks separately. The reads were then normalized to 10 million reads as in Figure 4. This was used for all the UCSC images in Figure 5 and S5.

For Fig 5H variant calling on the vcf gave genotypes of all the individuals in the cell village:
```
bcftools query \
  -i 'ID="rs3814707"' \
  -f '[%SAMPLE\t%TGT\n]' \
  fusion_iPSC.vcf.gz
```

GWAS enrichment analyses - Folder of code used to generate Figure 5F and Fig S5G.

3_SCENT.R - Performing SCENT to compare the peak of interest to genes within 50kb.
