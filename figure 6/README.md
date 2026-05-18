This folder contains the analysis of the iPSC-FAP subtypes, it contains code used to generate Figure 6 and Fig S6.

1_trajectory.ipynb - Pseudotime analyses used in both Fig 6B and Fig S6D-E.

2_milo.R - Neighborhood analyses on the iPSC-FAP subtypes used for Fig 6C and D.

For the binomial test, can use R's function binom.test:
```
binom.test(838, 1221, 0.6)
```
Where here it tests the number of adipogenic nuclei in the insulin sample, the total adipogenic nuclei, and the expected proportion of insulin nuclei (0.6, as insulin nuclei make up 0.6 of the full sample).

A similar method used in Figure 3 to separate by individual was used to split the ATAC bams by subcluster and call peaks separately. The reads were then normalized to 10 million reads as in Figure 4. This was used for all the UCSC images in Figure 6 and S6.

GWAS enrichment analyses - Folder of code used to generate Figure 6F.

3_SCENT.R - Performing SCENT to compare the peak of interest to genes within 50kb.
