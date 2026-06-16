This folder contains the analysis of the iPSC-FAP subtypes, it contains code used to generate Figure 6 and Fig S6.

1_trajectory.ipynb - Pseudotime analyses used in both Fig 6B and Fig S6D-E.

2_milo.R - Neighborhood analyses on the iPSC-FAP subtypes used for Fig 6C and D.

For the binomial test, can use R's function binom.test:
```
binom.test(838, 1221, 0.6)
```
Where here it tests the number of adipogenic nuclei in the insulin sample, the total adipogenic nuclei, and the expected proportion of insulin nuclei (0.6, as insulin nuclei make up 0.6 of the full sample).

A similar method used in Figure 4 to separate by individual was used to split the ATAC bams by subcluster and call peaks separately. The reads were then normalized to 10 million reads as in Figure 5. This was used for all the UCSC images in Figure 6 and S6.

For identifying subtype-specific peaks, we used the same method as in Figure 5 (when identifying basal vs insulin peaks) but used a model testing for the impact of subtype. Since there are 3 subtypes, we ran the model twice, once with adipogenic as the reference cell type and again releveling for progenitors as the reference. A peak was designated as subtype-specific if it showed a positive coefficient and Wald test P < 0.05 in both comparisons involving that subtype, and a non-significant Wald P (> 0.05) for the comparison between the other two subtypes. For example, a peak was called adipogenic-specific if it had positive coefficients and P < 0.05 for adipogenic vs. fibrogenic and adipogenic vs. progenitor comparisons, but P > 0.05 for fibrogenic vs. progenitor. 
Here's an example calling adipogenic peaks, this would then be modified for the other two subtypes:
```
adipogenic <- subset(all_comparisons, adipo_fibro_p_GLMMNB < 0.05 & adipo_prog_p_GLMMNB < 0.05 & prog_fibro_p_GLMMNB > 0.05 & adipo_fibro_estimate > 0 & adipo_prog_estimate > 0)
```

GWAS enrichment analyses - Folder of code used to generate Figure 6F.

3_SCENT.R - Performing SCENT to compare the peak of interest to genes within 50kb.
