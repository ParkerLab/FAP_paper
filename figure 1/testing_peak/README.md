This has how to get the p-values in Figure 1J.

1_split_bams.py - First, split the skeletal muscle FAP bam by sub-type, then use this script to split the sub-types by individual donor.

Then need consensus peaks, can generate from broadPeak files using this:
```
cat /nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/subclustering/subclustering_FAP_02_2025_bigwigs/work/dc/4dc319598716b304a0858c3271298e/preadipogenic-narrowpeaks/preadipogenic_peaks.narrowPeak /nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/subclustering/subclustering_FAP_02_2025_bigwigs/work/34/0a063dfbdd803422870ff99d19b939/progenitors-narrowpeaks/progenitors_peaks.narrowPeak /nfs/turbo/umms-scjp/arushiv/projects/muscle-sn/analyses_hg38/subclustering/subclustering_FAP_02_2025_bigwigs/work/ba/b25b62160590267cd26d1a832e4340/fibrous-narrowpeaks/fibrous_peaks.narrowPeak > subcluster_peaks.bed

sort -k1,1 -k2,2n subcluster_peaks.bed | bedtools merge > consensus_peaks.bed

awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"} \
     {print "peak_"NR,$1,$2,$3,"."}' consensus_peaks.bed > consensus_peaks.saf
```

Next use featureCounts to get counts per individual:
```
featureCounts -T 4 -p -a consensus_peaks.saf -F SAF -o peak_counts_by_indiv.txt [list of bams per subtype per donor]
```

2_deseq2.R - Finally, use DESeq2 to test if the peak of interest is significantly different between sub-types.
