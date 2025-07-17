For Figure 6, using the tool scDRS as described here: https://martinjzhang.github.io/scDRS/index.html 

1_get_h5ad.R - First need to convert the Seurat object to an h5ad object.

Then can use scDRS to compute enrichment scores using their command line interface:

```
scdrs compute-score \
    --h5ad-file /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/RNA_clustering/rna_final.h5ad\
    --h5ad-species human\
    --gs-file /scratch/scjp_root/scjp1/christav/fap_village_multiome/results/scDRS/magma_data_fig/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs\
    --gs-species human\
    --out-folder compute_score
```

2_plot_umap.ipynb - Then can plot the enrichment scores for each nuclei on the UMAP using this notebook.

3_calc_enrichment_per_cluster.R - Finally can calculate the average enrichment score per iPSC-FAP subcluster and plot as a bar plot.
