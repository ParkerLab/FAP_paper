This folder contains instructions to generate Fig 5F and Fig S5G. It uses a NextFlow pipeline to generate the fgwas results.

config.yaml - NextFlow config file, directs to broadPeak file being tested for enrichment.

In a subfolder labelled scripts/ need the files:

trait_metadata.tsv - TSV with the following columns for each GWAS: trait	description	traitname	variable_type	source

main.nf - NextFlow pipeline.

nextflow.config - NextFlow config file.

Within scripts/bin/ need the following files:

munge_sumstats.py, plot_single.R, compile_single.py

Will also need the subfolder data/sumstats_selected/ containing summary statistics from the GWAS of interest.

Can then run this pipeline with:
```
workflow='fgwas_peaks'
nfdir='/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/subcluster_peaks/fgwas_adipogenic/scripts/'
nextflow -C ${nfdir}/nextflow.config run ${nfdir}/main.nf -params-file config.yaml -with-trace trace.txt -resume
```
