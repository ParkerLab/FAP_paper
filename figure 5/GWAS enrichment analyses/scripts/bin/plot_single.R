#!/usr/bin/env Rscript

library(ggplot2)
library(glue)
library(cowplot)
ggplot2::theme_set(theme_cowplot())
library(tidyr)
library(dplyr)

args = commandArgs(TRUE)

savep = function(p, prefix, h=8, w=4){
    png(glue("{prefix}.png"), height=h, width=w, units="in", res = 150)
    print(p)
    dev.off()

    pdf(glue("{prefix}.pdf"), height=h, width=w)
    print(p)
    dev.off()
}
  
makeplot <- function(d){
    p <- ggplot(d, aes(y=fold, x=annotation)) +
        geom_point(aes(fill=sig), size=2, shape=21) +
        geom_errorbar(aes(ymax=ci_hi, ymin=ci_lo), width=0.5, size=0.2) +
        facet_wrap(~traitname, scales = "free_x") +
        labs(y="Fold enrichment", x="Annotation") +
        theme_bw() +
        theme(strip.text.y = element_text(size = 7),
              axis.text.x=element_text(size=9),
              panel.background=element_rect(fill="white", colour="black"),
              legend.position="bottom") +
        geom_hline(yintercept=1, size=0.5) +
        coord_flip() 
        # coord_flip(ylim = c(-2, 300))

    return(p)
}

t = read.csv("/scratch/scjp_root/scjp1/christav/fap_village_multiome/results/subcluster_peaks/fgwas_progenitors/scripts/trait_metadata.tsv", sep='\t', header=T, as.is=T)
d = read.csv(args[1], sep='\t', header=T, as.is=T)
d$sig = as.character(d$sig)

d$log2_fold = log2(exp(d$estimate))
d$log2_ci_lo = log2(exp(d$CI_lo))
d$log2_ci_hi = log2(exp(d$CI_hi))

d$fold = exp(d$estimate)
d$ci_lo = exp(d$CI_lo)
d$ci_hi = exp(d$CI_hi)
d = merge(d, t, by="trait")

print(head(d))
rename_list = c(
    "da-linked-peaks"="da-linked-peaks 0.338175Mb",
    "da-peaks"="da-peaks 1.93737Mb",
    "linked-peaks"="linked-peaks 14.0014Mb",
    "linked-tiles"="linked-tiles 82.0145Mb",
    "all-peaks"="all-peaks 106.501Mb",
    "sm-unlinked-peaks"="sm-unlinked-peaks 14.0014Mb",
    "sm-unlinked-tiles"="sm-unlinked-tiles 82.0145Mb",
    "all-unlinked-peaks"="all-unlinked-peaks 14.3161Mb",
    "all-unlinked-tiles"="all-unlinked-tiles 234.757Mb")

d$annotation <- rename_list[d$annotation]

# Convert to factor and sort
d$annotation <- factor(d$annotation, levels = rename_list)

# Sort dataframe by the factor column
d <- d[order(d$annotation), ]
# d <- d %>%
#   mutate(
#     annotation = recode(annotation, !!!rename_list), # Rename elements
#     annotation = factor(annotation, levels = rename_list) # Set factor levels to control sort order
#   ) %>%  arrange(annotation) # Sort by the specified levels

p = makeplot(d)

height = 4
savep(p, args[2], height, 5)

# #SUBSET plot
# keep = c( "type1.casnp-peaks-subsampled", "type1.other-peaks-subsampled", "type1.esnp-peaks-subsampled")
# d = d[d$annotation %in% keep,]
# d$annotation = gsub("-peaks-subsampled", "", d$annotation)
# d$annotation = factor(d$annotation, levels = c("type1.casnp", "type1.esnp", "type1.other"))
# p = makeplot(d)

# height = length(unique(d$annotation)) * 1.5
# savep(p, glue("{args[2]}.subset"), height, 5)
