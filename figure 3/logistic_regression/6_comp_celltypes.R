setwd("~/links/christav/fall22_multiome/log_reg")

cell_type_comp <- read.table("12345results-cell-type-specific-peaks.tsv", header = T, sep = "\t")

cell_type_comp <- subset(cell_type_comp, our_cell_type == "12345")
cell_type_comp$coef <- cell_type_comp$coef+2
cell_type_comp$norm_coeff <- cell_type_comp$coef/max(cell_type_comp$coef)
cell_type_comp <- cell_type_comp[order(-cell_type_comp$norm_coeff),]
rownames(cell_type_comp) <- NULL

cell_type_comp$ren_cell_type <- factor(cell_type_comp$ren_cell_type, levels=c("T_cell", "Neuronal", "Adipocyte", "Type_2x", "Endothelial", "Muscle_Fiber_Mixed", "Type_2a", "Type_1", "Macrophage", "Smooth_Muscle", "Satellite_Cell", "Neuromuscular_junction", "Mesenchymal_Stem_Cell"))
cell_type_comp[is.na(cell_type_comp)] <- 0
library(ggplot2)
ggplot(cell_type_comp, aes(x=ren_cell_type, fill=norm_coeff, y=our_cell_type)) + geom_tile(color = NA) + scale_fill_gradient(low="white", high="red") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Reference Cell Type") + ylab("Day of Differentiation") + labs(fill = "Normalized Coefficient")
ggplot(cell_type_comp, aes(x=ren_cell_type, fill=norm_coeff, y=our_cell_type)) + geom_tile(color = NA) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Reference Cell Type") + ylab("Day of Differentiation") + labs(fill = "Normalized Coefficient") +
  scale_fill_gradientn(values = scales::rescale(x=c(min(cell_type_comp$norm_coeff),0, max(cell_type_comp$norm_coeff)), to = c(0,1), from = c(min(cell_type_comp$norm_coeff), max(cell_type_comp$norm_coeff))),  colors = c('white', "white", "red"))
