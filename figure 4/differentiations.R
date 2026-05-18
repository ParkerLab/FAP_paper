#To compare the efficacy of the differentiations, I categorized any cells above the 95th percentile of the negative control as expressing the FAP marker gene.
neg_control <- read.csv("~/Christa/Research/Parker Lab/Flow Cytometry Data/2024-08-01/export_neg control_Data Source - 1.csv")
pos_sample <- read.csv("~/Christa/Research/Parker Lab/Flow Cytometry Data/2024-08-01/export_55_Data Source - 1.csv")

neg_control$sample <- 'neg_control'
pos_sample$sample <- "line_55"

flowcyto <- rbind(neg_control, pos_sample)

library(ggplot2)
ggplot(flowcyto, aes(FITC.A, fill = sample)) + geom_density(alpha = 0.2)
ggplot(flowcyto, aes(PE.A, fill = sample)) + geom_density(alpha = 0.2)
ggplot(flowcyto, aes(APC.A, fill = sample)) + geom_density(alpha = 0.2)

#get cutoff from neg_control
quantile(neg_control$FITC.A, 0.95)
#541

library(tidyr)
#get percentage from positive sample
pos_sample_pos = subset(pos_sample, FITC.A > 708.25)
nrow(pos_sample_pos)/nrow(pos_sample)*100
#72.99814

#put the data in a spreadsheet for the next step


#For Figure 3B, just a simple histogram of the FACS data.

#histogram of all 30 differentiations
scatter_data <- read.csv("~/Christa/Research/Parker Lab/Flow Cytometry Data/batch30data.csv", header = TRUE)

library(ggplot2)
ggplot(scatter_data, aes(x=FAPMarker)) + geom_histogram()
