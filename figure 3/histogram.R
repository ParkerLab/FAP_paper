#For Figure 3B, just a simple histogram of the FACS data.

#histogram of all 30 differentiations
scatter_data <- read.csv("~/Christa/Research/Parker Lab/Flow Cytometry Data/batch30data.csv", header = TRUE)

library(ggplot2)
ggplot(scatter_data, aes(x=FAPMarker)) + geom_histogram()
