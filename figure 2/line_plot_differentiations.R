#This code generates supp fig 2B/C, and figure 2C.

#single line, over course of differentiation, numbers were taken from flow cytometry data
gene <- c("Tra-1-60", "NT5E", "Tra-1-60", "NT5E", "Tra-1-60", "NT5E")
exp_level <- c(73.4, 7.25, 1.8, 11.5, 0.42, 74.8)
data <- data.frame(gene, exp_level)
data$Day <- c("Day -2", "Day -2", "Day 6", "Day 6", "Day 21", "Day 21")
data$Day <- factor(data$Day, levels=c("Day -2", "Day 6", "Day 21"))
library(ggplot2)
ggplot(data, aes(x=Day, y=exp_level, group = gene, color=gene)) + geom_line(size = 1.5) + geom_point(size = 4) + xlab("Day of Differentiation") + ylab("Percentage of Cells Expressing") + theme_bw() + scale_color_manual(values = c("Tra-1-60" = "orange", "NT5E" = "purple"))

#different single line, over course of differentiation, numbers were taken from flow cytometry data
gene <- c("Tra-1-60", "NT5E", "Tra-1-60", "NT5E", "Tra-1-60", "NT5E")
exp_level <- c(94.9, 6.96, 0.86, 0.45, 0.37, 99.8)
data <- data.frame(gene, exp_level)
data$Day <- c("Day -2", "Day -2", "Day 6", "Day 6", "Day 21", "Day 21")
data$Day <- factor(data$Day, levels=c("Day -2", "Day 6", "Day 21"))
library(ggplot2)
ggplot(data, aes(x=Day, y=exp_level, group = gene, color=gene)) + geom_line(size = 1.5) + geom_point(size = 4) + xlab("Day of Differentiation") + ylab("Percentage of Cells Expressing") + theme_bw() + scale_color_manual(values = c("Tra-1-60" = "orange", "NT5E" = "purple"))


#aggregate of all differentiations I have time course data for (Fig 2C)
diff <- read.csv("~/Christa/Research/Parker Lab/Flow Cytometry Data/all_differentiations.csv")
diff$day <- factor(diff$day, levels=c("Day -2", "Day 6", "Day 21"))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

diff2 <- data_summary(diff, varname="exp", 
                    groupnames=c("marker", "day"))
# Convert dose to a factor variable
diff2$day=as.factor(diff2$day)

ggplot(diff2, aes(x=day, y=exp, group = marker, color=marker)) + geom_line(size = 1.5) + geom_point(size = 4) + xlab("Day of Differentiation") + ylab("Percentage of Cells Expressing") + theme_bw() + scale_color_manual(values = c("Tra-1-60" = "orange", "NT5E" = "purple")) +
  geom_errorbar(aes(ymin=exp-sd, ymax=exp+sd, width=.2))
