#This gives Figure 3D and E
#Pre-processing of the glucose stimulated insulin uptake assay involved averaging background fluorescence off of the plate (I had a full row of empty wells) and subtracting that from each well.

setwd("C:/Users/chris/OneDrive/Documents/Christa/Research/Parker Lab/Plate Reader/08.2024")
library(readxl)
figure <- read_excel("figure.xlsx")
library(ggplot2)


#to adjust based on cell counts, set up a linear model (for a million cells, what is expected?):
cell_counts <- read.table("cell counts.txt", sep = "\t", header = T)

mod <- lm(Luminesence ~ Total.Cells, data = basal_ins_counts)
expectation_at_1M_cells <- predict(mod, data.frame(Total.Cells = c(1e6))) 
residuals <- resid(mod)
per_sample_estimation_for_1M_cells_basalins <- expectation_at_1M_cells + residuals
basal_ins_counts$corrected.lum <- per_sample_estimation_for_1M_cells_basalins
ggplot(basal_ins_counts, aes(x=Total.Cells, y=corrected.lum, colour = Line)) + geom_point() +
  geom_text(label=basal_ins_counts$Line) + ggtitle("Cell Counts vs. Corrected Basal + Insulin")
ggplot(basal_ins_counts, aes(x=Total.Cells, y=corrected.lum)) + geom_point() +
  geom_smooth(method=lm, se=FALSE) + ggtitle("Cell Counts vs. Corrected Basal + Insulin") +
  geom_text(x = 2000000, y = 10000000, label = lm_eqn(basal_ins_counts), parse = TRUE)

basal_counts <- merge(cell_counts, basal)
mod <- lm(Luminesence ~ Total.Cells, data = basal_counts)
expectation_at_1M_cells <- predict(mod, data.frame(Total.Cells = c(1e6))) 
residuals <- resid(mod)
per_sample_estimation_for_1M_cells_basal <- expectation_at_1M_cells + residuals
basal_counts$corrected.lum <- per_sample_estimation_for_1M_cells_basal
ggplot(basal_counts, aes(x=Total.Cells, y=corrected.lum, colour = Line)) + geom_point() +
  geom_text(label=basal_counts$Line) + ggtitle("Cell Counts vs. Corrected Basal")
ggplot(basal_counts, aes(x=Total.Cells, y=corrected.lum)) + geom_point() +
  geom_smooth(method=lm, se=FALSE) + ggtitle("Cell Counts vs. Corrected Basal") +
  geom_text(x = 2000000, y = 10000000, label = lm_eqn(basal_counts), parse = TRUE)

mod <- lm(Luminesence ~ Total.Cells, data = insulin_ins_counts)
expectation_at_1M_cells <- predict(mod, data.frame(Total.Cells = c(1e6))) 
residuals <- resid(mod)
per_sample_estimation_for_1M_cells_insulinins <- expectation_at_1M_cells + residuals
insulin_ins_counts$corrected.lum <- per_sample_estimation_for_1M_cells_insulinins
ggplot(insulin_ins_counts, aes(x=Total.Cells, y=corrected.lum, colour = Line)) + geom_point() +
  geom_text(label=insulin_ins_counts$Line) + ggtitle("Cell Counts vs. Corrected Insulin + Insulin")
ggplot(insulin_ins_counts, aes(x=Total.Cells, y=corrected.lum)) + geom_point() +
  geom_smooth(method=lm, se=FALSE) + ggtitle("Cell Counts vs. Corrected Insulin + Insulin") +
  geom_text(x = 1000000, y = 5000000, label = lm_eqn(insulin_ins_counts), parse = TRUE)

insulin_counts <- merge(cell_counts, insulin)
mod <- lm(Luminesence ~ Total.Cells, data = insulin_counts)
expectation_at_1M_cells <- predict(mod, data.frame(Total.Cells = c(1e6))) 
residuals <- resid(mod)
per_sample_estimation_for_1M_cells_insulin <- expectation_at_1M_cells + residuals
insulin_counts$corrected.lum <- per_sample_estimation_for_1M_cells_insulin
ggplot(insulin_counts, aes(x=Total.Cells, y=corrected.lum, colour = Line)) + geom_point() +
  geom_text(label=insulin_counts$Line) + ggtitle("Cell Counts vs. Corrected Insulin")
ggplot(insulin_counts, aes(x=Total.Cells, y=corrected.lum)) + geom_point() +
  geom_smooth(method=lm, se=FALSE) + ggtitle("Cell Counts vs. Corrected Insulin") +
  geom_text(x = 1000000, y = 4000000, label = lm_eqn(insulin_counts), parse = TRUE)
#saved as corr_counts

#then inverse rank normalization
basal_prenorm <- subset(corr_counts, Environment == "Basal")
basalins_prenorm <- subset(corr_counts, Environment == "Basal + Insulin")
basal_ins_prenorm <- rbind(basal_prenorm, basalins_prenorm)
basal_ins_prenorm <- na.omit(basal_ins_prenorm)
normalized_basal <- qnorm((rank(basal_ins_prenorm$corrected.lum,na.last="keep")-0.5)/sum(!is.na(basal_ins_prenorm$corrected.lum)))
basal_ins_prenorm$quantnorm_lum <- normalized_basal
insulin_prenorm <- subset(corr_counts, Environment == "Insulin")
insulinins_prenorm <- subset(corr_counts, Environment == "Insulin + Insulin")
insulin_ins_prenorm <- rbind(insulin_prenorm, insulinins_prenorm)
insulin_ins_prenorm <- na.omit(insulin_ins_prenorm)
normalized_insulin <- qnorm((rank(insulin_ins_prenorm$corrected.lum,na.last="keep")-0.5)/sum(!is.na(insulin_ins_prenorm$corrected.lum)))
insulin_ins_prenorm$quantnorm_lum <- normalized_insulin
#saved as new_corr_counts

#got the mean of counts and set to quantnorm_lum_mean, plot below is Figure 3D
ggplot(figure_mean, aes(x=Environment, y=quantnorm_lum, colour = Line)) + 
  geom_jitter(width = 0.1, alpha=0.5) + geom_point(aes(Environment, quantnorm_lum_mean), shape=95, size=20, alpha=0.5) +
  ylab("Normalized Lum") 

#a t-test was used to calculate significance:
t.test(basal$quantnorm_lum, basal_ins$quantnorm_lum, paired = TRUE, alternative = "two.sided")
t.test(insulin$quantnorm_lum, insulin_ins$quantnorm_lum, paired = TRUE, alternative = "two.sided")

#Finally set up a linear mixed model with donor metadata and plotting it (Figure 3E)
library(tidyr)
library(ggpubr)
#adding metadata to the luminescence data
metadata <- read_excel("../../iPSC lines info/fusion_ipsc_batch_assignment.sex-age-bmi.50ksim (2).xlsx")
#subset
metadata_sub <- metadata %>% separate(NYSCF_ID, c("Line", "A", "B", "C"))
metadata_final <- metadata_sub[,c("Line", "sex", "bmi", "ogtt", "age")]
#adding metadata
corr_counts_meta <- merge(corr_counts, metadata_final, by = "Line")
corr_counts_meta$age <- as.numeric(corr_counts_meta$age)
corr_counts_meta$bmi <- as.numeric(corr_counts_meta$bmi)
#linear mixed model
library(lme4)
library(lmerTest)

#first basal
mixed.lmer <- lmer(quantnorm_lum ~ age + sex + ogtt + bmi + (1|Line), data = basal_ins)
summary(mixed.lmer)
#save as basal_ins_model to plot it
basal_ins_model$Trait <- factor(basal_ins_model$Trait, levels=c("ogttT2D", "sexM", "bmi", "age"))
basal_ins_model$color_group <- with(basal_ins_model, ifelse(
  pval < 0.05 & Estimate > 0, "Up (pval < 0.05)",
  ifelse(pval < 0.05 & Estimate < 0, "Down (pval < 0.05)", "NS")
))
ggplot(basal_ins_model, aes(x=Estimate, y=Trait, color = color_group)) + geom_point() +
  geom_errorbarh(aes(xmin = (Estimate - Std_Error),xmax = (Estimate + Std_Error))) +
  scale_color_manual(values = c("Up (pval < 0.05)" = "blue",
                                "Down (pval < 0.05)" = "red",
                                "NS" = "gray")) +
  labs(
    title = "Basal Insulin Model",
    x = "Estimate",
    y = "Trait",
    color = "Significance"
  ) +
  theme_minimal()

#now insulin
mixed.lmer <- lmer(quantnorm_lum ~ age + sex + ogtt + bmi + (1|Line), data = insulin_ins)
summary(mixed.lmer)
#save as insulin_ins_model to plot it
insulin_ins_model$Trait <- factor(insulin_ins_model$Trait, levels=c("ogttT2D", "sexM", "bmi", "age"))
insulin_ins_model$color_group <- with(insulin_ins_model, ifelse(
  pval < 0.05 & Estimate > 0, "Up (pval < 0.05)",
  ifelse(pval < 0.05 & Estimate < 0, "Down (pval < 0.05)", "NS")
))
ggplot(insulin_ins_model, aes(x=Estimate, y=Trait, color = color_group)) + geom_point() +
  geom_errorbarh(aes(xmin = (Estimate - Std_Error),xmax = (Estimate + Std_Error))) +
  scale_color_manual(values = c("Up (pval < 0.05)" = "blue",
                                "Down (pval < 0.05)" = "red",
                                "NS" = "gray")) +
  labs(
    title = "Insulin Insulin Model",
    x = "Estimate",
    y = "Trait",
    color = "Significance"
  ) +
  theme_minimal()
