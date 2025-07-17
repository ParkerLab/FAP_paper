#comparing different states, generating the Figure 4B plot
census_seq_ins <- read.table("census_seq_output_insulin.txt", header = TRUE, sep = "\t")
census_seq_basal <- read.table("census_seq_output_basal.txt", header = TRUE, sep = "\t")
census_seq_basal$Environment <- "BASAL"
census_seq_ins$Environment <- "INSULIN"

census_seq_plot <- rbind(census_seq_basal, census_seq_ins)

ggplot(census_seq_plot, aes(x=Environment, y=REPRESENTATION, colour = Environment)) +
  geom_point() + facet_wrap(vars(DONOR)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
