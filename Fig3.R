library(dplyr)
specific_richness <- estimate_richness(ceiling(otu_table(PEARL)), measures = c("Observed", "Shannon"))
# Remove rows with SampleID matching "Blank_1.RPKM"
# Remove rows with "Blank_1.RPKM" and "Blank_2.RPKM"
specific_richness <- specific_richness %>%
  filter(!rownames(specific_richness) %in% c("Blank_1.RPKM", "Blank_2.RPKM"))
rownames(specific_richness)
# Add SampleID column
specific_richness$SampleID <- rownames(specific_richness)
samples_df$SampleID  <- rownames(samples_df)

samples_df_richness <- merge(samples_df,specific_richness, by = "SampleID")

samples_df_richness$Observed <-as.numeric(samples_df_richness$Observed)
# Define the desired order of months
time_order <- c("T1", "T2", "T3", "B", "W1", "W3", "M4", "M8", "M12", "M16", "M20", "M24")

# Reorder the levels of the "Months" factor variable
samples_df_richness$Months <- factor(samples_df_richness$Months, levels = time_order)