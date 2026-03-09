library(dplyr)
specific_richness <- estimate_richness(ceiling(otu_table(PEARL)), measures = c("Observed", "Shannon"))