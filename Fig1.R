library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(phyloseq)
library(tidyverse)
library(phyloseq)
library(patchwork)
##################################
samples_df <- colData # Supplementary table 1
tax_mat <- vOTUtax #Supplementary table 2
otu_mat <- as.data.frame(vOTUabundance)#Supplementary table 4
# Takes first column and passes to rownames
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("contig_id") 
# Data formatting...
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("contig_id")
samples_df <- colData %>% 
  tibble::column_to_rownames("Sample")
otu_mat <- as.data.frame(otu_mat)
tax_mat <- as.matrix(tax_mat)
# Producing the  phyloseq object
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

PEARL1 <- phyloseq(OTU, TAX, samples)
PEARL <-prune_taxa(taxa_sums(PEARL1) > 0, PEARL1)
tax_table_df <- as.data.frame(as.matrix(tax_table(PEARL)))  # Convert to data frame
tax_table(PEARL) <- tax_table(as.matrix(tax_table_df))  # Convert back to matrix

# === Load sequencing stats ===
seq_stats_raw <- read_tsv("/seqfu_output.txt", col_types = cols())

seq_stats <- seq_stats_raw %>%
  filter(str_detect(File, "_R1.fastq.gz$")) %>%
  mutate(
    Filename = basename(File),
    SampleID = str_remove(Filename, "_R1.fastq.gz"),
    SampleType = if_else(str_detect(SampleID, "^E\\d{3}B(_|$)"), "Baby", "Mother"),
    PairID = str_extract(SampleID, "^E\\d{3}"),
    Month = str_extract(SampleID, "(?<=_).+"),
    SeqCount = as.numeric(`#Seq`)
  )

# === Load PEARL metadata ===
sample_metadata_df <- data.frame(sample_data(PEARL))
sample_metadata_df$SampleID <- rownames(sample_metadata_df)

# Add Pair, Month, clean up fields
sample_metadata_df <- sample_metadata_df %>%
  mutate(
    Pair = gsub("^(E\\d+).*", "\\1", SampleID),
    Months = gsub(".*_(\\w+)$", "\\1", SampleID),
    Feeding_type = trimws(Feeding_type),
    Delivery_mode = as.character(Delivery)
  )

# Feeding labels
feeding_labels <- c(
  "Exclusively_breastfed"     = "Breastfed",
  "breastfeeding_and_solids"  = "Breastfed + Solids",
  "Exclusively_formula"       = "Formula",
  "Formula_and_solids"        = "Formula + Solids",
  "Solids_only"               = "Solids only",
  "Not available"             = "Not Available",
  "Formula_and_Breastfeeding" = "Breastfed + Formula"
)

# Feeding + Delivery metadata
feeding_annotations <- sample_metadata_df %>%
  filter(Type == "Baby", !is.na(Feeding_type)) %>%
  mutate(
    Feeding_Label = dplyr::recode(Feeding_type, !!!feeding_labels, .default = "Other")
  ) %>%
  select(SampleID, Pair, Months, Feeding_Label, Delivery_mode) %>%
  distinct()

# === Join metadata ===
seq_stats <- seq_stats %>%
  mutate(PairID = as.character(PairID), Month = as.character(Month))

# Join on Baby samples
seq_stats_baby <- seq_stats %>%
  filter(SampleType == "Baby") %>%
  left_join(feeding_annotations, by = c("SampleID"))

# Mother samples
seq_stats_mother <- seq_stats %>%
  filter(SampleType == "Mother") %>%
  mutate(Feeding_Label = NA, Delivery_mode = NA)

# Combine
seq_stats <- bind_rows(seq_stats_baby, seq_stats_mother)

# Order months
time_order <- c("T1", "T2", "T3", "B", "W1", "W3", "M4", "M8", "M12", "M16", "M20", "M24")
seq_stats <- seq_stats %>%
  mutate(
    Month = factor(Month, levels = time_order),
    PairID = factor(PairID),
    SampleType = factor(SampleType, levels = c("Mother", "Baby"))
  )
# Add an identifier to the data so we know which panel we're in
seq_stats <- seq_stats %>%
  mutate(PanelPosition = if_else(SampleType == "Baby", "right", "left"))

# === p1_mother: Reads - Mother only with y-axis ===
p1_mother <- seq_stats %>%
  filter(SampleType == "Mother") %>%
  ggplot(aes(x = Month, y = PairID)) +
  geom_point(aes(size = log10(SeqCount), fill = SampleType), shape = 21, stroke = 0, alpha = 0.7) +
  scale_size_continuous(name = "log10(# Sequences)", range = c(0.5, 4)) +
  scale_fill_manual(values = c("Mother" = "cyan4")) +
  labs(title = "Sequencing Depth", x = "Timepoint", y = "Pair ID") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_blank(),
    legend.position = "none"
  )

# === p1_baby: Reads - Baby only WITHOUT y-axis ===
p1_baby <- seq_stats %>%
  filter(SampleType == "Baby") %>%
  ggplot(aes(x = Month, y = PairID)) +
  geom_point(aes(size = log10(SeqCount), fill = SampleType), shape = 21, stroke = 0, alpha = 0.7) +
  scale_size_continuous(name = "log10(# Sequences)", range = c(0.5, 4)) +
  scale_fill_manual(values = c("Baby" = "red1")) +
  labs(title = NULL, x = "Timepoint", y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )

# === p2: Feeding Mode ===
feedingmode <- seq_stats %>%
  filter(SampleType == "Baby", !is.na(Feeding_Label)) %>%
  ggplot(aes(x = Month, y = PairID)) +
  geom_point(aes(color = Feeding_Label), size = 2, alpha = 0.9) +
  scale_color_manual(
    name = "Feeding Mode",
    values = c(
      "Breastfed" = "green",
      "Breastfed + Solids" = "darkgreen",
      "Formula" = "purple",
      "Formula + Solids" = "purple4",
      "Solids only" = "orange",
      "Not Available" = "grey50",
      "Breastfed + Formula" = "orchid",
      "Other" = "black"
    )
  ) +
  facet_wrap(~ "Feeding Mode", ncol = 1) +
  labs(x = "Timepoint", y = NULL, title = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# === p3: Delivery Mode ===
deliverymode <- seq_stats %>%
  filter(SampleType == "Baby", !is.na(Delivery_mode)) %>%
  mutate(Facet = "Delivery Mode") %>%
  ggplot(aes(x = "Delivery", y = PairID)) +
  geom_point(aes(color = Delivery_mode), size = 1.5, alpha = 0.9) +
  scale_color_manual(
    name = "Delivery Mode",
    values = c(
      "Normal delivery" = "dodgerblue",
      "C- Section" = "firebrick",
      "Emergency C- Section" = "goldenrod",
      "Not available" = "grey50"
    )
  ) +
  facet_wrap(~ Facet, ncol = 1) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# === Combine all ===
p1 <- (
  (p1_mother + p1_baby + plot_layout(widths = c(1, 1))) | feedingmode | deliverymode
) + plot_layout(guides = "collect",widths = c(2.5, 1, 0.1))

p1


#####################################################

library(phyloseq)
library(ggplot2)
library(dplyr)

# Extract vOTU metadata from tax_table
votu_metadata <- as.data.frame(tax_table(PEARL))
votu_metadata$phatyp_length <- as.numeric(as.character(votu_metadata$phatyp_length))

vOTUtaxonomy$CV_contig_length == vOTUtaxonomy$phatyp_length
all(vOTUtaxonomy$CV_contig_length == vOTUtaxonomy$phatyp_length, na.rm = TRUE)
# Plot length distribution (histogram)

# Compute mean and median
length_mean <- mean(log10(votu_metadata$phatyp_length), na.rm = TRUE)
length_median <- median(log10(votu_metadata$phatyp_length), na.rm = TRUE)

# Plot
p2<- ggplot(votu_metadata, aes(x = log10(phatyp_length))) +
  geom_histogram(binwidth = 0.01, fill = "black", color = "black") +
  geom_vline(xintercept = length_mean, color = "blue", linetype = "dotted", size = 1) +
  geom_vline(xintercept = length_median, color = "red", linetype = "dotted", size = 1) +
  xlim(3, 6) +
  labs(
    x = "log10 vOTU Length in bp",
    y = "vOTU count"
  ) +
  theme_minimal()
library(dplyr)

realm_means <- votu_metadata %>%
  filter(GN_Realm %in% names(realm_colors), !is.na(phatyp_length)) %>%
  mutate(Length = as.numeric(phatyp_length)) %>%
  group_by(GN_Realm) %>%
  summarise(
    log10_mean = log10(mean(Length, na.rm = TRUE)),
    log10_median = log10(median(Length, na.rm = TRUE)),
    .groups = "drop"
  )
p2<- ggplot(votu_metadata, aes(x = log10(phatyp_length))) +
  geom_histogram(binwidth = 0.01, fill = "black", color = "black") +
  # Add one vline per realm
  geom_vline(
    data = realm_means,
    aes(xintercept = log10_median, color = GN_Realm),
    linetype = "dashed",
    size = 1
  ) +
  scale_color_manual(values = realm_colors, name = "GN Realm (Median)") +
  xlim(3, 6) +
  labs(
    #title = "vOTU Length Distribution with Realm-wise Medians",
    x = "log10 vOTU Length (bp)",
    y = "vOTU Count"
  ) +
  theme_minimal()+
  theme(
    legend.position = "none",
  )


library(dplyr)
library(ggplot2)
library(forcats)

# Custom colors
realm_colors <- c(
  "Duplodnaviria"   = "#ff1493",
  "Monodnaviria"    = "#9370db",
  "Riboviria"       = "#00ced1",
  "Varidnaviria"    = "#000080",
  "Unclassified"    = "#00ff00"
)

# Extract realm counts
votu_metadata <- as.data.frame(tax_table(PEARL))
realm_counts <- votu_metadata %>%
  count(GN_Realm, name = "Count") %>%
  filter(!is.na(GN_Realm)) %>%
  mutate(GN_Realm = fct_reorder(GN_Realm, Count))

# Plot
p3 <- ggplot(realm_counts, aes(x = GN_Realm, y = Count, fill = GN_Realm)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = realm_colors) +
  labs(
    y = "vOTU count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
library(patchwork)
combined_plot <-p1/(p2 | p3) +
  plot_layout(heights = c(2, 1))  # make p1 taller than p2/p3

ggsave("combined_figure1.svg", plot = combined_plot, width = 10, height = 10, dpi = 300)
###########################
library(phyloseq)
library(dplyr)

# Extract tax_table as a data frame
votu_metadata <- as.data.frame(tax_table(PEARL))

# Ensure proper types
votu_metadata$Length <- votu_metadata$CV_contig_length # adjust if different column name

# Filter and compute mean length per realm
votu_metadata %>%
  filter(GN_Realm %in% c("Riboviria", "Monodnaviria","Varidnaviria","Unclassified","Duplodnaviria"), !is.na(Length)) %>%
  group_by(GN_Realm) %>%
  summarise(
    mean_length = mean(Length),
    median_length = median(Length),
    count = n()
  )
# A tibble: 5 × 4
#GN_Realm      mean_length median_length count
#<chr>               <dbl>         <dbl> <int>
#1 Duplodnaviria      49575.        36953   2443
#2 Monodnaviria        5363.         5131    641
#3 Riboviria           4179.         2662.   104
#4 Unclassified       13878.        12419    337
#5 Varidnaviria       17128.        16344      6

####################################################