# Load libraries
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library(reshape2)
library("plyr"); packageVersion("plyr")
library(cowplot)
library(ggpattern)
library(ggthemes)
library(ggpubr)
library(ggrepel)
library(tidyverse)
library(vegan)
#####################################################################################################
specific_richness <- estimate_richness(ceiling(otu_table(PEARL)), measures = c("Observed", "Shannon"))
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
# Reshape the data for faceting multiple metrics
samples_df_richness_long <- melt(samples_df_richness, id.vars = c("SampleID", "Months", "Pair", "Type"),
                                 measure.vars = c("Observed", "Shannon"), variable.name = "Metric")
fig2A<- ggplot(samples_df_richness_long, aes(x = factor(Months), y = value)) +
  geom_boxplot(aes(fill = Type),outlier=FALSE,outlier.shape = NA) +
  geom_jitter(aes(color = Pair), width = 0.2) +
  labs(title = "Alpha_Diversity")+
  facet_grid(Metric ~ Type, scales = "free")+
  geom_line(aes(color = Pair, group = Pair), alpha = 0.7, 
            position = position_dodge(width = 0))+
  geom_smooth(aes(group = Type),colour='black', span=0.25, method="loess")+
  labs(title = "AlphaDiversity", x = "Months", y = "Alpha-Diversity Metric")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 20, angle = 45, hjust = 0.5),  # Adjust x-axis label style
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    legend.title = element_blank(),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) 

##############################
library(viridis)
library(dplyr)
library(ggplot2)
library(vegan)
library(pairwiseAdonis)
virome_changes <-  subset_samples(PEARL)
virome_changes <- prune_samples(sample_sums(virome_changes) > 0, virome_changes)
# PLOTS 1. Ensure `Type` includes both Baby and Mother
Type_levels <- unique(virome_changes@sam_data$Type)  # Get actual Type values
Type_levels <- c("Baby", "Mother")  # Set explicit levels to force inclusion

# Ensure `Type` is a factor with correct levels
virome_changes@sam_data$Type <- factor(virome_changes@sam_data$Type, levels = Type_levels)

# Perform NMDS ordination on the full dataset (NOT a subset)
ordination_nmds <- ordinate(virome_changes, method = "NMDS", distance = "bray")

# Modify Pair labels (Keep only last two characters, e.g., "02" instead of "E002")
virome_changes@sam_data$Pair <- gsub("^E0", "", virome_changes@sam_data$Pair)

# Fixed Color Palette for `Type` (Baby vs. Mother)
type_colors <- c("Baby" = "red1", "Mother" = "cyan4")  

# Define fixed shapes for each Month (ensuring they all have valid shapes)
month_shapes <- c("T1" = 15, "T2" = 16, "T3" = 17, "B" = 18, "W1" = 8, "W3" = 9, 
                  "M4" = 21, "M8" = 22, "M12" = 23, "M16" = 24, "M20" = 25, "M24" = 7)  

# Ensure Both Baby & Mother Appear in the Ellipse
fig2B <- plot_ordination(virome_changes, ordination_nmds, color = "Pair", shape = "Months") + 
  geom_point(aes(fill = Type, shape = Months, color = Pair), size = 3, alpha = 0.7) +  
  scale_shape_manual(values = month_shapes) +  # Assign valid shapes
  geom_point(size = 3) +                           
  
  # Fix Ellipse for Both Baby & Mother
  stat_ellipse(aes(group = Type, fill = Type), 
               geom = "polygon", 
               alpha = 0.25, 
               color = "black", linetype = "dashed",  # Dashed black outline for ellipses
               type = "t", 
               level = 0.95) +  
  
  # Use Color Scale for Baby & Mother
  scale_fill_manual(values = type_colors, name = "Type") +  
  
  # Restore Pair Legend
  scale_color_manual(values = c(
    "02" = "red",  
    "03" = "#377EB8",  
    "06" = "#4DAF4A",  
    "09" = "#984EA3",  
    "10" = "#FF7F00",  
    "12" = "#FFFF33",  
    "13" = "darkblue",
    "14" = "#ff1493",  
    "15" = "#F781BF",  
    "21" = "blue",  
    "26" = "red4",
    "27" = "green",
    "28" = "magenta",
    "37" = "cyan",
    "39" = "steelblue4"
  ), name = "Pair ID") +
  
  # Merge Shape and Color for Months in One Legend
  guides(
    color = guide_legend(title = "Pair ID"),
    fill = guide_legend(title = "Type"),  # Merge shape & color for Type (Baby vs. Mother)
    shape = guide_legend(title = "Months") # Shapes for Months
    # Separate legend for Pairs
  ) +
  
  # 🏷️ Improve Layout
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 0, face = "bold", hjust = 0.5)
  )

# Print the plot
print(fig2B)
############################################################################################
library(phyloseq)
library(ggplot2)
library(viridis)
library(ggnewscale)
library(tidyr)
library(dplyr)
library(magrittr)
## ── patidyr## ── palettes & shapes ────────────────────────────────────────────────────────
month_levels   <- c("M4","M8","M12","M16","M20","M24")#W3 has just 2 samples and ignored
purple_palette <- setNames(
  colorRampPalette(c("steelblue","mediumpurple","orchid4"))(length(month_levels)),
  month_levels
)
month_shapes <- c("M4"=8,"M8"=22,"M12"=23,"M16"=24,"M20"=25,"M24"=20)  # circle reused
pair_colors <- c(
  "02"="red","03"="#377EB8","06"="#4DAF4A","09"="#984EA3","10"="#FF7F00",
  "12"="#FFFF33","13"="darkblue","14"="black","15"="#F781BF","21"="blue",
  "26"="red4","27"="green","28"="magenta","37"="cyan","39"="steelblue4"
)

## ── subset Baby samples & ordination ─────────────────────────────────────────
baby <- PEARL %>%
  subset_samples(Type == "Baby") %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  subset_samples(Months %in% month_levels)

sample_data(baby)$Months <- factor(sample_data(baby)$Months, levels = month_levels)
sample_data(baby)$Pair   <- gsub("^E0","", sample_data(baby)$Pair)

## add variable stroke thickness (only M24 is thick)
sample_data(baby)$MonthStroke <- ifelse(sample_data(baby)$Months == "M24", 1.4, 0.7)
ord_nmds <- ordinate(baby, method = "NMDS", distance = "bray")

## ── plot ─────────────────────────────────────────────────────────────────────
fig2C <- plot_ordination(baby, ord_nmds, shape = "Months") +
  
  ## points
  geom_point(
    aes(fill = Pair, colour = Pair, stroke = MonthStroke),
    size = 3.2
  ) +
  scale_colour_manual(values = pair_colors, name = "Pair ID") +      # Pair legend
  scale_fill_manual(values = pair_colors, guide = "none") +       # hide fill legend
  scale_shape_manual(values = month_shapes, name = "Months") +
  
  ## second colour scale for ellipses
  new_scale_colour() +
  stat_ellipse(
    aes(group = Months, colour = Months),
    geom = "path", linetype = "solid", size = 1, level = .95
  ) +
  scale_colour_manual(values = purple_palette, name = "Months") +     # show ellipse legend
  
  ## override Months legend to show dashed line + thick stroke for M24
  guides(
    shape = guide_legend(
      title = "Months",
      override.aes = list(
        shape = unname(month_shapes),
        colour = purple_palette,
        fill = "white",
        size = 3
      ),
      order = 1
    ),
    colour = guide_legend(   # ellipse color line
      title = "Months",
      override.aes = list(
        linetype = "solid",
        size = 1.5
      ),
      order = 1
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title   = element_text(size = 12),
    plot.title   = element_text(size = 12, face = "bold", hjust = .5),
    legend.title = element_text(size = 12, face = "bold")
  )

print(fig2C)

###############################################################################
# Mothers – Beta-diversity over time #T1 has just 2 samples so not used
###############################################################################
library(phyloseq)
library(ggplot2)
library(viridis)
library(ggnewscale)

## 1 ─ palettes & shapes -------------------------------------------------------
month_levels   <- c("T2","T3","B","W1","W3","M4")
purple_palette <- setNames(
  colorRampPalette(c("steelblue","mediumpurple","orchid4"))(length(month_levels)),
  month_levels
)

## re-use circle for M4 but thicken its outline
month_shapes  <- c("T2"=21,"T3"=22,"B"=23,"W1"=24,"W3"=25,"M4"=8)
month_stroke  <- c("T2"=.6,"T3"=.6,"B"=.6,"W1"=.6,"W3"=.6,"M4"=1.2)

pair_colors <- c(
  "02"="red","03"="#377EB8","06"="#4DAF4A","09"="#984EA3","10"="#FF7F00",
  "12"="#FFFF33","13"="darkblue","14"="black","15"="#F781BF","21"="blue",
  "26"="red4","27"="green","28"="magenta","37"="cyan","39"="steelblue4"
)

## 2 ─ subset Mother samples & ordination -------------------------------------
library(magrittr)   # or just load dplyr/phyloseq which re-export %>%
moms <- PEARL %>%                               # start with the full object
  subset_samples(Type == "Mother") %>%          # keep mothers only
  prune_samples(sample_sums(.) > 0, .) %>%      # drop empty libraries
  subset_samples(Months %in% month_levels & Months != "T1")  # filter months
sample_data(moms)$Months <- factor(sample_data(moms)$Months, levels = month_levels)
sample_data(moms)$Pair   <- gsub("^E0","", sample_data(moms)$Pair)

## attach the variable outline width
sample_data(moms)$MonthStroke <- month_stroke[ sample_data(moms)$Months ]
ord_nmds <- ordinate(moms, method = "NMDS", distance = "bray")
## 3 ─ plot -------------------------------------------------------------------
fig2D <- plot_ordination(moms, ord_nmds, shape = "Months") +
  # points: fill/shape = Months, outline colour = Pair, stroke varies
  geom_point(
    aes(fill = Pair, colour = Pair, stroke = MonthStroke),
    size = 3.2, alpha = .83
  ) +
  scale_colour_manual(values = pair_colors, name = "Pair ID") +   # Pair legend
  scale_fill_manual(values = pair_colors, guide = "none") +    # hide fill legend
  scale_shape_manual(values = month_shapes,  name = "Months") +   # shape legend
  
  # new colour scale for ellipses so they don't inherit Pair colours
  new_scale_colour() +
  stat_ellipse(
    aes(group = Months, colour = Months),
    geom     = "path",
    linetype = "solid",
    size     = 1,
    level    = .95
  ) +
  scale_colour_manual(values = purple_palette, guide = "none") +  # no extra legend
  
  # customise the Months legend to show dashed line + correct stroke widths
  ## override Months legend to show dashed line + thick stroke for M24
  guides(
    shape = guide_legend(
      title = "Months",
      override.aes = list(
        shape = unname(month_shapes),
        colour = purple_palette,
        fill = "white",
        size = 3
      ),
      order = 1
    ),
    colour = guide_legend(   # ellipse color line
      title = "Months",
      override.aes = list(
        linetype = "solid",
        size = 1.5
      ),
      order = 1
    )
    ## Pair ID legend (order 2) comes from scale_colour_manual above
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title   = element_text(size = 18),
    plot.title   = element_text(size = 18, face = "bold", hjust = .5),
    legend.title = element_text(size = 16, face = "bold")
  )

print(fig2D)

#############################################################################################
PEARL_baby <- subset_samples(PEARL, Type == "Baby"  & !Months =="W3") # Just 2 samples in W3
PEARL_baby <- prune_samples(sample_sums(PEARL_baby) > 0, PEARL_baby)
PEARL_baby <- prune_taxa(taxa_sums(PEARL_baby) > 0, PEARL_baby)
dist_matrix_baby <- phyloseq::distance(PEARL_baby, method = "bray")
# Extract metadata
metadata_baby <- as(sample_data(PEARL_baby), "data.frame")
# Ensure row names match the distance matrix
valid_samples <- rownames(as.matrix(dist_matrix_baby))
metadata_baby <- metadata_baby[valid_samples, , drop = FALSE]
# Check metadata structure
print(dim(metadata_baby))
#Betadispersion Baby
dispersion_test_baby <- betadisper(dist_matrix_baby, metadata_baby$Months) 
#Order months chronologically
disp_df <- data.frame(
  Distance = dispersion_test_baby$distances,
  Month = dispersion_test_baby$group
  )
month_levels <- c( "M4", "M8", "M12", "M16", "M20", "M24")
disp_df$Month <- factor(disp_df$Month, levels = month_levels)

fig2E <- ggplot(disp_df, aes(x = Month, y = Distance)) +
geom_boxplot(fill = "red1", color = "black", alpha = 0.7) +
geom_jitter(width = 0.2, alpha = 0.5) +
theme_minimal(base_size = 14) +
  labs(
    title = "Beta Diversity Dispersion Across Months (Baby Samples)",
    x = "Month",
    y = "Distance to Centroid (Bray-Curtis)"
    )
#############################################################################################
PEARL_mother <- subset_samples(PEARL, Type == "Mother" & !Months=="T1") # T1 has just 2 samples 
# Remove empty samples (no reads)
PEARL_mother <- prune_samples(sample_sums(PEARL_mother) > 0, PEARL_mother)
# Remove taxa that are completely absent
PEARL_mother <- prune_taxa(taxa_sums(PEARL_mother) > 0, PEARL_mother)
# Compute Bray-Curtis distance matrix
dist_matrix_mother <- phyloseq::distance(PEARL_mother, method = "bray")
# Extract metadata
metadata_mother <- as(sample_data(PEARL_mother), "data.frame")
# Ensure row names match the distance matrix
valid_samples_mother <- rownames(as.matrix(dist_matrix_mother))
metadata_mother <- metadata_mother[valid_samples_mother, , drop = FALSE]
# Check metadata structure
print(dim(metadata_mother))
#Order months chronologically
month_levels <- c( "T2", "T3", "B", "W1", "W3", "M4")
dispersion_test_mother <- betadisper(dist_matrix_mother, metadata_mother$Months)
# Visualize dispersion differences using a boxplot
dispersion_df <- data.frame(
  DistanceToCentroid = dispersion_test_mother$distances,
  Month = dispersion_test_mother$group
)

dispersion_df$Month <- factor(dispersion_test_mother$group, levels = month_levels)
fig2F<- ggplot(dispersion_df, aes(x = Month, y = DistanceToCentroid)) +
  geom_boxplot(fill = "cyan4", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, color = "black") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold")
) +
  labs(
   title = "Beta Diversity Dispersion (PERMDISP) - Mothers",
   x = "Month",
   y = "Distance to Centroid (Bray-Curtis)"
  )
#########################################################################################
#ASSOCIATED STATISTICS
###############Alpha diversity GLMM and LMM as the data is not in complete blocks,longitudinal and skewed. Repeated measures test wont work (Friedmann)################################
samples_df_richness$Months <- factor(samples_df_richness$Months, levels = time_order)
samples_df_richness$Pair   <- factor(samples_df_richness$Pair)
samples_df_richness$Type   <- factor(samples_df_richness$Type)  # Baby / Mother
library(glmmTMB)
baby_data <- subset(samples_df_richness, Type == "Baby")
model_obs_baby <- glmmTMB(
  Observed ~ Months + (1 | Pair),
  family = nbinom2,
  data = baby_data
)
summary(model_obs_baby)
library(car)
anova(model_obs_baby, type="II")
drop1(model_obs_baby, test="Chisq")

library(lme4)
library(lmerTest)

model_shannon_baby <- lmer(
  Shannon ~ Months + (1 | Pair),
  data = baby_data
)

anova(model_shannon_baby)
summary(model_shannon_baby)
###########################################################

mother_data <- subset(samples_df_richness, Type == "Mother")
model_obs_mother <- glmmTMB(
  Observed ~ Months + (1 | Pair),
  family = nbinom2,
  data = mother_data
)

summary(model_obs_mother)
anova(model_obs_mother)

model_shannon_mother <- lmer(
  Shannon ~ Months + (1 | Pair),
  data = mother_data
)

anova(model_shannon_mother)
summary(model_shannon_mother)

library(emmeans)

# For infant richness
emm <- emmeans(model_obs_baby, pairwise ~ Months, type = "response")
emm$contrasts  # All pairwise comparisons
###############################################################################################################################################
#Betadiversity dispersion test
dispersion_test_baby <- betadisper(dist_matrix_baby, metadata_baby$Months)
permutest(dispersion_test_baby)
tukey_disp <- TukeyHSD(dispersion_test_baby)
print(tukey_disp)
dispersion_test_mother <- betadisper(dist_matrix_mother, metadata_mother$Months)
permutest(dispersion_test_mother)
tukey_disp <- TukeyHSD(dispersion_test_mother)
print(tukey_disp)
########### permanova#####
# Run PERMANOVA with Months as the grouping factor for baby
adonis_result_baby <- adonis2(dist_matrix_baby ~ Months, data = metadata_baby, permutations = 999)
print(adonis_result_baby)
library(pairwiseAdonis)
pairwise_result <- pairwise.adonis2(dist_matrix_baby ~ Months, data = metadata_baby, permutations = 999)
print(pairwise_result)
#For p.value adjustment FDR corrections Extract all pairwise results (excluding $parent_call and $p.adjusted)
results_list <- pairwise_result[!names(pairwise_result) %in% c("parent_call", "p.adjusted")]
# Extract p-values and comparisons
pair_names <- names(results_list)
pvals <- sapply(results_list, function(x) x["Model", "Pr(>F)"])
#Apply FDR correction
pvals_fdr <- p.adjust(pvals, method = "fdr")
#Combine into a data frame
pairwise_df <- data.frame(
  Comparison = pair_names,
  p.value = pvals,
  p.adjusted = pvals_fdr,
  row.names = NULL
  )
# View result
print(pairwise_df)
# Run PERMANOVA with Months as the grouping factor for Mothers
adonis_result_mother <- adonis2(dist_matrix_mother ~ Months, data = metadata_mother, permutations = 999)
print(adonis_result_mother)
pairwise_result_mother <- pairwise.adonis2(dist_matrix_mother ~ Months, data = metadata_mother, permutations = 999)
print(pairwise_result_mother)#No significant p-values seen here, so p.val adj not required



#######################################################################################################
ibrary(viridis)
library(phyloseq)
otu_long <- psmelt(PEARL)
# Verify structure
print(head(otu_long))

bubble_data <- otu_long %>%
  group_by(vOTU, Months, GN_Realm, Pair, Type) %>%
  summarise_at(vars(Abundance), list(~sum(., na.rm = TRUE))) %>%
  ungroup()
########################################################SUPPLEMENTARY FIG -4 EARLY HOSTS###############
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)

# -----------------------------
# SETTINGS
# -----------------------------
early_months <- c("W3","M4")
month_order  <- c("T1","T2","T3","B","W1","W3","M4","M8","M12","M16","M20","M24")

# -----------------------------
# HOST MAP
# -----------------------------
recent_UHGV_ammended_taxonomy<- read.csv(file.choose(),header=TRUE)
recent_UHGV_ammended_taxonomy$genome_id <- fix_contig(recent_UHGV_ammended_taxonomy$genome_id)

recent_UHGV_clean_host <- left_join(
  recent_UHGV_ammended_taxonomy,
  votu_host_clean,
  by = c("genome_id" = "contig_fix")
)

host_map <- as.data.frame(recent_UHGV_ammended_taxonomy) %>%
  mutate(
    vOTU = genome_id,
    Host_Label = coalesce(iPHOP_host, ICTV_Broad_Host_Domain)#For figure purpose used genus level info, where not available MAG derived genus level info used, and if thats not available, ICTV broad host is used.
  ) %>%
  select(vOTU, Host_Label)

# -----------------------------
# INFANT W3/M4 RECORDS WITH HOSTS
# -----------------------------
inf_early <- bubble_data %>%
  filter(Type == "Baby", Abundance > 0, Months %in% early_months) %>%
  left_join(host_map, by = c("OTU" = "vOTU")) %>%
  filter(!is.na(Host_Label),
         !str_detect(tolower(Host_Label), "unclassified"),
         !str_detect(tolower(Host_Label), "not available")) %>%
  mutate(
    Months = factor(trimws(Months), levels = month_order),
    Host_Label = str_to_sentence(Host_Label)  # sentence case
  ) %>%
  select(Pair, Months, OTU, Host_Label) %>%
  distinct()

# -----------------------------
# UNIQUE vOTUs PER HOST × MONTH
# -----------------------------
host_counts <- inf_early %>%
  group_by(Host_Label, Months) %>%
  summarise_at(vars(OTU), ~ n_distinct(.)) %>%
  rename_at(vars(OTU), ~ "unique_vOTUs") %>%
  ungroup()

# Order hosts by total across W3+M4
host_order <- host_counts %>%
  group_by(Host_Label) %>%
  summarise_at(vars(unique_vOTUs), ~ sum(.)) %>%
  arrange(desc(unique_vOTUs)) %>%
  pull(Host_Label)

plot_df <- host_counts %>%
  mutate(
    Months = factor(as.character(Months), levels = early_months),
    Host_Label = factor(Host_Label, levels = host_order)
  )

# -----------------------------
# HEATMAP
# -----------------------------
ggplot(plot_df, aes(x = Months, y = Host_Label, fill = unique_vOTUs)) +
  geom_tile(color = "grey85") +
  scale_fill_viridis_c(name = "Unique vOTUs") +
  labs(
    title = "Infant virome hosts at early timepoints (W3 & M4)",
    x = "Month", y = "Predicted host"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic")  # italic labels
  )

###################To put in paper, how many early vOTUs and hoe many hosts ?
library(dplyr)
library(stringr)

# -----------------------------
# SETTINGS
# -----------------------------
early_months <- c("W3","M4")
month_order  <- c("T1","T2","T3","B","W1","W3","M4","M8","M12","M16","M20","M24")

# -----------------------------
# HOST MAP (ensure 1 row per vOTU)
# -----------------------------
# recent_UHGV_clean_host must already exist; if not, build it as you did earlier
host_map <- as.data.frame(recent_UHGV_ammended_taxonomy) %>%
  mutate(
    vOTU       = genome_id,
    Host_Label = coalesce(iPHOP_host, ICTV_Broad_Host_Domain)#Iphop family level prediction is most reliable, where iphop prediction not available Broad host domain is used
  ) %>%
  select(vOTU, Host_Label) %>%
  distinct(vOTU, .keep_all = TRUE)

# -----------------------------
# 1) ALL INFANT vOTUs AT W3/M4 (NO HOST FILTER)
# -----------------------------
early_votus_all <- bubble_data %>%
  filter(Type == "Baby", Abundance > 0, Months %in% early_months) %>%
  select(Pair, Months, OTU) %>%
  distinct() %>%
  mutate(Months = factor(Months, levels = month_order))

# overall unique vOTUs across W3+M4
votu_totals_overall <- early_votus_all %>%
  summarise_at(vars(OTU), ~ n_distinct(.)) %>%
  rename_at(vars(OTU), ~ "total_unique_vOTUs")

# per-month unique vOTUs
votu_per_month <- early_votus_all %>%
  group_by(Months) %>%
  summarise_at(vars(OTU), ~ n_distinct(.)) %>%
  rename_at(vars(OTU), ~ "unique_vOTUs") %>%
  ungroup()

print(votu_totals_overall)
print(votu_per_month)

# -----------------------------
# 2) HOST COUNTS (ONLY AMONG KNOWN/VALID HOSTS)
# -----------------------------
early_with_hosts <- early_votus_all %>%
  left_join(host_map, by = c("OTU" = "vOTU")) %>%
  mutate(
    Host_Label_clean = Host_Label %>%
      str_replace_all("_", " ") %>%
      str_to_sentence()
  )
#########################

# overall unique hosts across W3+M4 (exclude unknowns/unclassified)
host_totals_overall <- early_with_hosts %>%
  filter(!is.na(Host_Label_clean),
         !str_detect(tolower(Host_Label_clean), "unclassified"),
         !str_detect(tolower(Host_Label_clean), "not available"),
         Host_Label_clean != "") %>%
  summarise_at(vars(Host_Label_clean), ~ n_distinct(.)) %>%
  rename_at(vars(Host_Label_clean), ~ "total_unique_hosts")

# per-month unique hosts (same exclusions)
host_per_month <- early_with_hosts %>%
  filter(!is.na(Host_Label_clean),
         !str_detect(tolower(Host_Label_clean), "unclassified"),
         !str_detect(tolower(Host_Label_clean), "not available"),
         Host_Label_clean != "") %>%
  group_by(Months) %>%
  summarise_at(vars(Host_Label_clean), ~ n_distinct(.)) %>%
  rename_at(vars(Host_Label_clean), ~ "unique_hosts") %>%
  ungroup()

print(host_totals_overall)
print(host_per_month)
###########################Supplementary Fig 5- Bracken Abundance ###############
# ── libraries ───────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(grDevices)     # for palette.colors()
})

# ── parameters ──────────────────────────────────────────────
top_N   <- 50                    # show 50 most abundant genera
infile  <- file.choose()               # interactively pick all_bracken_long_unified.tsv
out_png <- "genus_stacked_unified.png" # output plot file

# ── 1. read table & add Genus column ────────────────────────
long <- read_tsv(infile, show_col_types = FALSE) %>%
  mutate(
    new_est_reads = as.numeric(gsub(",", "", New_est_reads)),   # ensure numeric
    Genus         = str_extract(TaxName, "^[A-Za-z0-9_]+")      # first token
  ) %>%
  filter(!is.na(Genus))         # <-- remove rows where genus extraction failed

# 1. Create unified_frac: per-sample share of all classified reads
long <- long %>%
  group_by(Sample) %>%
  mutate(unified_frac = New_est_reads / sum(New_est_reads,na.rm=TRUE)) %>%
  ungroup()

# ── 2. aggregate unified fractions to genus per sample ─────
genus <- long %>% 
  group_by(Sample, Genus) %>% 
  summarise_at(
    vars(unified_frac),          # columns to summarise
    ~ sum(., na.rm = TRUE)       # summarising function
  ) %>% 
  ungroup()                      # equivalent to `.groups = "drop"`

# ── 3. keep TOP N genera; others → "Other" ──────────────────
top_genera <- genus %>% 
  group_by(Genus) %>% 
  summarise_at(                     # <─ superseded but still available
    vars(unified_frac),             # column(s) to summarise
    ~ mean(., na.rm = TRUE)         # summary function
  ) %>% 
  mutate(meanFrac = unified_frac) %>%   # give the new column a friendly name
  arrange(desc(meanFrac)) %>% 
  slice_head(n = top_N) %>%            # keep the first `top_N` genera
  pull(Genus)                           # return just the Genus vector

summarise_atgenus <- genus %>% 
  mutate(
    Genus = ifelse(Genus %in% top_genera, Genus, "Other")  # collapse low-abundance genera
  ) %>% 
  group_by(Sample, Genus) %>% 
  summarise_at(                       # ← superseded but still supported
    vars(unified_frac),               # column(s) to reduce
    ~ sum(., na.rm = TRUE)            # summary function
  ) %>% 
  mutate(Frac = unified_frac) %>%     # give the result a clearer name
  ungroup()                           # equivalent to .groups = "drop"


# ── AFTER genus aggregation, BEFORE choosing palette ─────────
# label each sample as "Baby" (contains "BM") or "Mother"
genus <- summarise_atgenus %>%
  mutate(
    SampleType = ifelse(
      grepl("B(?:M|W|Birth)", Sample, ignore.case = TRUE),  # BM, BW, or BBirth
      "Baby",
      "Mother"
    )
  )

# ensure samples stay grouped by type in the x-axis order you have
genus <- genus %>%
  arrange(SampleType, Sample)
genus$Genus <- factor(genus$Genus)

# ── 4. distinct colour palette (base R “Polychrome 36”) ─────
core_genera <- setdiff(levels(genus$Genus), "Other")
pal_core    <- if (length(core_genera) <= 36) {
  palette.colors(length(core_genera), palette = "Polychrome 36")
} else {
  colorRampPalette(palette.colors(36, "Polychrome 36"))(length(core_genera))
}
names(pal_core) <- core_genera
pal <- c(pal_core, Other = "grey70")   
# ensure legend ordered by overall abundance
level_order <- genus %>% 
  group_by(Genus) %>% 
  summarise_at(vars(Frac), sum) %>% 
  arrange(desc(Frac)) %>% 
  pull(Genus)

genus$Genus <- factor(genus$Genus, levels = level_order)
# ── 5. plot ─────────────────────────────────────────────────
# ── plot with faceting by SampleType ─────────────────────────
p <- ggplot(genus, aes(Sample, Frac, fill = Genus)) +
  geom_col(width = 1, colour = "black", size = 0.25) +
  scale_fill_manual(values = pal, name = "Genus") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~SampleType, ncol = 1, scales = "free_x") +       # <<< facet here
  labs(y = "Relative abundance (%)",
       x = NULL,
       title = paste("Top", top_N,
                     "genera in Baby vs Mother samples (Bracken, 100 % stacked)")) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.x = element_text(face = "bold"),
        legend.position = "right")

ggsave("genus_stacked_faceted.png", p, width = 12, height = 8, dpi = 300)
cat("Plot saved → genus_stacked_faceted.png\n")
###################SUPPLEMENTARY FIG -6 ###############################################
library(dplyr)
library(ggplot2)
library(stringr)

# Month ordering
month_order <- c( "W3", "M4", "M8", "M12", "M16", "M20", "M24")
#W3_baby_samples <- pers_df_all_with_type %>% filter(HostType=="Baby", Months =="W3")
# Step 1: Prepare persistent vOTUs from pers_df_all given in Fig4.R
baby_df <- bubble_data %>%
  filter(Abundance > 0, Type == "Baby") %>%
  distinct(Pair, OTU, Months)

# Step 2: Add taxonomy
recent_UHGV_ammended_taxonomy<-read.csv(file.choose(),header=TRUE)#Supplementary table 2



tax_data <- as.data.frame(recent_UHGV_ammended_taxonomy)
#tax_data$contig_fix <- tax_data$vOTU
baby_df$contig_fix <- fix_contig(baby_df$OTU)

baby_df <- baby_df %>%
  left_join(tax_data,
            by = c("contig_fix"="genome_id")) %>%
  mutate(
    Host_Genus_final  = coalesce(Final_Host_Genus, ICTV_Broad_Host_Domain),
    Host_Phylum_final = coalesce(Host_Phylum, ICTV_Broad_Host_Domain)
  ) %>%
  filter(!is.na(Family_Cluster_id), !is.na(Months)) %>%
  mutate(Months = factor(trimws(Months), levels = month_order))

library(dplyr)
library(ggplot2)
library(stringr)

# --- Step 3: First detection month per Family_Cluster_id ---
cluster_summary <- baby_df %>%
  group_by(Family_Cluster_id) %>%
  summarise_at(
    vars(Months),
    ~ levels(.)[min(as.numeric(.))]
  ) %>%
  ungroup()
names(cluster_summary)[2] <- "first_month"

# --- Step 4: Host lookup ---
host_lookup <- baby_df %>%
  group_by(Family_Cluster_id) %>%
  summarise_at(
    vars(Host_Genus_final),
    ~ paste0(unique(na.omit(.)), collapse = ", ")
  ) %>%
  ungroup()
names(host_lookup)[2] <- "host_families"

# --- Step 5: Join + collapse bare "Cluster###" into one group ---
cluster_summary <- cluster_summary %>%
  left_join(host_lookup, by = "Family_Cluster_id") %>%
  mutate(
    first_month = factor(first_month, levels = month_order),
    Display_ID = ifelse(
      str_detect(Family_Cluster_id, "^Cluster\\d+$"),
      "Other Clusters",
      as.character(Family_Cluster_id)
    ),
    Count = 1
  ) %>%
  # Collapse rows so "Other Clusters" counts as one per month
  group_by(first_month, Display_ID) %>%
  summarise_at(vars(Count), sum) %>%
  ungroup()

# Total counts per month (for labels above bars)
month_totals <- cluster_summary %>%
  group_by(first_month) %>%
  summarise_at(vars(Count), sum) %>%
  ungroup() %>%
  rename(Total = Count)

# --- Step 6: Plot ---
ggplot(cluster_summary, aes(x = first_month, y = Count)) +
  geom_col(color = "black", fill = "gray90") +
  geom_text(
    aes(label = Display_ID),
    position = position_stack(vjust = 0.5),
    size = 2.7,
    check_overlap = TRUE
  ) +
  geom_text(
    data = month_totals,
    aes(x = first_month, y = Total, label = Total),
    vjust = -0.35,
    size = 4,
    inherit.aes = FALSE
  ) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Month-wise Colonization of Persistent vOTU Clusters (Grouped)",
    x = "First Detection Month",
    y = "Number of vOTU Clusters"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
##################better aesthetic##########
library(dplyr)
library(ggplot2)
library(stringr)
# optional (makes labels avoid overlap more nicely)
suppressWarnings({ library(ggrepel) })

# 1) First detection month per Family_Cluster_id
first_appear <- baby_df %>%
  mutate(Months = factor(Months, levels = month_order)) %>%
  filter(!is.na(Family_Cluster_id), !is.na(Months)) %>%
  group_by(Family_Cluster_id) %>%
  summarise(first_month = levels(Months)[min(as.numeric(Months))], .groups = "drop") %>%
  mutate(first_month = factor(first_month, levels = month_order))

# 2) Collapse plain "Cluster###" (no parentheses) into "Other Clusters"; keep everything else
first_appear <- first_appear %>%
  mutate(Display_ID = ifelse(str_detect(Family_Cluster_id, "^Cluster\\d+$"),
                             "Other Clusters", as.character(Family_Cluster_id)))

# 3) Split labels: we will *list* only the named/annotated ones; count "Other Clusters"
labels_df <- first_appear %>%
  filter(Display_ID != "Other Clusters")

other_counts <- first_appear %>%
  filter(Display_ID == "Other Clusters") %>%
  count(first_month, name = "other_n")

# 4) For each month, give each label a simple vertical rank so they don’t overlap
labels_ranked <- labels_df %>%
  arrange(first_month, Display_ID) %>%
  group_by(first_month) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# 5) Per-month totals to print above the month
month_totals <- first_appear %>%
  count(first_month, name = "Total")

# 6) Plot: a clean first-appearance “timeline”
p <- ggplot() +
  # light baseline row so months read like a strip
  geom_hline(yintercept = 0, color = "grey90") +
  # dots for each first appearance (equal spacing by rank; no proportional heights)
  geom_point(
    data = labels_ranked,
    aes(x = first_month, y = rank),
    size = 1.8, alpha = 0.9
  ) +
  # labels for named families/annotated clusters
  { if (requireNamespace("ggrepel", quietly = TRUE))
    geom_text_repel(
      data = labels_ranked,
      aes(x = first_month, y = rank, label = Display_ID),
      size = 3, min.segment.length = 0.1, max.overlaps = 50, box.padding = 0.2, seed = 1
    ) else
      geom_text(
        data = labels_ranked,
        aes(x = first_month, y = rank, label = Display_ID),
        size = 3, vjust = -0.2, check_overlap = TRUE
      )
  } +
  # add a small “+n other” note if there were plain Cluster### that month
  geom_text(
    data = other_counts,
    aes(x = first_month, y = 0.5, label = paste0("+", other_n, " other")),
    size = 3.2, color = "grey20"
  ) +
  # total counts above each month
  geom_text(
    data = month_totals,
    aes(x = first_month, y = max(labels_ranked$rank %||% 0) + 1.5, label = Total),
    size = 4.2
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "First Month of Detection for vOTU Families / Annotated Clusters",
    x = "First Detection Month",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title.position = "plot"
  )

print(p)
###################SUPPLEMENTARY8######################################################
library(dplyr)
library(ggplot2)
library(patchwork)

# Add delivery mode
samples_df_richness_long <- samples_df_richness_long %>%
  mutate(Delivery = ifelse(Pair %in% c("E002","E009","E010","E037"), "C-section", "Vaginal"))

# --- Mothers ---
mother_months <- c("T1","T2","T3","B","W1","W3","M4")
inf_div_mom <- samples_df_richness_long %>%
  filter(Type == "Mother", Metric == "Observed") %>%
  mutate(Months = factor(Months, levels = mother_months))

A <- ggplot(inf_div_mom, aes(x = Months, y = value, fill = Delivery, color = Delivery)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, position = position_dodge(width = 0.8), alpha = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0,
                                             dodge.width = 0.8), size = 1.8, alpha = 0.9) +
  scale_fill_manual(values = c("C-section" = "#E41A1C", "Vaginal" = "#377EB8")) +
  scale_color_manual(values = c("C-section" = "#E41A1C", "Vaginal" = "#377EB8")) +
  labs(title = "Mothers: virome alpha diversity by delivery mode",
       x = "Month", y = "Observed vOTU richness", fill = "Delivery", color = "Delivery") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# --- Babies ---
baby_months <- c("W3","M4","M8","M12","M16","M20","M24")
inf_div_baby <- samples_df_richness_long %>%
  filter(Type == "Baby", Metric == "Observed") %>%
  mutate(Months = factor(Months, levels = baby_months))

B <- ggplot(inf_div_baby, aes(x = Months, y = value, fill = Delivery, color = Delivery)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, position = position_dodge(width = 0.8), alpha = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0,
                                             dodge.width = 0.8), size = 1.8, alpha = 0.9) +
  scale_fill_manual(values = c("C-section" = "#E41A1C", "Vaginal" = "#377EB8")) +
  scale_color_manual(values = c("C-section" = "#E41A1C", "Vaginal" = "#377EB8")) +
  labs(title = "Infants: virome alpha diversity by delivery mode",
       x = "Month", y = "Observed vOTU richness", fill = "Delivery", color = "Delivery") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# Combine
deliverymode_combined <- A | B + plot_layout(nrow = 1)
deliverymode_combined
############################################################################################
##################Supplementary-7 Baby–mother within the same pair share more similar viromes than random baby–mother across different pairs##############

library(phyloseq)
library(vegan)
library(dplyr)
library(tibble)
library(tidyr)

# Get binary OTU matrix
otu_bin <- otu_table(PEARL)
#otu_bin[otu_bin > 0] <- 1

# Remove empty samples (columns with all 0s)
non_empty_samples <- colSums(otu_bin) > 0
otu_bin_filtered <- otu_bin[, non_empty_samples]

# # Compute binary Sørensen dissimilarity (Bray-Curtis on presence/absence)
diss <- vegdist(t(otu_bin_filtered), method = "bray")  # this is binary Sørensen

# Convert to data frame and add SampleID
meta <- data.frame(sample_data(PEARL))
meta$SampleID <- rownames(meta)

# Filter metadata to include only samples used in dissimilarity matrix
meta_filtered <- meta %>% filter(SampleID %in% colnames(otu_bin_filtered))
# Convert Converts the square dissimilarity matrix to long-form dataframe
diss_df <- as.matrix(diss) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "Dissimilarity") %>% #Columns: Sample1, Sample2, Dissimilarity.
  filter(Sample1 < Sample2)#Sample1 < Sample2 ensures each pair is listed once (upper triangle only).

# Add metadata for both samples
diss_df <- diss_df %>%
  left_join(meta %>% select(Sample1 = SampleID, Type1 = Type, Pair1 = Pair), by = "Sample1") %>% #Joins metadata (Type and Pair) for both samples in the dissimilarity pair.
  left_join(meta %>% select(Sample2 = SampleID, Type2 = Type, Pair2 = Pair), by = "Sample2") %>%
  filter(Type1 != Type2) %>%  # Only Baby-Mother comparisons #Ensures that only Baby vs. Mother comparisons are retained, Discards Baby-Baby or Mother-Mother sample comparisons.
  mutate(Relationship = ifelse(Pair1 == Pair2, "Within_Pair", "Between_Pair")) #Ensures both samples are from the same pair (same baby–mother duo,within_pair) else, between_pair

# Statistical test
wilcox.test(Dissimilarity ~ Relationship, data = diss_df)#W = 790322, p-value = 0.002131Sorensons and W = 813397, p-value = 4.867e-05 bray-curtis

###### Sample balancing####Still significany p-val
set.seed(42)

# Subsample 332 rows randomly from Between_Pair to match Within_Pair
between_subsample <- diss_df %>% 
  filter(Relationship == "Between_Pair") %>% 
  sample_n(200)

within_subset <- diss_df %>% 
  filter(Relationship == "Within_Pair")

balanced_df <- bind_rows(within_subset, between_subsample)

# Re-run Wilcoxon on balanced data
wilcox.test(Dissimilarity ~ Relationship, data = balanced_df)

library(ggplot2)
library(dplyr)

# Make sure Similarity is present
diss_df <- diss_df %>%
  mutate(Similarity = 1 - Dissimilarity)

effect_sizes <- diss_df %>%
  group_by(Relationship) %>%
  summarise(
    median_sim = median(Similarity),
    IQR_sim = IQR(Similarity),
    mean_sim = mean(Similarity)
  )

effect_sizes

# A tibble: 2 × 4
#Relationship median_sim IQR_sim mean_sim
#<chr>             <dbl>   <dbl>    <dbl>
#  1 Between_Pair   0.000622 0.00445  0.00868
#2 Within_Pair    0.00113  0.0146   0.0291 

###Plot 

ggplot(diss_df, aes(x = Relationship, y = Similarity, fill = Relationship)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  scale_y_log10() +
  scale_fill_manual(values = c("Within_Pair" = "#1b9e77", "Between_Pair" = "#d95f02")) +
  theme_classic(base_size = 14) +
  ylab("Bray–Curtis similarity (log scale)") +
  xlab("")