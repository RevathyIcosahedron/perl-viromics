####FIGURE4########################################
#bubble data from figure3.R
#Defining persistence as below.
process_pair <- function(pair_id, bubble_data) {
  pair_data <- bubble_data %>% filter(Pair == pair_id)
  
  # Define vOTUs
  baby_present <- pair_data %>% filter(Type == "Baby", Abundance > 0) %>% pull(OTU) %>% unique()
  mother_present <- pair_data %>% filter(Type == "Mother", Abundance > 0) %>% pull(OTU) %>% unique()
  
  shared_votus <- intersect(baby_present, mother_present)
  baby_only <- setdiff(baby_present, mother_present)
  mother_only <- setdiff(mother_present, baby_present)
  
  # Define persistence function
  calc_persistence <- function(df, type, pair_id) { # each baby and each mother in each pair this is calculated
    df_filtered <- df %>%
      ungroup() %>%
      filter(Pair == pair_id, Type == type, Abundance > 0) %>%
      mutate(
        TimeGroup = case_when(
          type == "Mother" & Months %in% c("B", "W1", "W3") ~ "Early", #Pool Early Mother's samples as time interval between sampling is short
          TRUE ~ Months
        )
      )
    
    n_samples <- df_filtered %>%
      summarise(n = n_distinct(TimeGroup)) %>%
      pull(n)
    
    df_filtered %>%
      group_by(OTU) %>%
      summarise(n_timepoints = n_distinct(TimeGroup), .groups = "drop") %>%
      mutate(Persistent = n_timepoints >= ceiling(n_samples / 2)) # vOTU present atleast in 50% of time points collected
  }
  
  baby_pers <- calc_persistence(bubble_data, "Baby", pair_id)
  mother_pers <- calc_persistence(bubble_data, "Mother", pair_id)
  
  # Skip this pair entirely if baby or mother persistence is NULL
  if (is.null(baby_pers) | is.null(mother_pers)) return(NULL)
  
  pers_df <- full_join(baby_pers, mother_pers, by = "OTU", suffix = c("_Baby", "_Mother")) %>%
    mutate(
      Class = case_when(
        OTU %in% shared_votus ~ "Shared",
        OTU %in% baby_only ~ "Unique to Baby",
        OTU %in% mother_only ~ "Unique to Mother",
        TRUE ~ NA_character_
      ),
      Subclass = case_when(
        Class == "Unique to Baby" & Persistent_Baby ~ "unique-Persistent in Baby",
        Class == "Unique to Baby" & !Persistent_Baby ~ "unique-transient in Baby",
        Class == "Unique to Mother" & Persistent_Mother ~ "unique-Persistent in Mother",
        Class == "Unique to Mother" & !Persistent_Mother ~ "unique-transient in Mother",
        Class == "Shared" & Persistent_Baby & Persistent_Mother ~ "Persistent in Both",
        Class == "Shared" & Persistent_Baby & !Persistent_Mother ~ "Persistent in Baby, detected in Mother",
        Class == "Shared" & !Persistent_Baby & Persistent_Mother ~ "Persistent in Mother, detected in Baby",
        Class == "Shared" & !Persistent_Baby & !Persistent_Mother ~ "Shared but not Persistent in Either"
      ),
      Interpretation = case_when(
        Subclass == "Persistent in Both" ~ "Strongly Shared",
        Subclass %in% c("Persistent in Baby, detected in Mother", "Persistent in Mother, detected in Baby") ~ "Shared in either baby/mother and not persistent in the other",
        Subclass == "Shared but not Persistent in Either" ~ "Transiently Shared between baby and mother",
        Subclass == "unique-Persistent in Baby" ~ "Persistent Unique to Baby",
        Subclass == "unique-transient in Baby" ~ "Transient Unique to Baby",
        Subclass == "unique-Persistent in Mother" ~ "Persistent Unique to Mother",
        Subclass == "unique-transient in Mother" ~ "Transient Unique to Mother",
        TRUE ~ NA_character_
      ),
      Pair = pair_id
    )
  
  return(pers_df)
}
all_pairs <- unique(bubble_data$Pair)
# Remove pairs where only one baby sample was available, so persistence definition is not possible.
pers_df_all <- purrr::map_dfr(setdiff(unique(bubble_data$Pair), c("E013", "E028")), 
                              ~process_pair(.x, bubble_data))
###############Figure 4-B ############################
library(ggplot2)
library(dplyr)
library(scales)
library(ggpattern)
library(patchwork)

# 1. Interpretation setup
interpretation_levels <- c(
  "Strongly Shared",
  "Shared in either baby/mother and not persistent in the other",
  "Transiently Shared between baby and mother",
  "Persistent Unique to Mother",
  "Transient Unique to Mother",
  "Persistent Unique to Baby",
  "Transient Unique to Baby"
)

interpretation_colors <- c(
  "Strongly Shared" = "#0b1e77",
  "Shared in either baby/mother and not persistent in the other" = "#377EB8",
  "Transiently Shared between baby and mother" = "#7570b3",
  "Persistent Unique to Baby" = "#a6761d",
  "Transient Unique to Baby" = "#e6ab02",
  "Persistent Unique to Mother" = "#a6761d",
  "Transient Unique to Mother" = "#e6ab02"
)

# 2. Preprocess input data
df_base <- pers_df_all %>%
  filter(!is.na(Interpretation)) %>%
  distinct(Pair, OTU, Interpretation) %>%
  mutate(
    Interpretation = factor(Interpretation, levels = interpretation_levels),
    pattern = ifelse(Interpretation %in% c("Persistent Unique to Baby", "Transient Unique to Baby"), "circle", "none")
  )

# 3. Count plot
plot_count <- df_base %>%
  count(Pair, Interpretation, pattern, name = "count") %>%
  mutate(pattern_alpha = ifelse(pattern == "circle", 1, 0)) %>%
  ggplot(aes(x = Pair, y = count, fill = Interpretation, pattern = pattern)) +
  geom_bar_pattern(
    stat = "identity",
    position = "stack",
    aes(pattern_alpha = pattern_alpha),
    pattern_fill = "black",
    pattern_density = 0.1,
    pattern_spacing = 0.05,
    pattern_key_scale_factor = 0.5
  ) +
  scale_fill_manual(values = interpretation_colors) +
  scale_pattern_manual(values = c("none" = "none", "circle" = "circle")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none") +   # <-- remove legend
  labs(
    title = "vOTU Count by Interpretation",
    y = "vOTU Count",
    x = "Pair",
    caption = "Circle pattern indicates Baby-specific interpretations"
  )
# 4. Proportion plot
plot_prop <- df_base %>%
  count(Pair, Interpretation, pattern, name = "count") %>%
  group_by(Pair) %>%
  mutate(
    prop = count / sum(count),
    pattern_alpha = ifelse(pattern == "circle", 1, 0)
  ) %>%
  ggplot(aes(x = Pair, y = prop, fill = Interpretation, pattern = pattern)) +
  geom_bar_pattern(
    stat = "identity",
    position = "stack",
    aes(pattern_alpha = pattern_alpha),
    pattern_fill = "black",
    pattern_density = 0.1,
    pattern_spacing = 0.05,
    pattern_key_scale_factor = 0.5
  ) +
  scale_fill_manual(values = interpretation_colors) +
  scale_pattern_manual(values = c("none" = "none", "circle" = "circle")) +
  scale_y_continuous(labels = percent_format()) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none") +   # <-- remove legend
  labs(
    title = "vOTU Proportion by Interpretation",
    y = "Proportion of vOTUs",
    x = "Pair",
    caption = "Circle pattern indicates Baby-specific interpretations"
  )
# 5. Combine and export
final_plot <- plot_count / plot_prop +
  plot_layout(ncol = 1)

ggsave(
  filename = "shared_unique_patterns_combined_circle.pdf",
  plot = final_plot,
  width = 6,    # half A4 width
  height =12,
  units = "in"
)

###########################Figure 4-C###################################################
# Extract unique-only vOTUs with persistence status and host
library(dplyr)
library(tidyr)
library(ggplot2)


# Step 1: Prepare data using summarise_at with Baby exclusions
test_df <- pers_df_all %>%
  filter(Interpretation %in% c(
    "Persistent Unique to Baby", "Transient Unique to Baby",
    "Persistent Unique to Mother", "Transient Unique to Mother"
  )) %>%
  mutate(
    Host = ifelse(grepl("Baby", Interpretation), "Baby", "Mother"),
    Persistence = ifelse(grepl("Persistent", Interpretation), "Persistent", "Transient")
  ) %>%
  filter(!(Pair %in% c("E026", "E015") & Host == "Baby")) %>%  # Exclude only Baby from E026 and E015
  distinct(Pair, OTU, Host, Persistence) %>%
  group_by(Pair, Host, Persistence) %>%
  summarise_at(vars(OTU), funs(vOTU_count = n())) %>%
  ungroup()

library(ggplot2)

ggplot(test_df, aes(x = Persistence, y = vOTU_count, fill = Host)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_jitter(width = 0.15, size = 2) +
  facet_wrap(~ Host) +
  labs(title = "Unique vOTUs: Transient vs Persistent", y = "vOTU count") +
  theme_minimal()
#########################Stats
baby_df <- test_df %>% filter(Host == "Baby") %>%
  pivot_wider(names_from = Persistence, values_from = vOTU_count)
mother_df <- test_df %>% filter(Host == "Mother") %>%
  pivot_wider(names_from = Persistence, values_from = vOTU_count)

wilcox.test(baby_df$Transient, baby_df$Persistent, paired = FALSE)#V = 66, p-value = 0.0009766;Mann-Whitney- W = 117, p-value = 0.0002348
wilcox.test(baby_df$Persistent , mother_df$Persistent, paired =FALSE)#W = 18, p-value = 0.002131
wilcox.test(baby_df$Transient , mother_df$Transient, paired =FALSE)#W = 91, p-value = 0.2767
wilcox.test(mother_df$Transient , mother_df$Persistent, paired =FALSE)#W = 116.5, p-value = 0.1062
## ── 2.  MEDIAN + IQR FOR BABIES ────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(ggplot2)

#----- 1. summary stats -------------------------------------------------
baby_stats <- test_df %>%                       # start from the long table
  filter(Host == "Baby") %>%                    # only infant samples
  filter(!(Pair %in% c("E026", "E015"))) %>%    # drop two Baby pairs
  group_by(Persistence) %>%                     # Persistent / Transient
  summarise_at(vars(vOTU_count),                # <<<   summarise_at!
               list(median = ~median(., na.rm = TRUE),
                    q1     = ~quantile(., 0.25, na.rm = TRUE),
                    q3     = ~quantile(., 0.75, na.rm = TRUE),
                    iqr    = ~IQR(., na.rm = TRUE))) %>% 
  ungroup() %>% 
  mutate(label = sprintf("%.1f [%.1f–%.1f] (IQR)", median, q1, q3))

print(baby_stats)
#####################SUPPLEMENTARY Figure -9 & 10 #################################################
library(dplyr)
library(ggplot2)
library(viridis)

# Filter strongly shared from pers_df_all
strongly_shared_df <- pers_df_all %>%
  filter(Subclass%in% c("Persistent in Both"))
##filter(Subclass %in% c("Persistent in Baby, detected in Mother", "Shared but not Persistent in Either","Persistent in Mother, detected in Baby"))#Other shared vOTUs
unique(pers_df_all$Subclass)
recent_UHGV_ammended_taxonomy$genome_id <- fix_contig(recent_UHGV_ammended_taxonomy$genome_id)
tax_data <- as.data.frame(recent_UHGV_ammended_taxonomy) %>% mutate(
  Final_Host_Genus  = coalesce(Final_Host_Genus, ICTV_Broad_Host_Domain),
  Host_Phylum_final = coalesce(Host_Phylum,ICTV_Broad_Host_Domain)
)

# ============================
bubble_data_extended <- bubble_data %>%
  filter(Abundance > 0) %>%
  mutate(contig_fix = fix_contig(OTU)) %>%                          
  left_join(tax_data[, c("contig_fix", "Family_Cluster_id", "Final_Host_Genus")],
            by = c("contig_fix" = "contig_fix")) %>%                          
  select(Pair, OTU, Months, Type, Abundance, Family_Cluster_id, Final_Host_Genus,contig_fix) %>%
  distinct()

# Combine strongly shared vOTUs with abundance and taxonomy

persistent_vOTUs_host <- strongly_shared_df %>%
  left_join(
    bubble_data_extended %>%
      select(Pair, contig_fix, Months, Type, Abundance, Family_Cluster_id, Final_Host_Genus) %>%
      distinct(),
    by = c("Pair", "OTU"="contig_fix")                         
  ) %>%
  filter(!is.na(Abundance), !is.na(Type))

persistent_vOTUs_host <- persistent_vOTUs_host %>%
  mutate(Final_Host_Genus = ifelse(is.na(Final_Host_Genus), "Not available", Final_Host_Genus))

persistent_vOTUs_host$Type <- factor(
  persistent_vOTUs_host$Type,
  levels = c("Mother", "Baby", "Pairs"),
  ordered = TRUE
)

# Collapse bacterial families per cluster
collapsed_families <- persistent_vOTUs_host %>%
  mutate(Family.y = ifelse(is.na(Final_Host_Genus), "Not available", Final_Host_Genus)) %>%
  distinct(Family_Cluster_id, Final_Host_Genus) %>%
  group_by(Family_Cluster_id) %>%
  summarise_at(vars(Final_Host_Genus), ~ paste(sort(unique(.)), collapse = ", ")) %>%
  ungroup()
names(collapsed_families)[2] <- "Family_Collapsed"

# Join collapsed family and prepare plotting vars
persistent_vOTUs_host <- persistent_vOTUs_host %>%
  left_join(collapsed_families, by = "Family_Cluster_id") %>%
  mutate(
    ClusterID = Family_Cluster_id,
    Pair = gsub("^E0", "", Pair),
    Months = factor(Months, levels = c("T1", "T2", "T3", "B", "W1", "W3", 
                                       "M4", "M8", "M12", "M16", "M20", "M24")),
    Type = factor(Type, levels = c("Mother", "Baby", "Pairs"), ordered = TRUE),
    Final_Host_Genus = gsub(", ", ",\n", Family_Collapsed)
  ) %>%
  filter(!is.na(Abundance), !is.na(Type))

# Dot data for Baby-only panel (x = Pair)
pair_dot_data <- persistent_vOTUs_host %>%
  filter(Type == "Baby") %>%
  distinct(ClusterID, Pair, Final_Host_Genus) %>%
  mutate(Type = "Pairs", Pair = factor(Pair))
pair_dot_data$Type <- factor(
  pair_dot_data$Type,
  levels = c("Mother", "Baby", "Pairs"),
  ordered = TRUE
)

# Colors for pairs
pair_colors <- c(
  "02" = "red", "03" = "#377EB8", "06" = "#4DAF4A", "09" = "#984EA3", 
  "10" = "#FF7F00", "12" = "#FFFF33", "13" = "darkblue", "14" = "black", 
  "15" = "#F781BF", "21" = "blue", "26" = "red4", "27" = "green", 
  "28" = "magenta", "37" = "cyan", "39" = "steelblue4"
)

# Final plot (unchanged)
final_bubble_plot <- ggplot() +
  geom_point(data = persistent_vOTUs_host, aes(
    x = Months, y = ClusterID, size = log10(Abundance), fill = log10(Abundance)
  ), alpha = 0.8, shape = 21, color = "black") +
  geom_point(data = pair_dot_data, aes(x = Pair, y = ClusterID, color = Pair), size = 4, shape = 16) +
  scale_size_continuous(range = c(1,5)) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "log(Abundance)") +
  scale_color_manual(values = pair_colors, name = "Pair ID") +
  facet_grid(Final_Host_Genus ~ Type, scales = "free", space = "free_y", switch = "y") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10, face = "bold"),
    strip.text.y = element_text(size = 10, face = "bold", angle = 0),
    strip.text.x = element_text(size = 14, face = "bold"),
    panel.spacing = unit(1, "lines"),
    legend.position = "right",
    strip.placement = "outside",
    strip.text.y.left = element_text(size = 5, angle = 0, face = "bold"),
    strip.background.y.left = element_rect(fill = "gray95", color = NA),
    panel.spacing.y = unit(0.2, "lines")
  ) +
  labs(x = "Months/Pair ID", y = "Family Cluster ID", size = "log(Abundance)")


# Auto height
a4_width <- 11.7
bubble_height <- max(6, min(length(unique(persistent_vOTUs_host$ClusterID)) / 3, 20))
ggsave("All_shared_not_strongly_shared_vOTUs_bubble.pdf", final_bubble_plot,
       width = a4_width, height = bubble_height, dpi = 600, units = "in", limitsize = FALSE)

message("Saved A4 plot with strongly shared vOTUs. Height: ", bubble_height, " inches") 

# Reduced height version
bubble_height <- max(6, min(length(unique(persistent_vOTUs_host$ClusterID)) / 6, 14))
ggsave("Stronglyshared_vOTUs_BubblePlot_A4_newtax_data.pdf", final_bubble_plot,
       width = a4_width, height = bubble_height, dpi = 600, units = "in", limitsize = FALSE)
print(paste("Saved PDF with A4 width and reduced height:", bubble_height, "inches"))

#############################What is the actual numbers of sharing -supplementary table 5#############################################################################

library(dplyr)
library(tidyr)

# 1. Define categories
shared_categories <- c(
  "Strongly Shared",
  "Shared in either baby/mother and not persistent in the other",
  "Transiently Shared between baby and mother"
)

unique_categories <- c(
  "Persistent Unique to Mother",
  "Transient Unique to Mother",
  "Persistent Unique to Baby",
  "Transient Unique to Baby"
)

# 2. Summarize per Pair
pair_sharing_summary <- df_base %>%
  count(Pair, Interpretation, pattern, name = "count") %>%
  group_by(Pair) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  mutate(
    SharingStatus = case_when(
      Interpretation %in% shared_categories ~ "Shared",
      Interpretation %in% unique_categories ~ "Unique",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Pair, SharingStatus) %>%
  summarise(
    Count = sum(count),
    Prop = sum(prop),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = SharingStatus,
    values_from = c(Count, Prop),
    values_fill = 0
  ) %>%
  mutate(
    Total_vOTUs   = Count_Shared + Count_Unique,
    Shared_Prop   = Prop_Shared,
    Unique_Prop   = Prop_Unique,
    Shared_Count  = Count_Shared,
    Unique_Count  = Count_Unique,
    Total_Prop    = Shared_Prop + Unique_Prop
  ) %>%
  select(Pair, Total_vOTUs, Shared_Count, Unique_Count,
         Shared_Prop, Unique_Prop, Total_Prop) %>%
  mutate(
    Shared_Prop = round(Shared_Prop, 3),
    Unique_Prop = round(Unique_Prop, 3),
    Total_Prop  = round(Total_Prop, 3)
  ) %>%
  arrange(desc(Shared_Prop))

# 3. View concise result
pair_sharing_summary

write.csv(pair_sharing_summary, "Pair_Sharing_Summary.csv", row.names = FALSE)

