library(dplyr)
library(purrr)

# Get tumor types from names that match "_ArmCNV"
tumor_types <- names(data_table)[grepl("_ArmCNV$", names(data_table))] %>%
  gsub("_ArmCNV", "", .)  # Extract tumor type names

# Combine everything
all_scores <- map_dfr(tumor_types, function(tumor) {
  arm_cnv <- data_table[[paste0(tumor, "_ArmCNV")]]
  is_data <- data_table[[paste0(tumor, "_IS")]]

  # Only proceed if both are present
  if (is.null(arm_cnv) || is.null(is_data)) return(NULL)

  # Get shared patients
  shared_samples <- intersect(rownames(arm_cnv), rownames(is_data))

  # Filter to shared patients
  arm_cnv <- arm_cnv[shared_samples, , drop = FALSE]

  # Calculate aneuploidy score
  aneuploidy_score <- rowSums(abs(arm_cnv), na.rm = TRUE)

  # Output data frame
  data.frame(
    Sample = shared_samples,
    Tumor_Type = tumor,
    Aneuploidy_Score = aneuploidy_score,
    stringsAsFactors = FALSE
  )
})

# Save to CSV
write.csv(all_scores, "All_Tumor_Aneuploidy_Scores.csv", row.names = FALSE)



### Plotting for patient 
library(tidyverse)

# Load working table (named list of data frames)
data_table <- readRDS("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/data/TCGA_working_tables.rds")

# Define tumor group mapping
tumor_group_map <- list(
  Adeno = c("LUAD", "LIHC", "ACC", "MESO"),
  GI = c("COADREAD.MSI", "COADREAD.MSS", "ESCA.AD", "STAD", "PAAD"),
  Gyn = c("UCEC.MSI", "UCEC.MSS", "OV", "CESC"),
  Kidney = c("KICH", "KIRC", "KIRP"),
  `NC-Derived` = c("BLCA", "BRCA.neg", "BRCA.pos", "PRAD", "GBM", "LGG", "SKCM", "UVM"),
  Squamous = c("LUSC", "HNSC.HPVneg", "ESCA.SC"),
  Other = c("PCPG", "SARC", "TGCT", "UCS")
)

# Create a flat lookup table: TumorType → TumorGroup
tumor_to_group <- tibble::enframe(tumor_group_map, name = "TumorGroup", value = "TumorType") %>%
  tidyr::unnest(TumorType)

# Loop through each tumor type, extract ArmCNV matrix, calculate altered arms per sample
aneuploidy_df <- purrr::map_dfr(tumor_to_group$TumorType, function(tumor) {
  df_name <- paste0(tumor, "_ArmCNV")
  if (!df_name %in% names(data_table)) return(NULL)
  
  mat <- data_table[[df_name]]
  if (nrow(mat) == 0) return(NULL)

  altered_arm_counts <- rowSums(abs(mat) > 0.2, na.rm = TRUE)

  tibble(
    TumorType = tumor,
    Sample = names(altered_arm_counts),
    AlteredArms = altered_arm_counts
  )
})

# Join tumor group labels
aneuploidy_df <- left_join(aneuploidy_df, tumor_to_group, by = "TumorType")

# Average per tumor group
group_summary <- aneuploidy_df %>%
  group_by(TumorGroup) %>%
  summarize(
    MeanAlteredArms = mean(AlteredArms, na.rm = TRUE),
    SD = sd(AlteredArms, na.rm = TRUE),
    N = n()
  ) %>%
  arrange(desc(MeanAlteredArms))

# Save CSV and plot
write.csv(group_summary, "AneuploidyScore_TumorGroup.csv", row.names = FALSE)

ggplot(group_summary, aes(x = reorder(TumorGroup, MeanAlteredArms), y = MeanAlteredArms)) +
  geom_bar(stat = "identity", fill = "#666666") +
  coord_flip() +
  labs(
    title = "Figure S2: Avg # of Altered Arms per Patient (per Tumor Group)",
    x = "Tumor Group",
    y = "Avg # of Altered Arms"
  ) +
  theme_minimal(base_size = 14)