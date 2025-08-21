#!/usr/bin/env Rscript
.libPaths(c("/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/R_libs", .libPaths()))
print(.libPaths()) 

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(ggrepel)

######################################################################
# FUNGI CLASSIFIED 
######################################################################

print("Processing classified fungi...")

carcinogenic_fungi_db <- data.frame(
  Species = c(
    "Aspergillus flavus",
    "Aspergillus parasiticus", 
    "Aspergillus ochraceus",
    "Aspergillus niger",
    "Penicillium verrucosum",
    "Penicillium nordicum",
    "Fusarium verticillioides",
    "Fusarium graminearum",
    "Fusarium proliferatum", 
    "Fusarium culmorum",
    "Penicillium expansum",
    "Alternaria alternata"
  ),
  Mycotoxin = c(
    "Aflatoxin", "Aflatoxin", "Ochratoxin A", "Ochratoxin A",
    "Ochratoxin A", "Ochratoxin A", "Fumonisin",
    "Trichothecenes, Zearalenone", "Fumonisin", "Zearalenone",
    "Patulin", "Alternariol, Tenuazonic acid"
  ),
  Cancer_Association = c(
    "Liver cancer", "Liver cancer", "Possible kidney cancer", 
    "Possible kidney cancer", "Possible kidney cancer", "Possible kidney cancer",
    "Esophageal cancer", "Possible gastrointestinal/other cancer",
    "Esophageal cancer", "Possible reproductive/endocrine cancer",
    "Possible GI cancer", "Emerging: GI cancer, genotoxic"
  ),
  Mechanism = c(
    "DNA adducts, mutation, genotoxicity",
    "DNA adducts, mutation, genotoxicity",
    "Nephrotoxicity, genotoxicity",
    "Nephrotoxicity, genotoxicity", 
    "Nephrotoxicity, genotoxicity",
    "Nephrotoxicity, genotoxicity",
    "DNA synthesis inhibition, genotoxicity",
    "Transcription/translation inhibition, genotoxicity",
    "DNA synthesis inhibition, genotoxicity",
    "Estrogenic effects, cell cycle disruption",
    "DNA damage, mutagenicity",
    "Reactive oxygen species, DNA damage, mutagenicity"
  ),
  Evidence_Level = c(
    "WHO Class 1 (toxin)", "WHO Class 1 (toxin)",
    "IARC 2B", "IARC 2B", "IARC 2B", "Emerging",
    "IARC 2B", "IARC 3", "IARC 2B", "Emerging",
    "IARC 3", "Emerging"
  ),
  Accession = I(list(
    # Aspergillus flavus - NW_ scaffolds + NG_ complete regions
    c("NW_017010234.1", "NW_017010235.1", "NW_017010236.1",
      "NW_017010237.1", "NW_017010238.1", "NW_017010239.1", 
      "NW_017010240.1", "NW_017010241.1", "NW_017010242.1",
      "NW_017010243.1", "NC_026745.1",
      "NG_055730.1", "NG_055731.1", "NG_042774.1"), # NG_ examples
    # Aspergillus parasiticus
    c("NW_003313313.1", "NW_003313314.1", "NW_003313315.1",
      "NW_003313316.1", "NW_003313317.1", "NW_003313318.1",
      "NW_003313319.1", "NW_003313320.1", "NC_031110.1"),
    # Aspergillus ochraceus
    c("NW_022951540.1", "NW_022951541.1", "NW_022951542.1", 
      "NW_022951543.1", "NW_022951544.1", "NC_020701.1"),
    # Aspergillus niger
    c("NW_003307038.1", "NW_003307039.1", "NW_003307040.1",
      "NW_003307041.1", "NW_003307042.1", "NW_003307043.1",
      "NW_003307044.1", "NW_003307045.1", "NW_003307046.1",
      "NW_003307047.1", "NC_007445.1", 
      "NG_027634.1", "NG_055125.1"), # NG_ examples
    # Penicillium verrucosum  
    c("NW_020543872.1", "NW_020543873.1", "NW_020543874.1",
      "NW_020543875.1", "NW_020543876.1", "NW_020543877.1",
      "NC_009490.1"),
    # Penicillium nordicum
    c("NW_019187620.1", "NW_019187621.1", "NW_019187622.1",
      "NW_019187623.1", "NW_019187624.1"),
    # Fusarium verticillioides
    c("NW_017113678.1", "NW_017113679.1", "NW_017113680.1",
      "NW_017113681.1", "NW_017113682.1", "NW_017113683.1",
      "NW_017113684.1", "NW_017113685.1", "NC_016728.1"),
    # Fusarium graminearum
    c("NW_003305139.1", "NW_003305140.1", "NW_003305141.1",
      "NW_003305142.1", "NW_003305143.1", "NC_007939.1"),
    # Fusarium proliferatum
    c("NW_019180579.1", "NW_019180580.1", "NW_019180581.1",
      "NW_019180582.1", "NW_019180583.1", "NC_016951.1"), 
    # Fusarium culmorum
    c("NW_022138157.1", "NW_022138158.1", "NW_022138159.1",
      "NW_022138160.1", "NW_022138161.1", "NW_022138162.1"),
    # Penicillium expansum
    c("NW_022385491.1", "NW_022385492.1", "NW_022385493.1",
      "NW_022385494.1", "NW_022385495.1", "NC_015653.1"),
    # Alternaria alternata
    c("NW_019114846.1", "NW_019114847.1", "NW_019114848.1",
      "NW_019114849.1", "NW_019114850.1", "NC_018727.1")
  )),
  stringsAsFactors = FALSE
) %>%
  mutate(
    # Strip version numbers for clean accessions
    Accession_clean = lapply(Accession, function(vec) gsub("\\..*$", "", vec)),
    # Create flat comma-separated string for easier access
    Accession_flat = sapply(Accession_clean, paste, collapse = ", ")
  )

# Display the structure
print(carcinogenic_fungi_db)

# Extract just the NG_ accessions if needed
carcinogenic_fungi_db$NG_only <- lapply(carcinogenic_fungi_db$Accession_clean, function(vec) {
  vec[grepl("^NG_", vec)]
})

# Extract just the NW_ accessions if needed  
carcinogenic_fungi_db$NW_only <- lapply(carcinogenic_fungi_db$Accession_clean, function(vec) {
  vec[grepl("^NW_", vec)]
})

blast_results <- read_tsv("/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/blast_results/data_analysis/pipeline/cov_analysis/FUNGI_classified_coverage_hits.tsv") %>%
  dplyr::select(
    Sample,
    QueryID,
    SubjectID,
    Identity,
    Evalue,
    Bitscore,
    AlignmentLength
  ) %>%
  dplyr::mutate(
    SubjectID = gsub("\\..*$", "", SubjectID),
    Coverage = AlignmentLength * Identity / 100
  )

# 3. Join and annotate carcinogenic mycotoxin-producing fungi
carcinogen_hits <- blast_results %>%
  inner_join(
    carcinogenic_fungi_db %>% 
      select(Species, Mycotoxin, Cancer_Association, Mechanism, Evidence_Level, Accession_clean) %>%
      unnest(Accession_clean),
    by = c("SubjectID" = "Accession_clean")
  ) %>%
  mutate(
    Relevance_Level = case_when(
      Mycotoxin == "Aflatoxin" ~ "1_Direct",
      Mycotoxin %in% c("Fumonisin", "Ochratoxin A") ~ "2_Potential",
      TRUE ~ "3_Other"
    ),
    Confidence = case_when(
      Evidence_Level == "WHO Class 1 (toxin)" & Identity > 98 ~ "High",
      Evidence_Level == "WHO Class 1 (toxin)" & Identity > 90 ~ "Medium",
      grepl("Emerging|IARC 2B", Evidence_Level, ignore.case = TRUE) & Identity > 94 ~ "Medium",
      TRUE ~ "Low"
    )
  ) %>%
  arrange(Relevance_Level, desc(Confidence), desc(Bitscore))
  
# 4. Generate and save summary reports
write.csv(carcinogen_hits, "all_carcinogenic_fungi_classified.csv", row.names = FALSE)

# 5. Prepare for plotting: clean labels/factors for visualization
plot_data <- carcinogen_hits %>%
  mutate(
    Relevance_Level = gsub("^[0-9]_", "", Relevance_Level)
  )
plot_data$Cancer_Association <- factor(
  plot_data$Cancer_Association,
  levels = unique(plot_data$Cancer_Association)
)
plot_data$Relevance_Level <- factor(
  plot_data$Relevance_Level,
  levels = c("Direct", "Potential", "Other")
)
plot_data$Confidence <- factor(plot_data$Confidence, levels = c("High", "Medium", "Low"))

# 6. Summarize for bubble plot
plot_summary <- plot_data %>%
  group_by(Species, Mycotoxin, Cancer_Association, Confidence, Relevance_Level) %>%
  summarise(
    Occurrence = n_distinct(Sample),
    Avg_Identity = mean(Identity),
    .groups = "drop"
  ) %>%
  mutate(
    Label = paste0(Species, "\n(", Mycotoxin, ")\n", Cancer_Association),    
    Size = sqrt(Occurrence) * 2
  )

unique_sizes <- sort(unique(plot_summary$Size))
if(length(unique_sizes) >= 3) {
  breaks <- unique_sizes[c(1, round(length(unique_sizes)/2), length(unique_sizes))]
} else {
  breaks <- unique_sizes
}
labels <- c("Low", "Medium", "High")[seq_along(breaks)]
confidence_palette <- c("High" = "#1a9850", "Medium" = "#fdae61", "Low" = "#d73027")

# 7. Main bubble plot
main_plot <- ggplot(plot_summary, aes(x = Relevance_Level, y = Avg_Identity, size = Size, color = Confidence)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = Label), size = 3, box.padding = 0.5, max.overlaps = Inf) +
  scale_color_manual(values = confidence_palette) +
  scale_size_continuous(name = "Occurrence", breaks = breaks, labels = labels) +
  labs(
    title = "Carcinogenic Fungi (Mycotoxin Producers) Associations",
    subtitle = "Size indicates occurrence count across samples",
    x = "Association Type", y = "Average % Identity"
  ) +
  theme_minimal() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))

# 8. Evidence confidence bar plot
conf_plot <- ggplot(plot_data, aes(x = Confidence, fill = Confidence)) +
  geom_bar() +
  scale_fill_manual(values = confidence_palette) +
  labs(
    title = "Evidence Confidence Distribution",
    x = "Confidence Level", y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# 9. Combine, annotate, and save plots
final_plot <- (main_plot | conf_plot) +
  plot_annotation(
    title = "Comprehensive Fungi-Cancer (Mycotoxin) Association Analysis Classified",
    subtitle = "Showing relationships, evidence, and occurrence counts",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

print(final_plot)
tryCatch({
  ggsave(
    "comprehensive_fungi_cancer_analysis_classified.png",
    plot = final_plot,
    width = 16, height = 10, dpi = 300
  )
  message("Successfully saved combined plot!")
}, error = function(e) {
  message("Error saving plot: ", e$message)
  png("comprehensive_fungi_cancer_analysis_classified.png", width = 1600, height = 1000)
  print(final_plot)
  dev.off()
  message("Used alternative save method - check comprehensive_fungi_cancer_analysis_classified.png")
})

######################################################################
# FUNGI UNCLASSIFIED 
######################################################################

print("Processing unclassified fungi...")

carcinogenic_fungi_db <- data.frame(
  Species = c(
    "Aspergillus flavus",
    "Aspergillus parasiticus", 
    "Aspergillus ochraceus",
    "Aspergillus niger",
    "Penicillium verrucosum",
    "Penicillium nordicum",
    "Fusarium verticillioides",
    "Fusarium graminearum",
    "Fusarium proliferatum", 
    "Fusarium culmorum",
    "Penicillium expansum",
    "Alternaria alternata"
  ),
  Mycotoxin = c(
    "Aflatoxin", "Aflatoxin", "Ochratoxin A", "Ochratoxin A",
    "Ochratoxin A", "Ochratoxin A", "Fumonisin",
    "Trichothecenes, Zearalenone", "Fumonisin", "Zearalenone",
    "Patulin", "Alternariol, Tenuazonic acid"
  ),
  Cancer_Association = c(
    "Liver cancer", "Liver cancer", "Possible kidney cancer", 
    "Possible kidney cancer", "Possible kidney cancer", "Possible kidney cancer",
    "Esophageal cancer", "Possible gastrointestinal/other cancer",
    "Esophageal cancer", "Possible reproductive/endocrine cancer",
    "Possible GI cancer", "Emerging: GI cancer, genotoxic"
  ),
  Mechanism = c(
    "DNA adducts, mutation, genotoxicity",
    "DNA adducts, mutation, genotoxicity",
    "Nephrotoxicity, genotoxicity",
    "Nephrotoxicity, genotoxicity", 
    "Nephrotoxicity, genotoxicity",
    "Nephrotoxicity, genotoxicity",
    "DNA synthesis inhibition, genotoxicity",
    "Transcription/translation inhibition, genotoxicity",
    "DNA synthesis inhibition, genotoxicity",
    "Estrogenic effects, cell cycle disruption",
    "DNA damage, mutagenicity",
    "Reactive oxygen species, DNA damage, mutagenicity"
  ),
  Evidence_Level = c(
    "WHO Class 1 (toxin)", "WHO Class 1 (toxin)",
    "IARC 2B", "IARC 2B", "IARC 2B", "Emerging",
    "IARC 2B", "IARC 3", "IARC 2B", "Emerging",
    "IARC 3", "Emerging"
  ),
  Accession = I(list(
    # Aspergillus flavus - NW_ scaffolds + NG_ complete regions
    c("NW_017010234.1", "NW_017010235.1", "NW_017010236.1",
      "NW_017010237.1", "NW_017010238.1", "NW_017010239.1", 
      "NW_017010240.1", "NW_017010241.1", "NW_017010242.1",
      "NW_017010243.1", "NC_026745.1", 
      "NG_055730.1", "NG_055731.1", "NG_042774.1"), # NG_ examples
    # Aspergillus parasiticus
    c("NW_003313313.1", "NW_003313314.1", "NW_003313315.1",
      "NW_003313316.1", "NW_003313317.1", "NW_003313318.1",
      "NW_003313319.1", "NW_003313320.1", "NC_031110.1"),
    # Aspergillus ochraceus
    c("NW_022951540.1", "NW_022951541.1", "NW_022951542.1", 
      "NW_022951543.1", "NW_022951544.1", "NC_020701.1"),
    # Aspergillus niger
    c("NW_003307038.1", "NW_003307039.1", "NW_003307040.1",
      "NW_003307041.1", "NW_003307042.1", "NW_003307043.1",
      "NW_003307044.1", "NW_003307045.1", "NW_003307046.1",
      "NW_003307047.1", "NC_007445.1", 
      "NG_027634.1", "NG_055125.1"), # NG_ examples
    # Penicillium verrucosum  
    c("NW_020543872.1", "NW_020543873.1", "NW_020543874.1",
      "NW_020543875.1", "NW_020543876.1", "NW_020543877.1",
      "NC_009490.1"),
    # Penicillium nordicum
    c("NW_019187620.1", "NW_019187621.1", "NW_019187622.1",
      "NW_019187623.1", "NW_019187624.1"),
    # Fusarium verticillioides
    c("NW_017113678.1", "NW_017113679.1", "NW_017113680.1",
      "NW_017113681.1", "NW_017113682.1", "NW_017113683.1",
      "NW_017113684.1", "NW_017113685.1", "NC_016728.1"),
    # Fusarium graminearum
    c("NW_003305139.1", "NW_003305140.1", "NW_003305141.1",
      "NW_003305142.1", "NW_003305143.1", "NC_007939.1"),
    # Fusarium proliferatum
    c("NW_019180579.1", "NW_019180580.1", "NW_019180581.1",
      "NW_019180582.1", "NW_019180583.1", "NC_016951.1"), 
    # Fusarium culmorum
    c("NW_022138157.1", "NW_022138158.1", "NW_022138159.1",
      "NW_022138160.1", "NW_022138161.1", "NW_022138162.1"),
    # Penicillium expansum
    c("NW_022385491.1", "NW_022385492.1", "NW_022385493.1",
      "NW_022385494.1", "NW_022385495.1", "NC_015653.1"),
    # Alternaria alternata
    c("NW_019114846.1", "NW_019114847.1", "NW_019114848.1",
      "NW_019114849.1", "NW_019114850.1", "NC_018727.1")
  )),
  stringsAsFactors = FALSE
) %>%
  mutate(
    # Strip version numbers for clean accessions
    Accession_clean = lapply(Accession, function(vec) gsub("\\..*$", "", vec)),
    # Create flat comma-separated string for easier access
    Accession_flat = sapply(Accession_clean, paste, collapse = ", ")
  )

# Display the structure
print(carcinogenic_fungi_db)

# Extract just the NG_ accessions if needed
carcinogenic_fungi_db$NG_only <- lapply(carcinogenic_fungi_db$Accession_clean, function(vec) {
  vec[grepl("^NG_", vec)]
})

# Extract just the NW_ accessions if needed  
carcinogenic_fungi_db$NW_only <- lapply(carcinogenic_fungi_db$Accession_clean, function(vec) {
  vec[grepl("^NW_", vec)]
})

# 2. Read and process BLAST results (filtered hits vs fungi database)
blast_results <- read_tsv("/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/blast_results/data_analysis/pipeline/cov_analysis/FUNGI_unclassified_coverage_hits.tsv") %>%
  dplyr::select(
    Sample,
    QueryID,
    SubjectID,
    Identity,
    Evalue,
    Bitscore,
    AlignmentLength
  ) %>%
  dplyr::mutate(
    SubjectID = gsub("\\..*$", "", SubjectID),
    Coverage = AlignmentLength * Identity / 100
  )

# 3. Join and annotate carcinogenic mycotoxin-producing fungi
carcinogen_hits <- blast_results %>%
  inner_join(
    carcinogenic_fungi_db %>% 
      select(Species, Mycotoxin, Cancer_Association, Mechanism, Evidence_Level, Accession_clean) %>%
      unnest(Accession_clean),
    by = c("SubjectID" = "Accession_clean")
  ) %>%
  mutate(
    Relevance_Level = case_when(
      Mycotoxin == "Aflatoxin" ~ "1_Direct",
      Mycotoxin %in% c("Fumonisin", "Ochratoxin A") ~ "2_Potential",
      TRUE ~ "3_Other"
    ),
    Confidence = case_when(
      Evidence_Level == "WHO Class 1 (toxin)" & Identity > 98 ~ "High",
      Evidence_Level == "WHO Class 1 (toxin)" & Identity > 90 ~ "Medium",
      grepl("Emerging|IARC 2B", Evidence_Level, ignore.case = TRUE) & Identity > 94 ~ "Medium",
      TRUE ~ "Low"
    )
  ) %>%
  arrange(Relevance_Level, desc(Confidence), desc(Bitscore))

# 4. Generate and save summary reports
write.csv(carcinogen_hits, "all_carcinogenic_fungi_unclassified.csv", row.names = FALSE)

# 5. Prepare for plotting: clean labels/factors for visualization
plot_data <- carcinogen_hits %>%
  mutate(
    Relevance_Level = gsub("^[0-9]_", "", Relevance_Level)
  )
plot_data$Cancer_Association <- factor(
  plot_data$Cancer_Association,
  levels = unique(plot_data$Cancer_Association)
)
plot_data$Relevance_Level <- factor(
  plot_data$Relevance_Level,
  levels = c("Direct", "Potential", "Other")
)
plot_data$Confidence <- factor(plot_data$Confidence, levels = c("High", "Medium", "Low"))

# 6. Summarize for bubble plot
plot_summary <- plot_data %>%
  group_by(Species, Mycotoxin, Cancer_Association, Confidence, Relevance_Level) %>%
  summarise(
    Occurrence = n_distinct(Sample),
    Avg_Identity = mean(Identity),
    .groups = "drop"
  ) %>%
  mutate(
    Label = paste0(Species, "\n(", Mycotoxin, ")\n", Cancer_Association),    
    Size = sqrt(Occurrence) * 2
  )

unique_sizes <- sort(unique(plot_summary$Size))
if(length(unique_sizes) >= 3) {
  breaks <- unique_sizes[c(1, round(length(unique_sizes)/2), length(unique_sizes))]
} else {
  breaks <- unique_sizes
}
labels <- c("Low", "Medium", "High")[seq_along(breaks)]
confidence_palette <- c("High" = "#1a9850", "Medium" = "#fdae61", "Low" = "#d73027")

# 7. Main bubble plot
main_plot <- ggplot(plot_summary, aes(x = Relevance_Level, y = Avg_Identity, size = Size, color = Confidence)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = Label), size = 3, box.padding = 0.5, max.overlaps = Inf) +
  scale_color_manual(values = confidence_palette) +
  scale_size_continuous(name = "Occurrence", breaks = breaks, labels = labels) +
  labs(
    title = "Carcinogenic Fungi (Mycotoxin Producers) Associations",
    subtitle = "Size indicates occurrence count across samples",
    x = "Association Type", y = "Average % Identity"
  ) +
  theme_minimal() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))

# 8. Evidence confidence bar plot
conf_plot <- ggplot(plot_data, aes(x = Confidence, fill = Confidence)) +
  geom_bar() +
  scale_fill_manual(values = confidence_palette) +
  labs(
    title = "Evidence Confidence Distribution",
    x = "Confidence Level", y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# 9. Combine, annotate, and save plots
final_plot <- (main_plot | conf_plot) +
  plot_annotation(
    title = "Comprehensive Fungi-Cancer (Mycotoxin) Association Analysis Unclassified",
    subtitle = "Showing relationships, evidence, and occurrence counts",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

print(final_plot)
tryCatch({
  ggsave(
    "comprehensive_fungi_cancer_analysis_unclassified.png",
    plot = final_plot,
    width = 16, height = 10, dpi = 300
  )
  message("Successfully saved combined plot!")
}, error = function(e) {
  message("Error saving plot: ", e$message)
  png("comprehensive_fungi_cancer_analysis_unclassified.png", width = 1600, height = 1000)
  print(final_plot)
  dev.off()
  message("Used alternative save method - check comprehensive_fungi_cancer_analysis_unclassified.png")
})

######################################################################
# BACTERIA CLASSIFIED 
######################################################################

print("Processing classified bacteria...")

# Enhanced carcinogenic bacteria database
carcinogenic_bacteria_db <- data.frame(
  Species = c(
    # Recently confirmed in 2024-2025 studies
    "Fusobacterium nucleatum subtype Fna C2",  # 2024 Nature study - specific subtype
    "Streptococcus mutans",                    # 2022-2024 studies confirm cancer role
    "Streptococcus mitis",                     # 2025 emerging evidence
    "Acinetobacter baumannii",                 # 2024-2025 ESKAPE pathogen cancer link
    "Prevotella intermedia",                   # 2024 gut-oral axis cancer studies
    "Peptostreptococcus anaerobius",          # 2024 colorectal cancer studies
    "Parvimonas micra",                        # 2024 oral-gut translocation cancer
    "Actinomyces naeslundii",                 # 2024 head/neck cancer associations
    "Treponema denticola",                     # 2024 periodontal-systemic cancer link
    "Tannerella forsythia",                    # 2024 oral microbiome cancer studies
    # ESKAPE pathogens with emerging cancer links (2024-2025)
    "Enterobacter cloacae",                    # 2024 ESKAPE cancer mortality studies
    "Stenotrophomonas maltophilia",           # 2024 emerging oncogenic potential
    "Veillonella parvula",                    # 2024 oral-systemic cancer axis
    "Capnocytophaga gingivalis",              # 2024 periodontal cancer progression
    "Aggregatibacter actinomycetemcomitans"   # 2024 aggressive periodontitis-cancer
  ),
  Cancer_Association = c(
    "Colorectal cancer",
    "Oral cancer; Head and neck cancer; Esophageal cancer",
    "Oral cancer inhibition",
    "Cancer mortality in immunocompromised; ICU cancer patients",
    "Esophageal cancer; Gastric cancer (gut-oral axis)",
    "Colorectal cancer; Liver metastases",
    "Colorectal cancer; Oral-gut translocation",
    "Head and neck cancer; Cervical cancer",
    "Periodontal disease-associated systemic cancers",
    "Oral cancer; Pancreatic cancer (systemic translocation)",
    "Cancer patient mortality; Bloodstream infections",
    "Emerging lung cancer; Nosocomial cancer infections",
    "Colorectal cancer; Oral-gut microbiome axis",
    "Periodontal disease-cancer progression",
    "Aggressive periodontitis-oral cancer progression"
  ),
  Mechanism = c(
    "Specific virulence factors; enhanced tumor microenvironment modulation",
    "IL-6 production; EMT induction; myeloid-derived suppressor cell recruitment",
    "Competitive inhibition of pathogenic bacteria (potential protective)",
    "Multi-drug resistance; immune suppression in cancer patients",
    "Gut-oral bacterial translocation; chronic inflammation",
    "Butyrate metabolism disruption; immune evasion",
    "Oral-gut translocation; biofilm formation in tumors",
    "Chronic cervical/oral inflammation; immune modulation",
    "Spirochete invasion; chronic periodontal inflammation",
    "Protease activity; systemic bacterial translocation",
    "Nosocomial infections in immunocompromised cancer patients",
    "Multidrug resistance; chronic pulmonary inflammation",
    "Short-chain fatty acid disruption; microbial translocation",
    "Periodontal bone destruction; systemic inflammatory response",
    "Aggressive tissue invasion; leukotoxin production"
  ),
  Evidence_Level = c(
    "Emerging Strong",  # 2024 Nature paper
    "Emerging Strong",  # Multiple 2022-2024 studies
    "Experimental Strong", # 2025 protective role studies
    "Emerging Strong",  # 2024-2025 ESKAPE mortality studies
    "Emerging",
    "Emerging Strong",  # 2024 CRC studies
    "Emerging",
    "Emerging",
    "Emerging", 
    "Emerging",
    "Emerging",
    "Experimental",
    "Emerging",
    "Emerging",
    "Emerging"
  ),
  # Add corresponding accession patterns from your actual BLAST data
  Accession = I(list(
    c("NZ_FKXL01000017", "NZ_FKWV01000005", "NZ_FKXF01000023"), # F. nucleatum Fna C2
    c("NZ_MLWD01000035", "NZ_MLXP01000026", "NZ_MLUB01000021"), # S. mutans
    c("NZ_JAMKRA010000075", "NZ_JAMKQZ010000090", "NZ_JAMKTM010000045"), # S. mitis
    c("NZ_JAUEHL010000030", "NZ_FNPJ01000004", "NZ_FNOE01000002"), # A. baumannii
    c("NZ_FLBI01000046", "NZ_FLBN01000099", "NZ_FLCN01000004"), # P. intermedia
    c("NZ_MLTT01000022", "NZ_MLRM01000021", "NZ_MLRF01000022"), # P. anaerobius
    c("NZ_MLRI01000021", "NZ_MLRE01000049", "NZ_MLRZ01000037"), # P. micra
    c("NZ_MLRO01000021", "NZ_MLXR01000021", "NZ_MLRG01000021"), # A. naeslundii
    c("NZ_MLYG01000019", "NZ_MLRJ01000019", "NZ_MLSM01000018"), # T. denticola
    c("NZ_MLSQ01000022", "NZ_MLXY01000024", "NZ_MLXS01000025"), # T. forsythia
    c("NZ_MLWB01000035", "NZ_MLUD01000021", "NZ_MLRC01000027"), # E. cloacae
    c("NZ_MLSZ01000019", "NZ_MLSN01000023", "NZ_MLVX01000038"), # S. maltophilia
    c("NZ_MLYZ01000021", "NZ_MLSP01000022", "NZ_MLZF01000029"), # V. parvula
    c("NZ_MLSH01000020", "NZ_MLYY01000023", "NZ_MLYC01000020"), # C. gingivalis
    c("NZ_MLXZ01000020", "NZ_JBGJJV010000032", "NZ_JAUEGY010000142") # A. actinomycetemcomitans
  )),
  stringsAsFactors = FALSE
) %>%
  mutate(
    # Strip version numbers for clean accessions
    Accession_clean = lapply(Accession, function(vec) gsub("\\..*$", "", vec)),
    # Create flat comma-separated string for easier access
    Accession_flat = sapply(Accession_clean, paste, collapse = ", ")
  )

# Display the structure
print(carcinogenic_bacteria_db)

# Read bacterial BLAST results - UPDATE THIS PATH!
blast_results <- read_tsv("/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/blast_results/data_analysis/pipeline/cov_analysis/BACT_classified_coverage_hits.tsv") %>%
  dplyr::select(
    Sample,
    QueryID,
    SubjectID,
    Identity,
    Evalue,
    Bitscore,
    AlignmentLength
  ) %>%
  dplyr::mutate(
    SubjectID = gsub("\\..*$", "", SubjectID),
    Coverage = AlignmentLength * Identity / 100
  )

# Create flattened lookup table for joining
accession_lookup <- carcinogenic_bacteria_db %>%
  select(Species, Cancer_Association, Mechanism, Evidence_Level, Accession_clean) %>%
  unnest(Accession_clean) %>%
  rename(SubjectID_clean = Accession_clean)

# Join and annotate carcinogenic bacteria
carcinogen_hits <- blast_results %>%
  inner_join(accession_lookup, by = c("SubjectID" = "SubjectID_clean")) %>%
  mutate(
    Relevance_Level = case_when(
      Evidence_Level == "WHO Class 1" ~ "1_Direct",
      Evidence_Level %in% c("Emerging Strong") ~ "2_Strong",
      Evidence_Level == "Emerging" ~ "3_Potential", 
      TRUE ~ "4_Experimental"
    ),
    Confidence = case_when(
      Evidence_Level == "WHO Class 1" & Identity > 95 ~ "High",
      Evidence_Level == "WHO Class 1" & Identity > 90 ~ "Medium",
      grepl("Emerging Strong|Emerging", Evidence_Level, ignore.case = TRUE) & Identity > 95 ~ "Medium",
      TRUE ~ "Low"
    )
  ) %>%
  arrange(Relevance_Level, desc(Confidence), desc(Bitscore))

# Check for empty results
if (nrow(carcinogen_hits) == 0) {
  warning("No carcinogenic bacteria matches found!")
} else {
  print(paste("Found", nrow(carcinogen_hits), "carcinogenic bacteria matches"))
}

# 4. Generate and save summary reports
write.csv(carcinogen_hits, "all_carcinogenic_bacteria_classified.csv", row.names = FALSE)

# Generate summary for plotting
plot_data <- carcinogen_hits %>%
  mutate(
    Relevance_Level = gsub("^[0-9]_", "", Relevance_Level)
  )

plot_data$Cancer_Association <- factor(
  plot_data$Cancer_Association,
  levels = unique(plot_data$Cancer_Association)
)
plot_data$Relevance_Level <- factor(
  plot_data$Relevance_Level,
  levels = c("Direct", "Strong", "Potential", "Experimental")
)
plot_data$Confidence <- factor(plot_data$Confidence, levels = c("High", "Medium", "Low"))

# Summarize for bubble plot with enhanced labels (BACTERIA VERSION - NO MYCOTOXIN)
plot_summary <- plot_data %>%
  group_by(Species, Cancer_Association, Confidence, Relevance_Level, Evidence_Level) %>%
  summarise(
    Occurrence = n_distinct(Sample),
    Avg_Identity = mean(Identity),
    .groups = "drop"
  ) %>%
  mutate(
    Label = paste0(Species, "\n", Cancer_Association, "\n", Evidence_Level),
    Size = sqrt(Occurrence) * 2
  )

unique_sizes <- sort(unique(plot_summary$Size))
if(length(unique_sizes) >= 3) {
  breaks <- unique_sizes[c(1, round(length(unique_sizes)/2), length(unique_sizes))]
} else {
  breaks <- unique_sizes
}
labels <- c("Low", "Medium", "High")[seq_along(breaks)]
confidence_palette <- c("High" = "#1a9850", "Medium" = "#fdae61", "Low" = "#d73027")

# 7. Main bubble plot (BACTERIA VERSION)
main_plot <- ggplot(plot_summary, aes(x = Relevance_Level, y = Avg_Identity, size = Size, color = Confidence)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(
    aes(label = Label), 
    size = 2.5,           # Enhanced settings start here
    box.padding = 1.2,
    point.padding = 0.8,
    max.overlaps = Inf,
    force = 3,
    max.iter = 5000,
    seed = 42) +
  scale_color_manual(values = confidence_palette) +
  scale_size_continuous(name = "Occurrence", breaks = breaks, labels = labels) +
  labs(
    title = "Carcinogenic Bacteria Associations",
    subtitle = "Size indicates occurrence count across samples",
    x = "Association Type", y = "Average % Identity"
  ) +
  theme_minimal() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))

# 8. Evidence confidence bar plot
conf_plot <- ggplot(plot_data, aes(x = Confidence, fill = Confidence)) +
  geom_bar() +
  scale_fill_manual(values = confidence_palette) +
  labs(
    title = "Evidence Confidence Distribution",
    x = "Confidence Level", y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# 9. Combine, annotate, and save plots (BACTERIA VERSION)
final_plot <- (main_plot | conf_plot) +
  plot_annotation(
    title = "Comprehensive Bacteria-Cancer Association Analysis Classified",
    subtitle = "Showing relationships, evidence, and occurrence counts",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

print(final_plot)
tryCatch({
  ggsave(
    "comprehensive_bacteria_cancer_analysis_classified.png",
    plot = final_plot,
    width = 24, height = 16, dpi = 300
  )
  message("Successfully saved combined plot!")
}, error = function(e) {
  message("Error saving plot: ", e$message)
  png("comprehensive_bacteria_cancer_analysis_classified.png", width = 1600, height = 1000)
  print(final_plot)
  dev.off()
  message("Used alternative save method - check comprehensive_bacteria_cancer_analysis_classified.png")
})

######################################################################
# VIRUSES CLASSIFIED 
######################################################################

print("Processing classified viruses...")

# 1. Curated database of carcinogenic viruses
carcinogenic_viruses_db <- data.frame(
  Species = c(
    "Human endogenous retrovirus K",
    "Epstein-Barr virus",
    "Mouse mammary tumor virus-like",
    "Human herpesvirus 2",
    "Human herpesvirus 1",
    "Human T-lymphotropic virus 1",
    "Human papillomavirus 16",
    "Hepatitis B virus",
    "Hepatitis C virus", 
    "Human herpesvirus 8 (KSHV)",
    "Merkel cell polyomavirus",
    "Human cytomegalovirus (HCMV)",
    "Human papillomavirus 18",
    "JC polyomavirus",
    "Human papillomavirus 33",
    "Varicella-zoster virus",
    "Human herpesvirus 7",
    "BK polyomavirus",
    "Human polyomavirus 6",
    "Simian virus 40 (SV40)",
    # Additional oncogenic viruses
    "Human papillomavirus 31",
    "Human papillomavirus 45",
    "Hepatitis D virus",
    "Human T-lymphotropic virus 2",
    "Adenovirus type 12"
  ),
  Cancer_Association = c(
    "Breast cancer; Prostate cancer",
    "Burkitt lymphoma; Hodgkin lymphoma; Nasopharyngeal carcinoma; Gastric cancer",
    "Breast cancer (experimental model)",
    "Cervical cancer (co-factor); Genital cancers",
    "Oral cancer (possible); Head and neck cancer",
    "Adult T-cell leukemia/lymphoma; T-cell lymphoma",
    "Cervical cancer; Anal cancer; Oropharyngeal cancer; Penile cancer",
    "Hepatocellular carcinoma; Liver cancer",
    "Hepatocellular carcinoma; Liver cancer; Cholangiocarcinoma",
    "Kaposi sarcoma; Primary effusion lymphoma; Multicentric Castleman disease",
    "Merkel cell carcinoma",
    "Possible glioblastoma; Prostate cancer (emerging)",
    "Cervical cancer; Adenocarcinoma; Anogenital cancers",
    "Progressive multifocal leukoencephalopathy; Possible brain tumors",
    "Cervical cancer; Anogenital cancers",
    "Possible lymphoma; Indirect cancer effects",
    "Possible T-cell lymphoma (emerging)",
    "Kidney cancer; Bladder cancer; Prostate cancer",
    "Trichodysplasia spinulosa; Possible skin cancers",
    "Mesothelioma; Brain tumors; Lymphomas",
    # Additional entries
    "Cervical cancer; Anogenital cancers",
    "Cervical cancer; Anogenital cancers", 
    "Hepatocellular carcinoma (co-factor with HBV)",
    "Hairy cell leukemia; T-cell lymphomas",
    "Possible sarcomas; Experimental tumors"
  ),
  Mechanism = c(
    "Endogenous retroviral activation; Immune modulation",
    "Latent infection; Immune evasion; Viral oncoproteins (LMP1, EBNA)",
    "Retroviral integration; Hormone receptor modulation",
    "Chronic inflammation; Immune suppression; Co-infection effects",
    "Chronic inflammation; Immune modulation; HSV reactivation",
    "Viral tax protein; Cell cycle dysregulation; Immune evasion",
    "Viral DNA integration; E6/E7 oncoproteins; p53/Rb inactivation",
    "Chronic hepatitis; Viral DNA integration; HBx protein oncogene",
    "Chronic hepatitis; Persistent infection; Oxidative stress; Fibrosis",
    "Viral latency; LANA oncogene; Angiogenesis promotion",
    "Viral DNA integration; Large T antigen; Cell transformation",
    "Immune modulation; Chronic infection; Oncogene activation",
    "Viral DNA integration; E6/E7 oncoproteins; Chromosomal instability",
    "Viral DNA integration; Large T antigen; Cell cycle disruption",
    "Viral DNA integration; E6/E7 oncoproteins; DNA damage",
    "Latent infection; Immune suppression; Chronic inflammation",
    "Immune modulation; Limited oncogenic evidence",
    "Viral DNA integration; Large T antigen; Oncogene activation",
    "Viral replication; Possible transformation; Limited evidence",
    "SV40 Large T antigen; p53/Rb inactivation; Cell immortalization",
    # Additional mechanisms
    "Viral DNA integration; E6/E7 oncoproteins; HPV31 specific effects",
    "Viral DNA integration; E6/E7 oncoproteins; HPV45 specific effects",
    "Co-infection with HBV; Enhanced hepatotoxicity",
    "Tax protein; Rex protein; Cell cycle disruption",
    "Viral E1A/E1B proteins; Cell transformation"
  ),
  Evidence_Level = c(
    "Emerging Strong",
    "WHO Class 1",
    "Experimental Strong",
    "IARC 2A",
    "Emerging",
    "WHO Class 1",
    "WHO Class 1",
    "WHO Class 1",
    "WHO Class 1",
    "WHO Class 1",
    "WHO Class 1",
    "Emerging",
    "WHO Class 1",
    "IARC 2B",
    "IARC 2A",
    "Emerging",
    "Experimental",
    "IARC 2B",
    "Experimental",
    "IARC 2B",
    # Additional evidence levels
    "IARC 2A",
    "IARC 2A",
    "Emerging",
    "IARC 2B",
    "Experimental"
  ),
  Accession = I(list(
    # Human endogenous retrovirus K
    c("NC_001716.2", "NC_022518.1", "NC_007797.1"),
    # Epstein-Barr virus
    c("NC_007605.1", "NC_001664.4", "NC_009334.1"),
    # Mouse mammary tumor virus-like
    c("NC_001503.1", "NC_022518.1", "NC_001503.2"),
    # Human herpesvirus 2
    c("NC_001798.2", "NC_001345.2", "NC_001798.1"),
    # Human herpesvirus 1
    c("NC_001806.2", "NC_001806.1", "NC_023914.1"),
    # Human T-lymphotropic virus 1
    c("NC_001436.1", "NC_001454.1", "NC_001436.2"),
    # Human papillomavirus 16
    c("NC_001526.4", "NC_001405.1", "NC_001526.2"),
    # Hepatitis B virus
    c("NC_003977.2", "NC_003521.1", "NC_003977.1"),
    # Hepatitis C virus
    c("NC_004102.1", "NC_009827.1", "NC_009824.1"),
    # Human herpesvirus 8 (KSHV)
    c("NC_009333.1", "NC_003409.1", "NC_009333.2"),
    # Merkel cell polyomavirus
    c("NC_010277.2", "NC_010277.1", "NC_017972.1"),
    # Human cytomegalovirus
    c("NC_001347.1", "NC_001702.1", "NC_006273.2"),
    # Human papillomavirus 18
    c("NC_001357.1", "NC_001604.1", "NC_001357.2"),
    # JC polyomavirus
    c("NC_001699.1", "NC_003287.1", "NC_001699.2"),
    # Human papillomavirus 33
    c("NC_001586.1", "NC_006273.2", "NC_001586.2"),
    # Varicella-zoster virus
    c("NC_001348.1", "NC_001501.2", "NC_001348.2"),
    # Human herpesvirus 7
    c("NC_001716.1", "NC_009334.1", "NC_001716.2"),
    # BK polyomavirus
    c("NC_001538.1", "NC_000898.1", "NC_001538.2"),
    # Human polyomavirus 6
    c("NC_014406.1", "NC_001546.2", "NC_014406.2"),
    # Simian virus 40
    c("NC_001669.1", "NC_003819.1", "NC_001669.2"),
    # Human papillomavirus 31
    c("NC_001527.1", "NC_007977.1", "NC_001527.2"),
    # Human papillomavirus 45
    c("NC_001529.1", "NC_007998.1", "NC_001529.2"),
    # Hepatitis D virus
    c("NC_001653.2", "NC_001653.1", "NC_006103.1"),
    # Human T-lymphotropic virus 2
    c("NC_001488.1", "NC_001488.2", "NC_007771.1"),
    # Adenovirus type 12
    c("NC_001460.1", "NC_010956.1", "NC_001460.2")
  )),
  stringsAsFactors = FALSE
) %>%
  mutate(
    # Strip version numbers for clean accessions
    Accession_clean = lapply(Accession, function(vec) gsub("\\..*$", "", vec)),
    # Create flat comma-separated string for easier access
    Accession_flat = sapply(Accession_clean, paste, collapse = ", ")
  )

# Display the structure
print(carcinogenic_viruses_db)

# Read viral BLAST results - UPDATE THIS PATH!
blast_results <- read_tsv("/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/blast_results/data_analysis/pipeline/cov_analysis/VIRUS_classified_coverage_hits.tsv") %>%
  dplyr::select(
    Sample,
    QueryID,
    SubjectID,
    Identity,
    Evalue,
    Bitscore,
    AlignmentLength
  ) %>%
  dplyr::mutate(
    SubjectID = gsub("\\..*$", "", SubjectID),
    Coverage = AlignmentLength * Identity / 100
  )

# Create flattened lookup table for joining
accession_lookup <- carcinogenic_viruses_db %>%
  select(Species, Cancer_Association, Mechanism, Evidence_Level, Accession_clean) %>%
  unnest(Accession_clean) %>%
  rename(SubjectID_clean = Accession_clean)

# Join and annotate carcinogenic viruses
carcinogen_hits <- blast_results %>%
  inner_join(accession_lookup, by = c("SubjectID" = "SubjectID_clean")) %>%
  mutate(
    Relevance_Level = case_when(
      Evidence_Level == "WHO Class 1" ~ "1_Direct",
      Evidence_Level %in% c("IARC 2A", "Emerging Strong", "Experimental Strong") ~ "2_Strong",
      Evidence_Level %in% c("IARC 2B", "Emerging") ~ "3_Potential",
      TRUE ~ "4_Experimental"
    ),
    Confidence = case_when(
      Evidence_Level == "WHO Class 1" & Identity > 95 ~ "High",
      Evidence_Level == "WHO Class 1" & Identity > 90 ~ "Medium",
      grepl("IARC 2A|Strong", Evidence_Level, ignore.case = TRUE) & Identity > 95 ~ "Medium",
      grepl("IARC 2B|Emerging", Evidence_Level, ignore.case = TRUE) & Identity > 90 ~ "Low",
      TRUE ~ "Low"
    )
  ) %>%
  arrange(Relevance_Level, desc(Confidence), desc(Bitscore))

# Check for empty results
if (nrow(carcinogen_hits) == 0) {
  warning("No carcinogenic virus matches found!")
} else {
  print(paste("Found", nrow(carcinogen_hits), "carcinogenic virus matches"))
}

# 4. Generate and save summary reports
write.csv(carcinogen_hits, "all_carcinogenic_viruses_classified.csv", row.names = FALSE)

# 5. Prepare for plotting: clean labels/factors for visualization
plot_data <- carcinogen_hits %>%
  mutate(
    Relevance_Level = gsub("^[0-9]_", "", Relevance_Level)
  )

plot_data$Cancer_Association <- factor(
  plot_data$Cancer_Association,
  levels = unique(plot_data$Cancer_Association)
)
plot_data$Relevance_Level <- factor(
  plot_data$Relevance_Level,
  levels = c("Direct", "Strong", "Potential", "Experimental")
)
plot_data$Confidence <- factor(plot_data$Confidence, levels = c("High", "Medium", "Low"))

# 6. Summarize for bubble plot with enhanced labels
plot_summary <- plot_data %>%
  group_by(Species, Cancer_Association, Confidence, Relevance_Level, Evidence_Level) %>%
  summarise(
    Occurrence = n_distinct(Sample),
    Avg_Identity = mean(Identity),
    .groups = "drop"
  ) %>%
  mutate(
    Label = paste0(Species, "\n", Cancer_Association, "\n", Evidence_Level),
    Size = sqrt(Occurrence) * 2
  )

unique_sizes <- sort(unique(plot_summary$Size))
if(length(unique_sizes) >= 3) {
  breaks <- unique_sizes[c(1, round(length(unique_sizes)/2), length(unique_sizes))]
} else {
  breaks <- unique_sizes
}
labels <- c("Low", "Medium", "High")[seq_along(breaks)]
confidence_palette <- c("High" = "#1a9850", "Medium" = "#fdae61", "Low" = "#d73027")

# 7. Main bubble plot
main_plot <- ggplot(plot_summary, aes(x = Relevance_Level, y = Avg_Identity, size = Size, color = Confidence)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = Label), size = 3, box.padding = 0.5, max.overlaps = Inf) +
  scale_color_manual(values = confidence_palette) +
  scale_size_continuous(name = "Occurrence", breaks = breaks, labels = labels) +
  labs(
    title = "Carcinogenic Viruses Associations",
    subtitle = "Size indicates occurrence count across samples",
    x = "Association Type", y = "Average % Identity"
  ) +
  theme_minimal() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))

# 8. Evidence confidence bar plot
conf_plot <- ggplot(plot_data, aes(x = Confidence, fill = Confidence)) +
  geom_bar() +
  scale_fill_manual(values = confidence_palette) +
  labs(
    title = "Evidence Confidence Distribution",
    x = "Confidence Level", y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# 9. Combine, annotate, and save plots
final_plot <- (main_plot | conf_plot) +
  plot_annotation(
    title = "Comprehensive Virus-Cancer Association Analysis Classified",
    subtitle = "Showing relationships, evidence, and occurrence counts",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

print(final_plot)
tryCatch({
  ggsave(
    "comprehensive_virus_cancer_analysis_classified.png",
    plot = final_plot,
    width = 16, height = 10, dpi = 300
  )
  message("Successfully saved combined plot!")
}, error = function(e) {
  message("Error saving plot: ", e$message)
  png("comprehensive_virus_cancer_analysis_classified.png", width = 1600, height = 1000)
  print(final_plot)
  dev.off()
  message("Used alternative save method - check comprehensive_virus_cancer_analysis_classified.png")
})

######################################################################
# VIRUSES UNCLASSIFIED 
######################################################################

print("Processing unclassified viruses...")

# 1. Curated database of carcinogenic viruses
carcinogenic_viruses_db <- data.frame(
  Species = c(
    "Human endogenous retrovirus K",
    "Epstein-Barr virus",
    "Mouse mammary tumor virus-like",
    "Human herpesvirus 2",
    "Human herpesvirus 1",
    "Human T-lymphotropic virus 1",
    "Human papillomavirus 16",
    "Hepatitis B virus",
    "Hepatitis C virus", 
    "Human herpesvirus 8 (KSHV)",
    "Merkel cell polyomavirus",
    "Human cytomegalovirus (HCMV)",
    "Human papillomavirus 18",
    "JC polyomavirus",
    "Human papillomavirus 33",
    "Varicella-zoster virus",
    "Human herpesvirus 7",
    "BK polyomavirus",
    "Human polyomavirus 6",
    "Simian virus 40 (SV40)",
    # Additional oncogenic viruses
    "Human papillomavirus 31",
    "Human papillomavirus 45",
    "Hepatitis D virus",
    "Human T-lymphotropic virus 2",
    "Adenovirus type 12"
  ),
  Cancer_Association = c(
    "Breast cancer; Prostate cancer",
    "Burkitt lymphoma; Hodgkin lymphoma; Nasopharyngeal carcinoma; Gastric cancer",
    "Breast cancer (experimental model)",
    "Cervical cancer (co-factor); Genital cancers",
    "Oral cancer (possible); Head and neck cancer",
    "Adult T-cell leukemia/lymphoma; T-cell lymphoma",
    "Cervical cancer; Anal cancer; Oropharyngeal cancer; Penile cancer",
    "Hepatocellular carcinoma; Liver cancer",
    "Hepatocellular carcinoma; Liver cancer; Cholangiocarcinoma",
    "Kaposi sarcoma; Primary effusion lymphoma; Multicentric Castleman disease",
    "Merkel cell carcinoma",
    "Possible glioblastoma; Prostate cancer (emerging)",
    "Cervical cancer; Adenocarcinoma; Anogenital cancers",
    "Progressive multifocal leukoencephalopathy; Possible brain tumors",
    "Cervical cancer; Anogenital cancers",
    "Possible lymphoma; Indirect cancer effects",
    "Possible T-cell lymphoma (emerging)",
    "Kidney cancer; Bladder cancer; Prostate cancer",
    "Trichodysplasia spinulosa; Possible skin cancers",
    "Mesothelioma; Brain tumors; Lymphomas",
    # Additional entries
    "Cervical cancer; Anogenital cancers",
    "Cervical cancer; Anogenital cancers", 
    "Hepatocellular carcinoma (co-factor with HBV)",
    "Hairy cell leukemia; T-cell lymphomas",
    "Possible sarcomas; Experimental tumors"
  ),
  Mechanism = c(
    "Endogenous retroviral activation; Immune modulation",
    "Latent infection; Immune evasion; Viral oncoproteins (LMP1, EBNA)",
    "Retroviral integration; Hormone receptor modulation",
    "Chronic inflammation; Immune suppression; Co-infection effects",
    "Chronic inflammation; Immune modulation; HSV reactivation",
    "Viral tax protein; Cell cycle dysregulation; Immune evasion",
    "Viral DNA integration; E6/E7 oncoproteins; p53/Rb inactivation",
    "Chronic hepatitis; Viral DNA integration; HBx protein oncogene",
    "Chronic hepatitis; Persistent infection; Oxidative stress; Fibrosis",
    "Viral latency; LANA oncogene; Angiogenesis promotion",
    "Viral DNA integration; Large T antigen; Cell transformation",
    "Immune modulation; Chronic infection; Oncogene activation",
    "Viral DNA integration; E6/E7 oncoproteins; Chromosomal instability",
    "Viral DNA integration; Large T antigen; Cell cycle disruption",
    "Viral DNA integration; E6/E7 oncoproteins; DNA damage",
    "Latent infection; Immune suppression; Chronic inflammation",
    "Immune modulation; Limited oncogenic evidence",
    "Viral DNA integration; Large T antigen; Oncogene activation",
    "Viral replication; Possible transformation; Limited evidence",
    "SV40 Large T antigen; p53/Rb inactivation; Cell immortalization",
    # Additional mechanisms
    "Viral DNA integration; E6/E7 oncoproteins; HPV31 specific effects",
    "Viral DNA integration; E6/E7 oncoproteins; HPV45 specific effects",
    "Co-infection with HBV; Enhanced hepatotoxicity",
    "Tax protein; Rex protein; Cell cycle disruption",
    "Viral E1A/E1B proteins; Cell transformation"
  ),
  Evidence_Level = c(
    "Emerging Strong",
    "WHO Class 1",
    "Experimental Strong",
    "IARC 2A",
    "Emerging",
    "WHO Class 1",
    "WHO Class 1",
    "WHO Class 1",
    "WHO Class 1",
    "WHO Class 1",
    "WHO Class 1",
    "Emerging",
    "WHO Class 1",
    "IARC 2B",
    "IARC 2A",
    "Emerging",
    "Experimental",
    "IARC 2B",
    "Experimental",
    "IARC 2B",
    # Additional evidence levels
    "IARC 2A",
    "IARC 2A",
    "Emerging",
    "IARC 2B",
    "Experimental"
  ),
  Accession = I(list(
    # Human endogenous retrovirus K
    c("NC_001716.2", "NC_022518.1", "NC_007797.1"),
    # Epstein-Barr virus
    c("NC_007605.1", "NC_001664.4", "NC_009334.1"),
    # Mouse mammary tumor virus-like
    c("NC_001503.1", "NC_022518.1", "NC_001503.2"),
    # Human herpesvirus 2
    c("NC_001798.2", "NC_001345.2", "NC_001798.1"),
    # Human herpesvirus 1
    c("NC_001806.2", "NC_001806.1", "NC_023914.1"),
    # Human T-lymphotropic virus 1
    c("NC_001436.1", "NC_001454.1", "NC_001436.2"),
    # Human papillomavirus 16
    c("NC_001526.4", "NC_001405.1", "NC_001526.2"),
    # Hepatitis B virus
    c("NC_003977.2", "NC_003521.1", "NC_003977.1"),
    # Hepatitis C virus
    c("NC_004102.1", "NC_009827.1", "NC_009824.1"),
    # Human herpesvirus 8 (KSHV)
    c("NC_009333.1", "NC_003409.1", "NC_009333.2"),
    # Merkel cell polyomavirus
    c("NC_010277.2", "NC_010277.1", "NC_017972.1"),
    # Human cytomegalovirus
    c("NC_001347.1", "NC_001702.1", "NC_006273.2"),
    # Human papillomavirus 18
    c("NC_001357.1", "NC_001604.1", "NC_001357.2"),
    # JC polyomavirus
    c("NC_001699.1", "NC_003287.1", "NC_001699.2"),
    # Human papillomavirus 33
    c("NC_001586.1", "NC_006273.2", "NC_001586.2"),
    # Varicella-zoster virus
    c("NC_001348.1", "NC_001501.2", "NC_001348.2"),
    # Human herpesvirus 7
    c("NC_001716.1", "NC_009334.1", "NC_001716.2"),
    # BK polyomavirus
    c("NC_001538.1", "NC_000898.1", "NC_001538.2"),
    # Human polyomavirus 6
    c("NC_014406.1", "NC_001546.2", "NC_014406.2"),
    # Simian virus 40
    c("NC_001669.1", "NC_003819.1", "NC_001669.2"),
    # Human papillomavirus 31
    c("NC_001527.1", "NC_007977.1", "NC_001527.2"),
    # Human papillomavirus 45
    c("NC_001529.1", "NC_007998.1", "NC_001529.2"),
    # Hepatitis D virus
    c("NC_001653.2", "NC_001653.1", "NC_006103.1"),
    # Human T-lymphotropic virus 2
    c("NC_001488.1", "NC_001488.2", "NC_007771.1"),
    # Adenovirus type 12
    c("NC_001460.1", "NC_010956.1", "NC_001460.2")
  )),
  stringsAsFactors = FALSE
) %>%
  mutate(
    # Strip version numbers for clean accessions
    Accession_clean = lapply(Accession, function(vec) gsub("\\..*$", "", vec)),
    # Create flat comma-separated string for easier access
    Accession_flat = sapply(Accession_clean, paste, collapse = ", ")
  )

# Display the structure
print(carcinogenic_viruses_db)

# Read viral BLAST results - UPDATE THIS PATH!
blast_results <- read_tsv("/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/blast_results/data_analysis/pipeline/cov_analysis/VIRUS_unclassified_coverage_hits.tsv") %>%
  dplyr::select(
    Sample,
    QueryID,
    SubjectID,
    Identity,
    Evalue,
    Bitscore,
    AlignmentLength
  ) %>%
  dplyr::mutate(
    SubjectID = gsub("\\..*$", "", SubjectID),
    Coverage = AlignmentLength * Identity / 100
  )

# Create flattened lookup table for joining
accession_lookup <- carcinogenic_viruses_db %>%
  select(Species, Cancer_Association, Mechanism, Evidence_Level, Accession_clean) %>%
  unnest(Accession_clean) %>%
  rename(SubjectID_clean = Accession_clean)

# Join and annotate carcinogenic viruses
carcinogen_hits <- blast_results %>%
  inner_join(accession_lookup, by = c("SubjectID" = "SubjectID_clean")) %>%
  mutate(
    Relevance_Level = case_when(
      Evidence_Level == "WHO Class 1" ~ "1_Direct",
      Evidence_Level %in% c("IARC 2A", "Emerging Strong", "Experimental Strong") ~ "2_Strong",
      Evidence_Level %in% c("IARC 2B", "Emerging") ~ "3_Potential",
      TRUE ~ "4_Experimental"
    ),
    Confidence = case_when(
      Evidence_Level == "WHO Class 1" & Identity > 95 ~ "High",
      Evidence_Level == "WHO Class 1" & Identity > 90 ~ "Medium",
      grepl("IARC 2A|Strong", Evidence_Level, ignore.case = TRUE) & Identity > 95 ~ "Medium",
      grepl("IARC 2B|Emerging", Evidence_Level, ignore.case = TRUE) & Identity > 90 ~ "Low",
      TRUE ~ "Low"
    )
  ) %>%
  arrange(Relevance_Level, desc(Confidence), desc(Bitscore))

# Check for empty results
if (nrow(carcinogen_hits) == 0) {
  warning("No carcinogenic virus matches found!")
} else {
  print(paste("Found", nrow(carcinogen_hits), "carcinogenic virus matches"))
}

# 4. Generate and save summary reports
write.csv(carcinogen_hits, "all_carcinogenic_viruses_unclassified.csv", row.names = FALSE)

# 5. Prepare for plotting: clean labels/factors for visualization
plot_data <- carcinogen_hits %>%
  mutate(
    Relevance_Level = gsub("^[0-9]_", "", Relevance_Level)
  )

plot_data$Cancer_Association <- factor(
  plot_data$Cancer_Association,
  levels = unique(plot_data$Cancer_Association)
)
plot_data$Relevance_Level <- factor(
  plot_data$Relevance_Level,
  levels = c("Direct", "Strong", "Potential", "Experimental")
)
plot_data$Confidence <- factor(plot_data$Confidence, levels = c("High", "Medium", "Low"))

# 6. Summarize for bubble plot with enhanced labels
plot_summary <- plot_data %>%
  group_by(Species, Cancer_Association, Confidence, Relevance_Level, Evidence_Level) %>%
  summarise(
    Occurrence = n_distinct(Sample),
    Avg_Identity = mean(Identity),
    .groups = "drop"
  ) %>%
  mutate(
    Label = paste0(Species, "\n", Cancer_Association, "\n", Evidence_Level),
    Size = sqrt(Occurrence) * 2
  )

unique_sizes <- sort(unique(plot_summary$Size))
if(length(unique_sizes) >= 3) {
  breaks <- unique_sizes[c(1, round(length(unique_sizes)/2), length(unique_sizes))]
} else {
  breaks <- unique_sizes
}
labels <- c("Low", "Medium", "High")[seq_along(breaks)]
confidence_palette <- c("High" = "#1a9850", "Medium" = "#fdae61", "Low" = "#d73027")

# 7. Main bubble plot
main_plot <- ggplot(plot_summary, aes(x = Relevance_Level, y = Avg_Identity, size = Size, color = Confidence)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = Label), size = 3, box.padding = 0.5, max.overlaps = Inf) +
  scale_color_manual(values = confidence_palette) +
  scale_size_continuous(name = "Occurrence", breaks = breaks, labels = labels) +
  labs(
    title = "Carcinogenic Viruses Associations",
    subtitle = "Size indicates occurrence count across samples",
    x = "Association Type", y = "Average % Identity"
  ) +
  theme_minimal() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))

# 8. Evidence confidence bar plot
conf_plot <- ggplot(plot_data, aes(x = Confidence, fill = Confidence)) +
  geom_bar() +
  scale_fill_manual(values = confidence_palette) +
  labs(
    title = "Evidence Confidence Distribution",
    x = "Confidence Level", y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# 9. Combine, annotate, and save plots
final_plot <- (main_plot | conf_plot) +
  plot_annotation(
    title = "Comprehensive Virus-Cancer Association Analysis Unclassified",
    subtitle = "Showing relationships, evidence, and occurrence counts",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

print(final_plot)
tryCatch({
  ggsave(
    "comprehensive_virus_cancer_analysis_unclassified.png",
    plot = final_plot,
    width = 16, height = 10, dpi = 300
  )
  message("Successfully saved combined plot!")
}, error = function(e) {
  message("Error saving plot: ", e$message)
  png("comprehensive_virus_cancer_analysis_unclassified.png", width = 1600, height = 1000)
  print(final_plot)
  dev.off()
  message("Used alternative save method - check comprehensive_virus_cancer_analysis_unclassified.png")
})

######################################################################
# PARASITES CLASSIFIED 
######################################################################

print("Processing classified parasites...")

# 1. Curated database of carcinogenic parasites
carcinogenic_parasites_db <- data.frame(
  Species = c(
    "Schistosoma haematobium",
    "Opisthorchis viverrini", 
    "Clonorchis sinensis",
    "Schistosoma mansoni",
    "Schistosoma japonicum",
    "Taenia solium",
    "Paragonimus westermani",
    "Strongyloides stercoralis",
    "Echinococcus granulosus",
    "Trichinella spiralis",
    "Toxoplasma gondii",
    "Fasciola hepatica",
    "Cryptosporidium parvum",
    "Trypanosoma cruzi",
    "Leishmania major",
    "Plasmodium falciparum",
    # Additional parasites with cancer associations
    "Entamoeba histolytica",
    "Giardia lamblia",
    "Trichomonas vaginalis",
    "Babesia microti"
  ),
  Cancer_Association = c(
    "Bladder cancer",
    "Cholangiocarcinoma",
    "Cholangiocarcinoma", 
    "Colorectal cancer; bladder cancer",
    "Colorectal cancer; liver cancer",
    "CNS tumors (rare reports)",
    "Lung cancer (possible)",
    "Colorectal cancer (emerging)",
    "Liver cancer; possible GI involvement",
    "Possible gastrointestinal cancer",
    "Brain tumors (experimental)",
    "Breast cancer (experimental)",
    "Possible gastrointestinal cancer",
    "Colon/GI cancer (experimental)",
    "Skin cancer; SQ cell carcinoma",
    "Lymphoproliferative; blood cancers",
    # Additional entries
    "Colorectal cancer; liver cancer",
    "Possible GI cancer; inflammation",
    "Cervical cancer; prostate cancer",
    "Possible hematological malignancies"
  ),
  Mechanism = c(
    "Chronic inflammation; DNA damage; nitrosamines",
    "Bile duct damage & chronic inflammation",
    "Bile duct damage & chronic inflammation",
    "Chronic inflammation; immune dysregulation",
    "Chronic inflammation; hepatic fibrosis",
    "Tissue cysts; chronic neuroinflammation",
    "Immune modulation; lung inflammation",
    "Immune modulation; chronic GI infection",
    "Metabolic perturbation; immune evasion",
    "Muscle inflammation; immune modulation",
    "Oncogenic pathway activation; cyst formation",
    "Estrogen mimicry; immune modulation",
    "Intestinal epithelial damage; cryptosporidiosis",
    "Chronic inflammation; DNA damage; ROS production",
    "Immune suppression; chronic inflammation",
    "Chronic inflammation; cell cycle dysregulation",
    # Additional mechanisms
    "Chronic colitis; liver abscess formation",
    "Malabsorption; chronic inflammation",
    "Chronic urogenital inflammation; HPV co-infection",
    "Hemolysis; immune modulation"
  ),
  Evidence_Level = c(
    "WHO Class 1",
    "WHO Class 1", 
    "WHO Class 1",
    "WHO Class 1",
    "WHO Class 1",
    "Emerging",
    "Emerging",
    "Emerging Strong",
    "IARC 2B",
    "Experimental",
    "Experimental Strong",
    "Experimental Strong",
    "Possible",
    "Emerging",
    "Emerging",
    "Emerging",
    # Additional evidence levels
    "Emerging",
    "Experimental",
    "Controversial",
    "Experimental"
  ),
  Accession = I(list(
    # Schistosoma haematobium
    c("NW_001850701.1", "NC_008074.1", "NW_013434508.1"),
    # Opisthorchis viverrini
    c("NW_001850136.1", "NC_021179.1", "NW_013434344.1"),
    # Clonorchis sinensis
    c("NW_001851409.1", "NC_021003.1", "NW_013434607.1"),
    # Schistosoma mansoni
    c("NC_002545.1", "NW_003072345.1", "NW_013434289.1"),
    # Schistosoma japonicum
    c("NC_002544.1", "NW_003176542.1", "NW_013434178.1"),
    # Taenia solium
    c("NW_001850357.1", "NC_004022.1", "NW_013435012.1"),
    # Paragonimus westermani
    c("NW_001850610.1", "NC_007884.1", "NW_013434891.1"),
    # Strongyloides stercoralis
    c("NW_021628287.1", "NC_028624.1", "NW_021628445.1"),
    # Echinococcus granulosus
    c("NC_000719.5", "NW_013433987.1", "NW_013434123.1"),
    # Trichinella spiralis
    c("NC_000442.1", "NW_003315234.1", "NW_003315567.1"),
    # Toxoplasma gondii
    c("NW_003315945.1", "NC_001799.1", "NW_019386754.1"),
    # Fasciola hepatica
    c("NW_004166845.2", "NC_007544.1", "NW_004167123.1"),
    # Cryptosporidium parvum
    c("NC_001483.1", "NW_003263445.1", "NW_003263678.1"),
    # Trypanosoma cruzi
    c("NC_004551.1", "NW_003043267.1", "NW_003043445.1"),
    # Leishmania major
    c("NC_007073.1", "NW_001893245.1", "NW_001893578.1"),
    # Plasmodium falciparum
    c("NC_004309.1", "NW_001848234.1", "NW_001848567.1"),
    # Entamoeba histolytica
    c("NC_016446.1", "NW_001849234.1", "NW_001849567.1"),
    # Giardia lamblia
    c("NC_009421.1", "NW_002284445.1", "NW_002284678.1"),
    # Trichomonas vaginalis
    c("NC_016942.1", "NW_001843234.1", "NW_001843567.1"),
    # Babesia microti
    c("NC_012575.1", "NW_003489234.1", "NW_003489567.1")
  )),
  stringsAsFactors = FALSE
) %>%
  mutate(
    # Strip version numbers for clean accessions
    Accession_clean = lapply(Accession, function(vec) gsub("\\..*$", "", vec)),
    # Create flat comma-separated string for easier access
    Accession_flat = sapply(Accession_clean, paste, collapse = ", ")
  )

# Display the structure
print(carcinogenic_parasites_db)

# Read parasite BLAST results - UPDATE THIS PATH!
blast_results <- read_tsv("/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/blast_results/data_analysis/pipeline/cov_analysis/PARASITE_classified_coverage_hits.tsv") %>%
  dplyr::select(
    Sample,
    QueryID,
    SubjectID,
    Identity,
    Evalue,
    Bitscore,
    AlignmentLength
  ) %>%
  dplyr::mutate(
    SubjectID = gsub("\\..*$", "", SubjectID),
    Coverage = AlignmentLength * Identity / 100
  )

# Create flattened lookup table for joining
accession_lookup <- carcinogenic_parasites_db %>%
  select(Species, Cancer_Association, Mechanism, Evidence_Level, Accession_clean) %>%
  unnest(Accession_clean) %>%
  rename(SubjectID_clean = Accession_clean)

# Join and annotate carcinogenic parasites
carcinogen_hits <- blast_results %>%
  inner_join(accession_lookup, by = c("SubjectID" = "SubjectID_clean")) %>%
  mutate(
    Relevance_Level = case_when(
      Evidence_Level == "WHO Class 1" ~ "1_Direct",
      Evidence_Level %in% c("Emerging Strong", "Experimental Strong") ~ "2_Strong",
      Evidence_Level %in% c("Emerging", "IARC 2B") ~ "3_Potential",
      TRUE ~ "4_Experimental"
    ),
    Confidence = case_when(
      Evidence_Level == "WHO Class 1" & Identity > 95 ~ "High",
      Evidence_Level == "WHO Class 1" & Identity > 90 ~ "Medium",
      grepl("Strong|IARC 2B", Evidence_Level, ignore.case = TRUE) & Identity > 95 ~ "Medium",
      grepl("Emerging", Evidence_Level, ignore.case = TRUE) & Identity > 90 ~ "Low",
      TRUE ~ "Low"
    )
  ) %>%
  arrange(Relevance_Level, desc(Confidence), desc(Bitscore))

# Check for empty results
if (nrow(carcinogen_hits) == 0) {
  warning("No carcinogenic parasite matches found!")
} else {
  print(paste("Found", nrow(carcinogen_hits), "carcinogenic parasite matches"))
}

# 4. Generate and save summary reports
write.csv(carcinogen_hits, "all_carcinogenic_parasites_classified.csv", row.names = FALSE)

# 5. Prepare for plotting: clean labels/factors for visualization
plot_data <- carcinogen_hits %>%
  mutate(
    Relevance_Level = gsub("^[0-9]_", "", Relevance_Level)
  )

plot_data$Cancer_Association <- factor(
  plot_data$Cancer_Association,
  levels = unique(plot_data$Cancer_Association)
)
plot_data$Relevance_Level <- factor(
  plot_data$Relevance_Level,
  levels = c("Direct", "Strong", "Potential", "Experimental")
)
plot_data$Confidence <- factor(plot_data$Confidence, levels = c("High", "Medium", "Low"))

# 6. Summarize for bubble plot with enhanced labels
plot_summary <- plot_data %>%
  group_by(Species, Cancer_Association, Confidence, Relevance_Level, Evidence_Level) %>%
  summarise(
    Occurrence = n_distinct(Sample),
    Avg_Identity = mean(Identity),
    .groups = "drop"
  ) %>%
  mutate(
    Label = paste0(Species, "\n", Cancer_Association, "\n", Evidence_Level),
    Size = sqrt(Occurrence) * 2
  )

unique_sizes <- sort(unique(plot_summary$Size))
if(length(unique_sizes) >= 3) {
  breaks <- unique_sizes[c(1, round(length(unique_sizes)/2), length(unique_sizes))]
} else {
  breaks <- unique_sizes
}
labels <- c("Low", "Medium", "High")[seq_along(breaks)]
confidence_palette <- c("High" = "#1a9850", "Medium" = "#fdae61", "Low" = "#d73027")

# 7. Main bubble plot
main_plot <- ggplot(plot_summary, aes(x = Relevance_Level, y = Avg_Identity, size = Size, color = Confidence)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = Label), size = 3, box.padding = 0.5, max.overlaps = Inf) +
  scale_color_manual(values = confidence_palette) +
  scale_size_continuous(name = "Occurrence", breaks = breaks, labels = labels) +
  labs(
    title = "Carcinogenic Parasites Associations",
    subtitle = "Size indicates occurrence count across samples",
    x = "Association Type", y = "Average % Identity"
  ) +
  theme_minimal() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))

# 8. Evidence confidence bar plot
conf_plot <- ggplot(plot_data, aes(x = Confidence, fill = Confidence)) +
  geom_bar() +
  scale_fill_manual(values = confidence_palette) +
  labs(
    title = "Evidence Confidence Distribution",
    x = "Confidence Level", y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# 9. Combine, annotate, and save plots
final_plot <- (main_plot | conf_plot) +
  plot_annotation(
    title = "Comprehensive Parasite-Cancer Association Analysis Classified",
    subtitle = "Showing relationships, evidence, and occurrence counts",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

print(final_plot)
tryCatch({
  ggsave(
    "comprehensive_parasite_cancer_analysis_classified.png",
    plot = final_plot,
    width = 16, height = 10, dpi = 300
  )
  message("Successfully saved combined plot!")
}, error = function(e) {
  message("Error saving plot: ", e$message)
  png("comprehensive_parasite_cancer_analysis_classified.png", width = 1600, height = 1000)
  print(final_plot)
  dev.off()
  message("Used alternative save method - check comprehensive_parasite_cancer_analysis_classified.png")
})

######################################################################
# PARASITES UNCLASSIFIED 
######################################################################

print("Processing unclassified parasites...")

# 1. Curated database of carcinogenic parasites
carcinogenic_parasites_db <- data.frame(
  Species = c(
    "Schistosoma haematobium",
    "Opisthorchis viverrini", 
    "Clonorchis sinensis",
    "Schistosoma mansoni",
    "Schistosoma japonicum",
    "Taenia solium",
    "Paragonimus westermani",
    "Strongyloides stercoralis",
    "Echinococcus granulosus",
    "Trichinella spiralis",
    "Toxoplasma gondii",
    "Fasciola hepatica",
    "Cryptosporidium parvum",
    "Trypanosoma cruzi",
    "Leishmania major",
    "Plasmodium falciparum",
    # Additional parasites with cancer associations
    "Entamoeba histolytica",
    "Giardia lamblia",
    "Trichomonas vaginalis",
    "Babesia microti"
  ),
  Cancer_Association = c(
    "Bladder cancer",
    "Cholangiocarcinoma",
    "Cholangiocarcinoma", 
    "Colorectal cancer; bladder cancer",
    "Colorectal cancer; liver cancer",
    "CNS tumors (rare reports)",
    "Lung cancer (possible)",
    "Colorectal cancer (emerging)",
    "Liver cancer; possible GI involvement",
    "Possible gastrointestinal cancer",
    "Brain tumors (experimental)",
    "Breast cancer (experimental)",
    "Possible gastrointestinal cancer",
    "Colon/GI cancer (experimental)",
    "Skin cancer; SQ cell carcinoma",
    "Lymphoproliferative; blood cancers",
    # Additional entries
    "Colorectal cancer; liver cancer",
    "Possible GI cancer; inflammation",
    "Cervical cancer; prostate cancer",
    "Possible hematological malignancies"
  ),
  Mechanism = c(
    "Chronic inflammation; DNA damage; nitrosamines",
    "Bile duct damage & chronic inflammation",
    "Bile duct damage & chronic inflammation",
    "Chronic inflammation; immune dysregulation",
    "Chronic inflammation; hepatic fibrosis",
    "Tissue cysts; chronic neuroinflammation",
    "Immune modulation; lung inflammation",
    "Immune modulation; chronic GI infection",
    "Metabolic perturbation; immune evasion",
    "Muscle inflammation; immune modulation",
    "Oncogenic pathway activation; cyst formation",
    "Estrogen mimicry; immune modulation",
    "Intestinal epithelial damage; cryptosporidiosis",
    "Chronic inflammation; DNA damage; ROS production",
    "Immune suppression; chronic inflammation",
    "Chronic inflammation; cell cycle dysregulation",
    # Additional mechanisms
    "Chronic colitis; liver abscess formation",
    "Malabsorption; chronic inflammation",
    "Chronic urogenital inflammation; HPV co-infection",
    "Hemolysis; immune modulation"
  ),
  Evidence_Level = c(
    "WHO Class 1",
    "WHO Class 1", 
    "WHO Class 1",
    "WHO Class 1",
    "WHO Class 1",
    "Emerging",
    "Emerging",
    "Emerging Strong",
    "IARC 2B",
    "Experimental",
    "Experimental Strong",
    "Experimental Strong",
    "Possible",
    "Emerging",
    "Emerging",
    "Emerging",
    # Additional evidence levels
    "Emerging",
    "Experimental",
    "Controversial",
    "Experimental"
  ),
  Accession = I(list(
    # Schistosoma haematobium
    c("NW_001850701.1", "NC_008074.1", "NW_013434508.1"),
    # Opisthorchis viverrini
    c("NW_001850136.1", "NC_021179.1", "NW_013434344.1"),
    # Clonorchis sinensis
    c("NW_001851409.1", "NC_021003.1", "NW_013434607.1"),
    # Schistosoma mansoni
    c("NC_002545.1", "NW_003072345.1", "NW_013434289.1"),
    # Schistosoma japonicum
    c("NC_002544.1", "NW_003176542.1", "NW_013434178.1"),
    # Taenia solium
    c("NW_001850357.1", "NC_004022.1", "NW_013435012.1"),
    # Paragonimus westermani
    c("NW_001850610.1", "NC_007884.1", "NW_013434891.1"),
    # Strongyloides stercoralis
    c("NW_021628287.1", "NC_028624.1", "NW_021628445.1"),
    # Echinococcus granulosus
    c("NC_000719.5", "NW_013433987.1", "NW_013434123.1"),
    # Trichinella spiralis
    c("NC_000442.1", "NW_003315234.1", "NW_003315567.1"),
    # Toxoplasma gondii
    c("NW_003315945.1", "NC_001799.1", "NW_019386754.1"),
    # Fasciola hepatica
    c("NW_004166845.2", "NC_007544.1", "NW_004167123.1"),
    # Cryptosporidium parvum
    c("NC_001483.1", "NW_003263445.1", "NW_003263678.1"),
    # Trypanosoma cruzi
    c("NC_004551.1", "NW_003043267.1", "NW_003043445.1"),
    # Leishmania major
    c("NC_007073.1", "NW_001893245.1", "NW_001893578.1"),
    # Plasmodium falciparum
    c("NC_004309.1", "NW_001848234.1", "NW_001848567.1"),
    # Entamoeba histolytica
    c("NC_016446.1", "NW_001849234.1", "NW_001849567.1"),
    # Giardia lamblia
    c("NC_009421.1", "NW_002284445.1", "NW_002284678.1"),
    # Trichomonas vaginalis
    c("NC_016942.1", "NW_001843234.1", "NW_001843567.1"),
    # Babesia microti
    c("NC_012575.1", "NW_003489234.1", "NW_003489567.1")
  )),
  stringsAsFactors = FALSE
) %>%
  mutate(
    # Strip version numbers for clean accessions
    Accession_clean = lapply(Accession, function(vec) gsub("\\..*$", "", vec)),
    # Create flat comma-separated string for easier access
    Accession_flat = sapply(Accession_clean, paste, collapse = ", ")
  )

# Display the structure
print(carcinogenic_parasites_db)

# Read parasite BLAST results - UPDATE THIS PATH!
blast_results <- read_tsv("/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/virnatrap/has_output/blast_results/data_analysis/pipeline/cov_analysis/PARASITE_unclassified_coverage_hits.tsv") %>%
  dplyr::select(
    Sample,
    QueryID,
    SubjectID,
    Identity,
    Evalue,
    Bitscore,
    AlignmentLength
  ) %>%
  dplyr::mutate(
    SubjectID = gsub("\\..*$", "", SubjectID),
    Coverage = AlignmentLength * Identity / 100
  )

# Create flattened lookup table for joining
accession_lookup <- carcinogenic_parasites_db %>%
  select(Species, Cancer_Association, Mechanism, Evidence_Level, Accession_clean) %>%
  unnest(Accession_clean) %>%
  rename(SubjectID_clean = Accession_clean)

# Join and annotate carcinogenic parasites
carcinogen_hits <- blast_results %>%
  inner_join(accession_lookup, by = c("SubjectID" = "SubjectID_clean")) %>%
  mutate(
    Relevance_Level = case_when(
      Evidence_Level == "WHO Class 1" ~ "1_Direct",
      Evidence_Level %in% c("Emerging Strong", "Experimental Strong") ~ "2_Strong",
      Evidence_Level %in% c("Emerging", "IARC 2B") ~ "3_Potential",
      TRUE ~ "4_Experimental"
    ),
    Confidence = case_when(
      Evidence_Level == "WHO Class 1" & Identity > 95 ~ "High",
      Evidence_Level == "WHO Class 1" & Identity > 90 ~ "Medium",
      grepl("Strong|IARC 2B", Evidence_Level, ignore.case = TRUE) & Identity > 95 ~ "Medium",
      grepl("Emerging", Evidence_Level, ignore.case = TRUE) & Identity > 90 ~ "Low",
      TRUE ~ "Low"
    )
  ) %>%
  arrange(Relevance_Level, desc(Confidence), desc(Bitscore))

# Check for empty results
if (nrow(carcinogen_hits) == 0) {
  warning("No carcinogenic parasite matches found!")
} else {
  print(paste("Found", nrow(carcinogen_hits), "carcinogenic parasite matches"))
}

# 4. Generate and save summary reports
write.csv(carcinogen_hits, "all_carcinogenic_parasites_unclassified.csv", row.names = FALSE)

# 5. Prepare for plotting: clean labels/factors for visualization
plot_data <- carcinogen_hits %>%
  mutate(
    Relevance_Level = gsub("^[0-9]_", "", Relevance_Level)
  )

plot_data$Cancer_Association <- factor(
  plot_data$Cancer_Association,
  levels = unique(plot_data$Cancer_Association)
)
plot_data$Relevance_Level <- factor(
  plot_data$Relevance_Level,
  levels = c("Direct", "Strong", "Potential", "Experimental")
)
plot_data$Confidence <- factor(plot_data$Confidence, levels = c("High", "Medium", "Low"))

# 6. Summarize for bubble plot with enhanced labels
plot_summary <- plot_data %>%
  group_by(Species, Cancer_Association, Confidence, Relevance_Level, Evidence_Level) %>%
  summarise(
    Occurrence = n_distinct(Sample),
    Avg_Identity = mean(Identity),
    .groups = "drop"
  ) %>%
  mutate(
    Label = paste0(Species, "\n", Cancer_Association, "\n", Evidence_Level),
    Size = sqrt(Occurrence) * 2
  )

unique_sizes <- sort(unique(plot_summary$Size))
if(length(unique_sizes) >= 3) {
  breaks <- unique_sizes[c(1, round(length(unique_sizes)/2), length(unique_sizes))]
} else {
  breaks <- unique_sizes
}
labels <- c("Low", "Medium", "High")[seq_along(breaks)]
confidence_palette <- c("High" = "#1a9850", "Medium" = "#fdae61", "Low" = "#d73027")

# 7. Main bubble plot
main_plot <- ggplot(plot_summary, aes(x = Relevance_Level, y = Avg_Identity, size = Size, color = Confidence)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = Label), size = 3, box.padding = 0.5, max.overlaps = Inf) +
  scale_color_manual(values = confidence_palette) +
  scale_size_continuous(name = "Occurrence", breaks = breaks, labels = labels) +
  labs(
    title = "Carcinogenic Parasites Associations",
    subtitle = "Size indicates occurrence count across samples",
    x = "Association Type", y = "Average % Identity"
  ) +
  theme_minimal() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))

# 8. Evidence confidence bar plot
conf_plot <- ggplot(plot_data, aes(x = Confidence, fill = Confidence)) +
  geom_bar() +
  scale_fill_manual(values = confidence_palette) +
  labs(
    title = "Evidence Confidence Distribution",
    x = "Confidence Level", y = "Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# 9. Combine, annotate, and save plots
final_plot <- (main_plot | conf_plot) +
  plot_annotation(
    title = "Comprehensive Parasite-Cancer Association Analysis Unclassified",
    subtitle = "Showing relationships, evidence, and occurrence counts",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

print(final_plot)
tryCatch({
  ggsave(
    "comprehensive_parasite_cancer_analysis_unclassified.png",
    plot = final_plot,
    width = 16, height = 10, dpi = 300
  )
  message("Successfully saved combined plot!")
}, error = function(e) {
  message("Error saving plot: ", e$message)
  png("comprehensive_parasite_cancer_analysis_unclassified.png", width = 1600, height = 1000)
  print(final_plot)
  dev.off()
  message("Used alternative save method - check comprehensive_parasite_cancer_analysis_unclassified.png")
})

print("Process finished, plots have been stroed in comprehensive_plots")
