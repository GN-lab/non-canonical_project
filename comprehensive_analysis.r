#!/usr/bin/env Rscript

.libPaths(c("/data/rds/DMP/UCEC/EVOLIMMU/csalas_rds/gaurav_rds/R_libs", .libPaths()))

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(data.table)
    library(RColorBrewer)
    library(patchwork)
    library(ggrepel)
})

# Set up graphics device for headless environment
options(bitmapType = 'cairo')
Sys.setenv(DISPLAY = '')

# Create plots directory
plots_dir <- "comprehensive_plots"
dir.create(plots_dir, showWarnings = FALSE)

# ========================================
# LOAD AND PREPARE DATA
# ========================================

# Load coverage data
coverage_file <- "cov_analysis/coverage_depth_stats.tsv"
if (file.exists(coverage_file)) {
    coverage_data <- fread(coverage_file)
    cat("Columns in coverage data:", names(coverage_data), "\n")
    cat("Microbial groups found:", unique(coverage_data$MicrobialGroup), "\n")
} else {
    stop(paste("Coverage file not found:", coverage_file))
}

# Include ALL microbial groups
coverage_data <- coverage_data[MicrobialGroup %in% c("BACT", "VIRUS", "PARASITE", "FUNGI")]

# Add sample categories
coverage_data[, SampleCategory := fcase(
  grepl("^HT01", Sample), "Healthy",
  grepl("^HT02", Sample), "Familial", 
  grepl("^HT03|^HT04|^HT05", Sample), "BRCA_mutated",
  grepl("^HT06|^HT07|^HT08|^HT09|^HT10", Sample), "Tumor",
  grepl("^ES", Sample), "Surrounding",
  grepl("^CPCT", Sample), "Metastatic",
  default = "Other"
)]

# Color palettes
microbial_colors <- c(
  BACT = "#1f77b4",      # blue
  VIRUS = "#ff7f0e",     # orange  
  PARASITE = "#2ca02c",  # green
  FUNGI = "#d62728"      # red
)

category_colors <- brewer.pal(length(unique(coverage_data$SampleCategory)), "Set1")
names(category_colors) <- unique(coverage_data$SampleCategory)

# ========================================
# 1. CORE OVERVIEW PLOTS
# ========================================

# 1.1 Microbial Group Distribution (Overview)
p1 <- coverage_data[, .(Count = .N), by = MicrobialGroup] %>%
  ggplot(aes(x = reorder(MicrobialGroup, -Count), y = Count, fill = MicrobialGroup)) +
  geom_col() +
  geom_text(aes(label = scales::comma(Count)), vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = microbial_colors) +
  labs(title = "Total Sequences by Microbial Group",
       x = "Microbial Group", y = "Number of Sequences") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(plots_dir, "1_microbial_group_distribution.png"), p1, width = 8, height = 6, dpi = 300)

# 1.2 Sample Category Distribution
p2 <- coverage_data[, .(Count = .N), by = SampleCategory] %>%
  ggplot(aes(x = reorder(SampleCategory, -Count), y = Count, fill = SampleCategory)) +
  geom_col() +
  geom_text(aes(label = scales::comma(Count)), vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = category_colors) +
  labs(title = "Sequences by Sample Category", 
       x = "Sample Category", y = "Count") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(plots_dir, "2_sample_category_distribution.png"), p2, width = 8, height = 6, dpi = 300)

# ========================================
# 2. SAMPLE-MICROBE INTERACTION PLOTS
# ========================================

# 2.1 Sample vs Microbial Group Heatmap (ABSOLUTE COUNTS)
heatmap_data <- coverage_data[, .(Count = .N), by = .(Sample, MicrobialGroup, SampleCategory)]

p3 <- ggplot(heatmap_data, aes(x = MicrobialGroup, y = Sample, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), color = "black", size = 2.5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "steelblue", 
                      name = "Sequence\nCount") +
  labs(title = "Sample-Microbe Interaction Matrix",
       subtitle = "Absolute sequence counts per sample and microbial group",
       x = "Microbial Group", y = "Sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(plots_dir, "3_sample_microbe_heatmap.png"), p3, width = 10, height = 12, dpi = 300)

# 2.2 Stacked Bar by Sample Category
p4 <- coverage_data[, .(Count = .N), by = .(SampleCategory, MicrobialGroup)] %>%
  ggplot(aes(x = SampleCategory, y = Count, fill = MicrobialGroup)) +
  geom_col(position = "stack") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), 
            color = "white", fontface = "bold", size = 3) +
  scale_fill_manual(values = microbial_colors) +
  labs(title = "Sequence Distribution by Sample Category",
       x = "Sample Category", y = "Count", fill = "Microbial Group") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(plots_dir, "4_category_stacked_bar.png"), p4, width = 10, height = 8, dpi = 300)

# ========================================
# 3. ALIGNMENT ANALYSIS PLOTS
# ========================================

# 3.1 Alignment Depth Distribution
p5 <- ggplot(coverage_data, aes(x = NumAlignments, fill = MicrobialGroup)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = microbial_colors) +
  labs(title = "Alignment Depth Distribution by Microbial Group",
       x = "Number of Alignments", y = "Density", fill = "Microbial Group") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  scale_x_log10()

ggsave(file.path(plots_dir, "5_alignment_depth_distribution.png"), p5, width = 10, height = 6, dpi = 300)

# 3.2 Average Alignments per Microbial Group
p6 <- coverage_data[, .(AvgAlignments = mean(NumAlignments)), by = MicrobialGroup] %>%
  ggplot(aes(x = reorder(MicrobialGroup, -AvgAlignments), y = AvgAlignments, fill = MicrobialGroup)) +
  geom_col() +
  geom_text(aes(label = round(AvgAlignments, 1)), vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = microbial_colors) +
  labs(title = "Average Alignment Depth by Microbial Group",
       x = "Microbial Group", y = "Average Number of Alignments") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(plots_dir, "6_avg_alignments_by_group.png"), p6, width = 8, height = 6, dpi = 300)

# ========================================
# 4. COMPOSITION ANALYSIS PLOTS
# ========================================

# 4.1 Microbial Composition by Sample Category (Donut charts)
composition_data <- coverage_data[, .(Count = .N), by = .(SampleCategory, MicrobialGroup)]
composition_data[, Percent := Count/sum(Count)*100, by = SampleCategory]
composition_data[, ypos := cumsum(Percent) - 0.5 * Percent, by = SampleCategory]

p7 <- ggplot(composition_data, aes(x = 2, y = Percent, fill = MicrobialGroup)) +
  geom_col(color = "white", width = 1) +
  geom_text(aes(y = ypos, label = ifelse(Percent > 5, paste0(round(Percent, 1), "%"), "")), 
            color = "black", size = 3, fontface = "bold") +
  coord_polar(theta = "y", start = 0) +
  scale_fill_manual(values = microbial_colors) +
  facet_wrap(~SampleCategory, ncol = 3) +
  labs(title = "Microbial Composition by Sample Category",
       fill = "Microbial Group") +
  theme_void() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold")) +
  xlim(0.5, 2.5)

ggsave(file.path(plots_dir, "7_composition_donuts.png"), p7, width = 12, height = 10, dpi = 300)

# 4.2 Relative Abundance Heatmap
abundance_data <- coverage_data[, .(Count = .N), by = .(SampleCategory, MicrobialGroup)]
abundance_data[, RelativeAbundance := Count/sum(Count), by = SampleCategory]

p8 <- ggplot(abundance_data, aes(x = MicrobialGroup, y = SampleCategory, fill = RelativeAbundance)) +
  geom_tile(color = "white") +
  geom_text(aes(label = scales::percent(RelativeAbundance, accuracy = 0.1)), 
            color = "black", size = 4, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "steelblue", 
                      labels = scales::percent, name = "Relative\nAbundance") +
  labs(title = "Relative Microbial Abundance by Sample Category",
       x = "Microbial Group", y = "Sample Category") +
  theme_classic() +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(plots_dir, "8_relative_abundance_heatmap.png"), p8, width = 8, height = 6, dpi = 300)

# ========================================
# 5. TOP SAMPLES ANALYSIS
# ========================================

# 5.1 Top Samples by Microbial Group
top_samples <- coverage_data[, .(Count = .N), by = .(Sample, MicrobialGroup)] %>%
  .[order(-Count)] %>%
  head(20)

p9 <- ggplot(top_samples, aes(x = reorder(Sample, Count), y = Count, fill = MicrobialGroup)) +
  geom_col() +
  geom_text(aes(label = Count), hjust = -0.2, size = 3, fontface = "bold") +
  scale_fill_manual(values = microbial_colors) +
  labs(title = "Top 20 Samples by Sequence Count",
       x = "Sample", y = "Number of Sequences") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 14, face = "bold")) +
  coord_flip()

ggsave(file.path(plots_dir, "9_top_samples.png"), p9, width = 10, height = 8, dpi = 300)

# ========================================
# 6. INTERACTIVE-READY PLOTS (for reports)
# ========================================

# 6.1 Comprehensive overview plot
overview_plot <- (p1 + p2) / (p6 + p8) +
  plot_annotation(title = "Comprehensive Microbial Analysis Overview",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

ggsave(file.path(plots_dir, "10_comprehensive_overview.png"), overview_plot, width = 16, height = 12, dpi = 300)

# ========================================
# SUMMARY STATISTICS AND EXPORT
# ========================================

# Generate summary statistics
summary_stats <- coverage_data[, .(
  Total_Sequences = .N,
  Unique_Sequences = uniqueN(qseqid),
  Total_Alignments = sum(NumAlignments),
  Avg_Alignments_Per_Sequence = round(mean(NumAlignments), 2),
  Max_Alignments = max(NumAlignments)
), by = .(MicrobialGroup, SampleCategory)]

# Write summary files
fwrite(summary_stats, file.path(plots_dir, "summary_statistics.csv"))
fwrite(coverage_data, file.path(plots_dir, "full_coverage_data.csv"))

# Print summary
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Plots saved in directory:", plots_dir, "\n")
cat("\nTop 3 samples per microbial group:\n")

for(group in unique(coverage_data$MicrobialGroup)) {
  top <- coverage_data[MicrobialGroup == group, .(Count = .N), by = Sample][order(-Count)][1:3]
  cat("\n", group, ":\n")
  print(top)
}



cat("\nKey visualizations created:\n")
cat("1. Microbial group distribution\n")
cat("2. Sample category distribution\n") 
cat("3. Sample-microbe interaction heatmap\n")
cat("4. Category-wise stacked bars\n")
cat("5. Alignment depth analysis\n")
cat("6. Composition donut charts\n")
cat("7. Relative abundance heatmap\n")
cat("8. Top samples analysis\n")
cat("9. Comprehensive overview\n")

cat("\nNext steps: Consider functional annotation analysis and pathway enrichment!\n")