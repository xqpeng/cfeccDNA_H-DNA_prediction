# =============================================================================
# Triplex Visualization System 
# =============================================================================

library(tidyverse)
library(readr)
library(ggplot2)
library(ggpubr)

# --- 1. Path Configuration (Modified for Command Line Args) ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Error: Please provide <Results_Directory> and <Group_File_Path>")
}

RESULTS_DIR <- args[1]
GROUP_FILE  <- args[2]

OUTPUT_PLOT_DIR <- file.path(RESULTS_DIR, "Custom_Style_Plots_Positive_Final_v8_LargeFont")
if (!dir.exists(OUTPUT_PLOT_DIR)) dir.create(OUTPUT_PLOT_DIR)

# --- 2. Color Palette for Biological Groups ---
my_colors <- c(
  "Control" = "#525252",
  "HCC"      = "#D93636",
  "HNSCC"   = "#3478C4",
  "NPC"      = "#54A356",
  "CRC"      = "#9655AA",
  "LC"       = "#D6A629"
)

# Global font family initialization
if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
  font_family <- "Arial"
} else {
  font_family <- "sans"
}

# --- 3. Plot Theme Definition with Scaled Font Sizes ---
custom_theme <- theme_bw() + theme(
  text = element_text(family = font_family),
  plot.title = element_blank(), 
  
  axis.title = element_text(size = 24, family = font_family),
  
  axis.text.x = element_text(angle = 45, hjust = 1, size = 20, color = "black", family = font_family),
  axis.text.y = element_text(size = 22, color = "black", family = font_family),
  
  legend.position = "none", 
  
  panel.grid.major = element_line(size = 0.2), 
  panel.grid.minor = element_blank(),
  
  plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt")
)

# =============================================================================
# Module 1: Data Acquisition and Preprocessing
# =============================================================================

# 1. Load sample metadata and apply group filters
cat(">>> [1/5] Reading group metadata...\n")
readme_lines <- readLines(GROUP_FILE)
header_idx <- which(grepl("Sample", readme_lines) & grepl("Group", readme_lines))[1]
group_data <- read.table(text = readme_lines[header_idx:length(readme_lines)], header = TRUE, sep = "\t", quote = "")

group_data <- group_data %>% 
  filter(library == "NonBS-seq") %>%
  select(Sample, Group) %>%
  mutate(Sample = trimws(gsub("\\.BS$", "", Sample)), Group = trimws(Group)) %>%
  filter(Group %in% names(my_colors))

target_samples <- unique(group_data$Sample)

# 2. Extract total sequence counts for normalization denominator
cat(">>> [2/5] Extracting Total Reads from summary files...\n")
summary_files <- list.files(RESULTS_DIR, pattern = "analysis_summary.csv", recursive = TRUE, full.names = TRUE)

summary_list <- lapply(summary_files, function(x) {
  fname <- basename(dirname(x))
  sname <- gsub("analysis_", "", fname); sname <- gsub("^TS_", "", sname)
  if (sname %in% target_samples) {
    df <- tryCatch(read.csv(x), error = function(e) NULL)
    if (!is.null(df)) return(data.frame(Sample = sname, Total_Reads = nrow(df))) 
  }
  return(NULL)
})
df_total <- do.call(rbind, summary_list)

# 3. Retrieve identified triplex records for normalization numerator
cat(">>> [3/5] Extracting Triplex Count and Length from detailed files...\n")
detail_files <- list.files(RESULTS_DIR, pattern = "best_triplex_sequences_detailed.csv", recursive = TRUE, full.names = TRUE)

detail_list <- lapply(detail_files, function(x) {
  fname <- basename(dirname(x))
  sname <- gsub("analysis_", "", fname); sname <- gsub("^TS_", "", sname)
  if (sname %in% target_samples) {
    df <- tryCatch(read.csv(x), error = function(e) NULL)
    if (!is.null(df)) {
      df$Sample <- sname
      return(df) 
    }
  }
  return(NULL)
})
global_details <- do.call(rbind, detail_list)

# 4. Data merging and detection rate calculation
cat(">>> [4/5] Integrating data and calculating detection rates...\n")

if (is.null(global_details)) {
  df_count <- data.frame(Sample = target_samples, Found_Count = 0)
} else {
  df_count <- global_details %>% group_by(Sample) %>% summarise(Found_Count = n())
}

merged_data <- group_data %>%
  left_join(df_total, by = "Sample") %>%
  left_join(df_count, by = "Sample") %>%
  mutate(
    Total_Reads = ifelse(is.na(Total_Reads), 1, Total_Reads), 
    Found_Count = ifelse(is.na(Found_Count), 0, Found_Count),
    Rate = (Found_Count / Total_Reads) * 100
  )

# Harmonize group ordering for visualization
merged_data$Group <- factor(merged_data$Group, levels = names(my_colors))
if (!is.null(global_details)) {
  global_details <- global_details %>% inner_join(group_data, by = "Sample")
  global_details$Group <- factor(global_details$Group, levels = names(my_colors))
}

# =============================================================================
# Module 2: Statistical Visualization Generation
# =============================================================================
cat(">>> [5/5] Generating plots...\n")

# --- Detection Rate Comparison (Boxplot + Jitter + T-test) ---
p1 <- ggplot(merged_data, aes(x = Group, y = Rate, fill = Group, color = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6, key_glyph = "path") + 
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") + 
  labs(y = "Triplex detection rate (%)", x = NULL) + 
  custom_theme + 
  scale_fill_manual(values = my_colors) + 
  scale_color_manual(values = my_colors) +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "Control", vjust = 0.5, color="black", hide.ns = TRUE)

ggsave(file.path(OUTPUT_PLOT_DIR, "01_Detection_Rate_NoLegend.png"), p1, width = 23.44, height = 18.36, units = "cm", dpi = 300)

cat("\n=================================================\n")
cat("✅ Task Completed! Statistical plots generated.\n")
cat("✅ Font scale: 24pt (Titles), 20-22pt (Labels).\n")
cat("Output: ", OUTPUT_PLOT_DIR, "\n")
cat("=================================================\n")