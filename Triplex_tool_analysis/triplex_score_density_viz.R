# =============================================================================
# Triplex Score Distribution Analysis (Customized Density Visualization)
# =============================================================================

# 1. Package Initialization
if (!require("ggridges")) install.packages("ggridges")
library(tidyverse)
library(readr)
library(ggplot2)
library(ggridges)

# --- 1. Path Configuration (Modified for Command Line Args) ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Error: Please provide <Results_Directory> and <Group_File_Path>")
}

POS_DIR <- args[1]
GROUP_FILE <- args[2]

# Output directory for publication-quality score analysis
OUTPUT_DIR <- file.path(POS_DIR, "Combined_Score_Analysis_Final_Best")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# --- 2. Color Palette for Biological Groups ---
my_colors <- c(
  "Control"  = "#525252", 
  "HCC"      = "#D93636", 
  "HNSCC"    = "#3478C4", 
  "NPC"      = "#54A356", 
  "CRC"      = "#9655AA", 
  "LC"       = "#D6A629"
)

# --- 3. Global Font Settings (Arial) ---
if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
  font_family <- "Arial"
} else {
  font_family <- "sans"
}

# =============================================================================
# Module 1: Data Acquisition and Preprocessing
# =============================================================================
cat(">>> [1/3] Preparing datasets...\n")

# 1. Parse group metadata and filter for specific library types
readme_lines <- readLines(GROUP_FILE)
header_idx <- which(grepl("Sample", readme_lines) & grepl("Group", readme_lines))[1]
group_data <- read.table(text = readme_lines[header_idx:length(readme_lines)], 
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                         quote = "", fill = TRUE, check.names = FALSE)

group_data <- group_data %>% 
  filter(library == "NonBS-seq") %>%
  select(Sample, Group) %>%
  mutate(Sample = trimws(gsub("\\.BS$", "", Sample)), Group = trimws(Group))

# 2. Function to aggregate triplex score data from multiple sample directories
read_score_data <- function(dir_path) {
  files <- list.files(dir_path, pattern = "best_triplex_sequences_detailed.csv", recursive = TRUE, full.names = TRUE)
  if(length(files) == 0) return(NULL)
  
  data_list <- lapply(files, function(x) {
    df <- tryCatch(read.csv(x, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) {
      fname <- basename(dirname(x))
      sname <- gsub("analysis_", "", fname)
      sname <- gsub("^TS_", "", sname)
      df$Sample <- sname
      return(df[, c("Sample", "Score")]) 
    }
    return(NULL)
  })
  do.call(rbind, data_list[!sapply(data_list, is.null)])
}

pos_data <- read_score_data(POS_DIR)
if(is.null(pos_data)) stop("Error: No data found in the specified directory!")

# 3. Merge score metrics with group identifiers
all_data <- pos_data %>% inner_join(group_data, by = "Sample")

# 4. Define factor levels for consistent legend and plot ordering
all_data$Group <- factor(all_data$Group, levels = names(my_colors))

cat(">>> Data ready for processing:", nrow(all_data), "rows found.\n")

# =============================================================================
# Module 2: Visualization Generation (Density Plot)
# =============================================================================
cat(">>> [2/3] Generating Density Plot...\n")

# Determine X-axis upper bounds using 99.5th percentile for robust scaling
max_val <- quantile(all_data$Score, 0.995, na.rm=TRUE) * 1.1
max_score <- max(15, max_val) 

# Optimized theme for scientific publication
theme_plot1 <- theme_bw() + theme(
  text = element_text(family = font_family),
  plot.title = element_blank(), 
  
  # Axis titles and text configuration: plain face and specific font sizes
  axis.title = element_text(size = 24, family = font_family, face = "plain"),
  axis.text = element_text(size = 22, color = "black", family = font_family),
  
  # Legend positioning: inset top-right with custom alignment
  legend.position = c(0.99, 0.99), 
  legend.justification = c(1, 1),
  
  # Legend styling: white background with black solid border
  legend.background = element_rect(fill = "white", color = "black", size = 0.5),
  legend.key = element_rect(fill = "transparent", color = NA),
  legend.text = element_text(size = 20, family = font_family),
  legend.title = element_blank(),
  
  panel.grid.minor = element_blank(),
  
  plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt")
)

# Force bold/solid line style in the legend for better visibility
legend_style_line <- guides(
  color = guide_legend(override.aes = list(size = 1.5, linetype = 1, alpha = 1))
)

p1 <- ggplot(all_data, aes(x = Score, color = Group)) +
  # Density profile with outline-only display (alpha=0)
  geom_density(size = 1, alpha = 0, key_glyph = "path") + 
  
  # X-axis configuration: starting from 10 with 5-unit increments
  scale_x_continuous(
    limits = c(10, max_score), 
    expand = c(0, 0),
    breaks = seq(10, max_score, 5)
  ) +
  scale_color_manual(values = my_colors) +
  
  theme_plot1 +
  legend_style_line +
  
  labs(
    title = NULL,
    x = "Triplex score",
    y = "Density"
  )

# Export high-resolution file (300 DPI)
ggsave(file.path(OUTPUT_DIR, "01_Score_Density_Final_Custom.png"), p1, 
       width = 23.44, height = 18.36, units = "cm", dpi = 300)

cat("\n=================================================\n")
cat(">>>> Task Completed successfully.\n")
cat(">>>> Output path: ", OUTPUT_DIR, "\n")
cat(">>>> File generated: 01_Score_Density_Final_Custom.png\n")
cat("=================================================\n")