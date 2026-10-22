# =============================================================================
# Triplex Length Distribution Analysis (Publication-Ready Customized Version)
# =============================================================================

library(tidyverse)
library(readr)
library(ggplot2)

# --- 1. Path Configuration (Modified for Command Line Args) ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Error: Please provide <Results_Directory> and <Group_File_Path>")
}

POS_DIR <- args[1]
GROUP_FILE <- args[2]

# Output directory for aggregated length analysis
OUTPUT_DIR <- file.path(POS_DIR, "Combined_Length_Analysis_Best")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# --- 2. Color Scheme for Biological Groups ---
my_colors <- c(
  "Control"  = "#525252", 
  "HCC"      = "#D93636", 
  "HNSCC"    = "#3478C4", 
  "NPC"      = "#54A356", 
  "CRC"      = "#9655AA", 
  "LC"       = "#D6A629"
)

# --- 3. Global Font Configuration (Arial) ---
if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
  font_family <- "Arial"
} else {
  font_family <- "sans"
}

# --- 4. Data Loading and Preparation ---
cat(">>> [1/3] Reading group metadata...\n")
readme_lines <- readLines(GROUP_FILE)
header_idx <- which(grepl("Sample", readme_lines) & grepl("Group", readme_lines))[1]
group_data <- read.table(text = readme_lines[header_idx:length(readme_lines)], 
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                         quote = "", fill = TRUE, check.names = FALSE)

group_data <- group_data %>% 
  filter(library == "NonBS-seq") %>%
  select(Sample, Group) %>%
  mutate(Sample = trimws(gsub("\\.BS$", "", Sample)), Group = trimws(Group))

# Function to recursively aggregate length data from batch CSV results
read_length_data <- function(dir_path) {
  files <- list.files(dir_path, pattern = "best_triplex_sequences_detailed.csv", recursive = TRUE, full.names = TRUE)
  if(length(files) == 0) return(NULL)
  
  data_list <- lapply(files, function(x) {
    df <- tryCatch(read.csv(x, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) {
      fname <- basename(dirname(x))
      sname <- gsub("analysis_", "", fname)
      sname <- gsub("^TS_", "", sname)
      
      df$Sample <- sname
      if("Length" %in% colnames(df)) {
        return(df[, c("Sample", "Length")]) 
      }
    }
    return(NULL)
  })
  do.call(rbind, data_list[!sapply(data_list, is.null)])
}

cat(">>> [2/3] Extracting triplex length data from positive samples...\n")
pos_data <- read_length_data(POS_DIR)

if(is.null(pos_data)) stop("No positive data found!")

# Merging group labels with sequence metrics
all_data <- pos_data %>% inner_join(group_data, by = "Sample")
all_data$Group <- factor(all_data$Group, levels = c(names(my_colors)))

# --- 5. Statistical Visualization (Density-Length Plot) ---
cat(">>> [3/3] Generating density distribution plot...\n")

# Calculate X-axis upper limit based on the 99.5th percentile
max_len <- quantile(all_data$Length, 0.995, na.rm = TRUE) * 1.1

# Theme optimization for high-quality publication output
custom_theme <- theme_bw() + theme(
  text = element_text(family = font_family),
  plot.title = element_blank(), 
  
  # Axis titles and labels: large font, plain face
  axis.title = element_text(size = 24, family = font_family, face = "plain"),
  axis.text = element_text(size = 22, color = "black", family = font_family),
  
  # Legend placement: top-right inset with customized border and font size
  legend.position = c(0.99, 0.99), 
  legend.justification = c(1, 1),
  legend.background = element_rect(fill = "white", color = "black", size = 0.5),
  legend.key = element_rect(fill = "transparent", color = NA), 
  legend.title = element_blank(),
  legend.text = element_text(size = 20, family = font_family),
  
  panel.grid.minor = element_blank(),
  
  plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt")
)

# Overriding default legend aesthetics for better visibility
legend_style <- guides(
  color = guide_legend(
    override.aes = list(
      size = 1.5,       # Increased line thickness
      linetype = 1,     # Solid line
      alpha = 1         # Full opacity
    )
  )
)

p <- ggplot(all_data, aes(x = Length, color = Group)) +
  
  geom_density(size = 1, alpha = 0, key_glyph = "path") + 
  
  # X-axis scale configuration: starting from 10bp with 20bp intervals
  scale_x_continuous(
    limits = c(10, max_len),            
    expand = c(0, 0),
    breaks = seq(10, max_len, 20)       
  ) +
  
  scale_color_manual(values = my_colors) +
  
  custom_theme +
  legend_style + 
  
  labs(
    title = NULL,
    x = "Triplex length (bp)",
    y = "Density" 
  )

# Saving high-resolution output (300 DPI)
out_file <- file.path(OUTPUT_DIR, "Length_Distribution_Custom_BoxedLegend.png")
ggsave(out_file, p, width = 23.44, height = 18.36, units = "cm", dpi = 300)

cat("\n=================================================\n")
cat("✅ Task Completed! Visualization saved.\n")
cat("✅ Legend: Boxed with size 20 font.\n")
cat("✅ Style: Plain face axis titles.\n")
cat("Output: ", out_file, "\n")
cat("=================================================\n")