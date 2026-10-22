# =============================================================================
# Triplex Distribution Analysis (High-Performance Vectorized Optimization)
# Developed using data.table and stringi for efficient genomic data processing
# =============================================================================

# 1. Dependency Management and Initialization
if (!require("data.table")) install.packages("data.table")
if (!require("stringi")) install.packages("stringi")

library(data.table)
library(stringi)
library(tidyverse)
library(ggplot2)
library(Biostrings) 
library(parallel) 

# --- 1. Path Configuration (Modified for Command Line Args) ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Error: Please provide <Results_Directory> <Refined_FASTA_Dir> <Group_File_Path>")
}

POS_RES_DIR   <- args[1]
POS_FASTA_DIR <- args[2]
GROUP_FILE    <- args[3]

# Version 6 fast-processing output directory
OUTPUT_DIR <- file.path(POS_RES_DIR, "Combined_Distance_Analysis_Select5Prime_V6_Fast")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# Cache directory for intermediate rds objects
TEMP_DIR <- file.path(OUTPUT_DIR, "temp_chunks_fast")
if (!dir.exists(TEMP_DIR)) dir.create(TEMP_DIR)

# --- 2. Aesthetic and Visualization Definitions ---
my_colors <- c(
  "Control" = "#525252",
  "HCC"     = "#D93636",
  "HNSCC"   = "#3478C4",
  "NPC"     = "#54A356",
  "CRC"     = "#9655AA",
  "LC"      = "#D6A629"
)

if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
  font_family <- "Arial"
} else {
  font_family <- "sans"
}

# --- 3. Plot Theme (Optimized for Large-Scale Data Visualization) ---
custom_theme <- theme_bw() + theme(
  text = element_text(family = font_family),
  plot.title = element_blank(), 
  
  # Axis settings: 24pt labels, plain face
  axis.title = element_text(size = 24, family = font_family, face = "plain"), 
  
  # Axis tick text: 22pt
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22, color = "black", family = font_family),
  axis.text.y = element_text(size = 22, color = "black", family = font_family),
  
  # Legend positioning: Top-right inset with customized border
  legend.position = c(0.99, 0.99), 
  legend.justification = c(1, 1),
  legend.background = element_rect(fill = "white", color = "black", size = 0.5),
  legend.key = element_rect(fill = "transparent", color = NA),
  legend.title = element_blank(),
  
  # Legend font settings
  legend.text = element_text(size = 20, family = font_family),
  
  panel.grid.major = element_line(size = 0.2), 
  panel.grid.minor = element_blank(),
  
  plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt")
)

legend_fix <- guides(
  fill = "none", 
  color = guide_legend(override.aes = list(size = 1.5, linetype = 1, alpha = 1))
)

# =============================================================================
# Core Worker Function (Fully Vectorized Mapping Pipeline)
# =============================================================================
process_and_filter_folder <- function(folder_name, res_dir, fasta_dir) {
  library(data.table)
  library(stringi)
  library(Biostrings)
  
  # Limit internal data.table threading for parallel efficiency
  setDTthreads(1) 
  
  sample_name_core <- gsub("^analysis_", "", folder_name)
  
  fasta_filename <- paste0(sample_name_core, ".fasta")
  fasta_path <- file.path(fasta_dir, fasta_filename)
  csv_path <- file.path(res_dir, folder_name, "all_triplex_sequences_detailed.csv")
  out_csv_path <- file.path(res_dir, folder_name, "closest_5prime_triplex.csv")
  
  if (!file.exists(fasta_path) || !file.exists(csv_path)) return(NULL)
  
  # 1. High-speed I/O for sequence metadata
  dt <- tryCatch(fread(csv_path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  
  # 2. Genomic sequence data ingestion
  fasta_seqs <- tryCatch(readDNAStringSet(fasta_path), error = function(e) NULL)
  if (is.null(fasta_seqs)) return(NULL)
  
  ref_seqs_str <- as.character(fasta_seqs)
  ref_names <- sapply(strsplit(names(ref_seqs_str), "\\|"), `[`, 1)
  ref_dt <- data.table(Short_ID = ref_names, Original_Seq = ref_seqs_str)
  
  # 3. Vectorized join between query data and reference sequences
  dt[, Short_ID := tstrsplit(SequenceID, "\\|", keep = 1)]
  dt <- merge(dt, ref_dt, by = "Short_ID", all.x = FALSE) 
  
  if (nrow(dt) == 0) return(NULL)
  
  # 4. Rapid substring pattern matching using C++ based stringi
  dt[, Start_Pos := stringi::stri_locate_first_fixed(Original_Seq, Sequence)[, 1]]
  
  # 5. Selective Reverse Complement handling for unmapped reads
  na_indices <- which(is.na(dt$Start_Pos))
  
  if (length(na_indices) > 0) {
    seqs_to_rc <- DNAStringSet(dt$Sequence[na_indices])
    rc_seqs <- as.character(reverseComplement(seqs_to_rc))
    
    original_seqs_subset <- dt$Original_Seq[na_indices]
    found_pos <- stringi::stri_locate_first_fixed(original_seqs_subset, rc_seqs)[, 1]
    
    dt$Start_Pos[na_indices] <- found_pos
  }
  
  dt <- dt[!is.na(Start_Pos)]
  if (nrow(dt) == 0) return(NULL)
  
  # 6. Optimized physical distance calculations (0-based)
  dt[, `:=`(
    Len_Seq = nchar(Original_Seq),
    Len_Tri = nchar(Sequence)
  )]
  
  dt[, `:=`(
    Dist_To_5 = Start_Pos - 1,
    End_Pos = Start_Pos + Len_Tri - 1
  )]
  
  dt[, Dist_To_3 := Len_Seq - End_Pos]
  
  # Compute minimum distance to nearest genomic end
  dt[, Calc_MinDist := pmin(Dist_To_5, Dist_To_3)]
  
  # 7. Extract the triplex sequence closest to the 5' end per sequence ID
  filtered_dt <- dt[order(Dist_To_5), .SD[1L], by = Short_ID]
  
  # 8. Fast-writing of processed results
  fwrite(filtered_dt, out_csv_path)
  
  clean_sample <- gsub("^TS_", "", sample_name_core)
  return(data.frame(
    Sample = clean_sample, 
    MinDist = filtered_dt$Calc_MinDist, 
    stringsAsFactors = FALSE
  ))
}

# =============================================================================
# Multi-threaded Scheduler (18-Core Configuration)
# =============================================================================
run_batch_processing <- function(res_dir, fasta_dir, n_cores, batch_size) {
  res_folders <- list.dirs(res_dir, full.names = FALSE, recursive = FALSE)
  res_folders <- res_folders[grep("^analysis_", res_folders)]
  
  total_folders <- length(res_folders)
  if(total_folders == 0) return(NULL)
  
  chunks <- split(res_folders, ceiling(seq_along(res_folders) / batch_size))
  cat(sprintf("\n>>> Initializing: %d folders in %d batches\n", total_folders, length(chunks)))
  
  all_results <- list()
  
  for (i in seq_along(chunks)) {
    chunk_folders <- chunks[[i]]
    cat(sprintf("   Batch %d/%d (18 Threads)... ", i, length(chunks)))
    
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, {
      library(data.table)
      library(stringi)
      library(Biostrings)
    })
    
    batch_res_list <- parLapply(cl, chunk_folders, process_and_filter_folder, 
                                res_dir = res_dir, fasta_dir = fasta_dir)
    stopCluster(cl)
    
    batch_df <- do.call(rbind, batch_res_list)
    
    if (!is.null(batch_df) && nrow(batch_df) > 0) {
      temp_file <- file.path(TEMP_DIR, paste0("batch_", i, ".rds"))
      saveRDS(batch_df, temp_file)
      cat(sprintf("Done! Saved rows: %d\n", nrow(batch_df)))
      all_results[[i]] <- batch_df
    } else {
      cat("No valid data.\n")
    }
    gc() 
  }
  return(do.call(rbind, all_results))
}

# =============================================================================
# Main Pipeline Execution
# =============================================================================

n_cores <- 18 
BATCH_SIZE <- 20 

cat("\n>>> [1/3] Loading Group Metadata...\n")
readme_lines <- readLines(GROUP_FILE)
header_idx <- which(grepl("Sample", readme_lines) & grepl("Group", readme_lines))[1]
group_data <- read.table(text = readme_lines[header_idx:length(readme_lines)], header = TRUE, sep = "\t", quote = "", fill = TRUE)

group_data <- group_data %>% 
  filter(library == "NonBS-seq") %>% 
  select(Sample, Group) %>% 
  mutate(Sample = trimws(gsub("\\.BS$", "", Sample)), Group = trimws(Group)) %>%
  filter(Group %in% names(my_colors)) 

cat("\n>>> [2/3] Executing High-Performance Calculation (5'-end Proximity)...\n")
start_t <- Sys.time()
pos_result <- run_batch_processing(POS_RES_DIR, POS_FASTA_DIR, n_cores, batch_size = BATCH_SIZE)
print(difftime(Sys.time(), start_t))

if (!is.null(pos_result)) {
  all_data <- pos_result %>% inner_join(group_data, by = "Sample")
  all_data$Group <- factor(all_data$Group, levels = c(names(my_colors)))
  
  cat("\n>>> [3/3] Generating Visualization...\n")
  
  max_x <- 80 
  
  # Density distribution plot for triplex distance metrics
  p <- ggplot(all_data, aes(x = MinDist, color = Group)) +
    geom_density(size = 1, alpha = 0, key_glyph = "path") + 
    
    scale_x_continuous(limits = c(0, max_x), expand = c(0, 0), breaks = seq(0, 80, 10)) +
    scale_color_manual(values = my_colors) +
    
    custom_theme + 
    legend_fix + 
    
    labs(
      title = NULL,
      x = "Triplex distance to nearest 5' end (bp)",
      y = "Density"
    )
  
  out_file <- file.path(OUTPUT_DIR, "04_Distance_Distribution_Select5Prime_Fast.png")
  ggsave(out_file, p, width = 23.44, height = 18.36, units = "cm", dpi = 300)
  
  cat("\n Processing Complete!\n")
  cat(" Optimized for Speed: C++ string matching and vectorized math.\n")
  cat(" Primary Output: closest_5prime_triplex.csv\n")
  cat(" Figure Generated: ", out_file, "\n")
  
} else {
  cat(" Error: No data generated.\n")
}