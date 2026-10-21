# =============================================================================
# Triplex Analysis Tool: LS_5bp Filtering and Robust ID Matching
# =============================================================================

library(dplyr)
library(readr)
library(stringr)

# --- 1. Path Configuration ---
RESULTS_ROOT_DIR <- "D:\\triplex\\triplex(R)\\batch_results_Fast_Parallel"
ID_LIST_DIR <- "D:\\triplex\\Target_sequences\\LS_5bp"
STATS_OUTPUT_FILE <- file.path(RESULTS_ROOT_DIR, "LS_5bp_Statistics_Summary_Robust.csv")

# =============================================================================
# Helper Function: Aggressive ID Normalization
# =============================================================================
clean_ids_aggressively <- function(ids) {
  # 1. Remove metadata suffixes following the pipe character
  ids <- sub("\\|.*$", "", ids)
  # 2. Strip descriptions following white space
  ids <- sub("\\s+.*$", "", ids)
  # 3. Remove paired-end read indicators (/1 or /2)
  ids <- sub("/[12]$", "", ids)
  # 4. Strip leading FASTQ/FASTA format markers (@ or >)
  ids <- sub("^[@>]", "", ids)
  # 5. Trim leading and trailing whitespace
  ids <- trimws(ids)
  return(ids)
}

# =============================================================================
# Core Processing Logic: Sample-wise ID Extraction and Statistics
# =============================================================================
process_sample_ls_5bp <- function(sample_name) {
  
  # --- A. Resource Discovery ---
  id_file <- file.path(ID_LIST_DIR, paste0("LS_5bp_", sample_name))
  if (!file.exists(id_file) && file.exists(paste0(id_file, ".txt"))) {
    id_file <- paste0(id_file, ".txt")
  }
  
  all_result_dirs <- list.dirs(RESULTS_ROOT_DIR, recursive = FALSE)
  target_dir <- all_result_dirs[grep(paste0("analysis_.*", sample_name, "$"), basename(all_result_dirs))]
  
  if (!file.exists(id_file)) return(list(status="skip", msg="ID file missing"))
  if (length(target_dir) == 0) return(list(status="skip", msg="Result directory missing"))
  
  target_dir <- target_dir[1]
  
  # --- B. Target ID Loading and Normalization ---
  target_ids_raw <- tryCatch({
    lines <- readLines(id_file, warn = FALSE)
    lines[lines != ""] 
  }, error = function(e) NULL)
  
  if (length(target_ids_raw) == 0) return(list(status="skip", msg="ID file empty"))
  
  target_ids_clean <- clean_ids_aggressively(target_ids_raw)
  target_ids_unique <- unique(target_ids_clean) 
  
  # --- C. Results CSV Loading and Normalization ---
  summary_csv <- file.path(target_dir, "analysis_summary.csv")
  if (!file.exists(summary_csv)) return(list(status="skip", msg="Summary CSV missing"))
  
  df_summary <- read.csv(summary_csv, stringsAsFactors = FALSE)
  
  # Create a normalized ID column for cross-referencing while preserving raw IDs
  df_summary$CleanID <- clean_ids_aggressively(df_summary$SequenceID)
  
  # --- D. Intersection Matching ---
  matches <- df_summary$CleanID %in% target_ids_unique
  count_found_rows <- sum(matches)
  
  # Diagnostic feedback for failed matches
  if (count_found_rows == 0) {
    cat(sprintf("\n[Diagnosis] Sample %s match failed!\n", sample_name))
    cat("  > Cleaned ID examples from LS_5bp file:\n")
    print(head(target_ids_clean, 3))
    cat("  > Cleaned ID examples from Summary CSV:\n")
    print(head(df_summary$CleanID, 3))
    cat("  --------------------------------------\n")
    return(list(status="skip", msg="ID mismatch: No overlaps found"))
  }
  
  filtered_summary <- df_summary[matches, ]
  filtered_summary$CleanID <- NULL 
  
  # --- E. Detailed Sequence Extraction (Conditional) ---
  detail_csv <- file.path(target_dir, "all_triplex_sequences_detailed.csv")
  if (file.exists(detail_csv)) {
    df_detail <- read.csv(detail_csv, stringsAsFactors = FALSE)
    df_detail$CleanID <- clean_ids_aggressively(df_detail$SequenceID)
    filtered_detail <- df_detail[df_detail$CleanID %in% target_ids_unique, ]
    filtered_detail$CleanID <- NULL
    
    write.csv(filtered_detail, file.path(target_dir, "filtered_LS_5bp_detailed.csv"), row.names = FALSE)
  }
  
  write.csv(filtered_summary, file.path(target_dir, "filtered_LS_5bp_summary.csv"), row.names = FALSE)
  
  # --- F. Metric Calculation ---
  count_total <- length(unique(target_ids_clean))
  
  # Quantify triplex findings with 'Successful' analysis status
  count_triplex_found <- sum(filtered_summary$AnalysisStatus == "³É¹¦")
  
  rate <- round((count_triplex_found / count_total) * 100, 2)
  
  return(list(
    status = "success",
    msg = paste("Extracted", count_found_rows, "records"),
    data = data.frame(
      Sample = sample_name,
      Total_LS_5bp_Seqs = count_total,
      Found_Triplexes = count_triplex_found,
      Detection_Rate = rate,
      stringsAsFactors = FALSE
    )
  ))
}

# =============================================================================
# Batch Execution and Statistical Summary
# =============================================================================

id_files <- list.files(ID_LIST_DIR, pattern = "LS_5bp_")
cat("Found", length(id_files), "ID files. Initiating matching pipeline...\n")

global_stats <- data.frame()
count_ok <- 0
pb <- txtProgressBar(min = 0, max = length(id_files), style = 3)

for (i in seq_along(id_files)) {
  f <- id_files[i]
  sample_name <- gsub("^LS_5bp_", "", tools::file_path_sans_ext(f))
  
  res <- process_sample_ls_5bp(sample_name)
  
  if (res$status == "success") {
    global_stats <- rbind(global_stats, res$data)
    count_ok <- count_ok + 1
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# Export final statistical summary
if (nrow(global_stats) > 0) {
  global_stats <- global_stats %>% arrange(desc(Detection_Rate))
  write.csv(global_stats, STATS_OUTPUT_FILE, row.names = FALSE)
  
  cat("\n\n============================================\n")
  cat("Processing complete. Summaries saved for", count_ok, "samples.\n")
  cat("Master summary file:\n", STATS_OUTPUT_FILE, "\n")
  cat("============================================\n")
  print(head(global_stats))
} else {
  cat("\n\nNo samples matched successfully. Check diagnostic logs for ID format discrepancies.\n")
}