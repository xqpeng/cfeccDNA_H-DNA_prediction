# =============================================================================
# Triplex Toolkit: FASTA Sequence Filtering based on LS_5bp ID Lists
# =============================================================================

# Check and load required libraries
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Bioconductor package 'Biostrings' is missing. Install via: BiocManager::install('Biostrings')")
}
library(Biostrings)
library(tools)

# =============================================================================
# Command Line Argument Parsing
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript filet_5bp.R <Input_FASTA_Dir> <ID_List_Dir> <Output_Dir>")
}

# 1. Input FASTA Directory (Output from Step 1)
INPUT_FASTA_DIR <- args[1]

# 2. Filtering ID List Directory (LS_5bp)
ID_LIST_DIR <- args[2]

# 3. Output Directory (For filtered FASTA files)
OUTPUT_DIR <- args[3]

# Create output directory if it does not exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# =============================================================================
# Helper Functions (Robust Cleaning)
# =============================================================================

#' Clean and Standardize IDs (Aggressive Mode)
#' 
#' This function normalizes IDs to ensure strictly accurate matching by removing
#' suffixes, metadata, and special characters.
#' 
#' @param id_string The raw ID string.
#' @return A cleaned ID string.
clean_id <- function(id_string) {
  # 1. Remove content after the pipe '|' (Specific to FASTA headers with metadata)
  id_string <- sub("\\|.*$", "", id_string)
  # 2. Remove content after the first whitespace (Specific to ID lists or FASTQ comments)
  id_string <- sub("\\s+.*$", "", id_string)
  # 3. Remove paired-end indicators like /1 or /2
  id_string <- sub("/[12]$", "", id_string)
  # 4. Remove leading format markers like > or @
  id_string <- sub("^[@>]", "", id_string)
  # 5. Trim leading and trailing whitespace
  id_string <- trimws(id_string)
  return(id_string)
}

# =============================================================================
# Main Processing Logic
# =============================================================================

process_filtering <- function() {
  # Retrieve all ID list files starting with "LS_5bp_"
  id_files <- list.files(ID_LIST_DIR, pattern = "^LS_5bp_", full.names = TRUE)
  
  cat(sprintf("Detected %d ID list files. Starting processing...\n", length(id_files)))
  cat("----------------------------------------------------------------\n")
  
  count_success <- 0
  
  for (id_path in id_files) {
    id_filename <- basename(id_path)
    
    # 1. Parse Sample Name
    # Assumes filename format: LS_5bp_SampleName (e.g., LS_5bp_AB010 -> AB010)
    # Removes potential file extensions like .txt
    sample_name <- file_path_sans_ext(sub("^LS_5bp_", "", id_filename))
    
    # 2. Find Corresponding FASTA File
    # Assumes FASTA filename format: TS_SampleName.fasta
    # Uses regex anchor '$' to ensure exact suffix matching
    fasta_pattern <- paste0("TS_", sample_name, ".fasta$")
    fasta_files <- list.files(INPUT_FASTA_DIR, pattern = fasta_pattern, full.names = TRUE)
    
    if (length(fasta_files) == 0) {
      cat(sprintf("[SKIP] FASTA file not found for sample %s (Expected TS_%s.fasta)\n", sample_name, sample_name))
      next
    }
    
    # Use the first match found
    fasta_path <- fasta_files[1] 
    
    # 3. Read and Clean ID List
    # Use tryCatch to prevent crashing on empty/locked files
    raw_lines <- tryCatch({
      readLines(id_path, warn = FALSE)
    }, error = function(e) NULL)
    
    # Filter out empty lines
    raw_lines <- raw_lines[raw_lines != ""]
    
    if (length(raw_lines) == 0) {
      cat(sprintf("[SKIP] ID file is empty or unreadable: %s\n", id_filename))
      next
    }
    
    # Clean IDs and convert to a unique set for O(1) lookup performance
    target_ids <- clean_id(raw_lines)
    target_ids_set <- unique(target_ids)
    
    # 4. Read FASTA File
    dna_set <- tryCatch({
      readDNAStringSet(fasta_path)
    }, error = function(e) {
      cat(sprintf("[ERROR] Failed to read FASTA: %s\n", basename(fasta_path)))
      return(NULL)
    })
    
    if (is.null(dna_set)) next
    
    # 5. Execute Filtering
    # Clean FASTA sequence names to match the format of the ID list
    fasta_names_cleaned <- clean_id(names(dna_set))
    
    # Identify indices of sequences that exist in the target ID set
    keep_indices <- which(fasta_names_cleaned %in% target_ids_set)
    
    if (length(keep_indices) == 0) {
      cat(sprintf("[WARNING] Sample %s: No sequences retained (No ID matches).\n", sample_name))
      # Diagnostic print to help user debug if 0 matches occur
      cat(sprintf("   > Example ID from List:  '%s'\n", head(target_ids_set, 1)))
      cat(sprintf("   > Example ID from FASTA: '%s'\n", head(fasta_names_cleaned, 1)))
      next
    }
    
    filtered_dna_set <- dna_set[keep_indices]
    
    # 6. Save Filtered Output
    output_filename <- paste0(sample_name, "_filtered.fasta")
    output_path <- file.path(OUTPUT_DIR, output_filename)
    
    writeXStringSet(filtered_dna_set, output_path)
    
    cat(sprintf("[DONE] %s -> Retained %d / %d sequences -> %s\n", 
                sample_name, length(filtered_dna_set), length(dna_set), output_filename))
    
    count_success <- count_success + 1
  }
  
  cat("----------------------------------------------------------------\n")
  cat(sprintf("Processing Complete! Successfully processed %d samples.\nResults saved to: %s\n", count_success, OUTPUT_DIR))
}

# Execute Main Function
process_filtering()