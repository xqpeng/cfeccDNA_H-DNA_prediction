# =============================================================================
# Triplex Analysis Pipeline: Parallel Search with Resume Capability
# Implements high-throughput triplex DNA identification and structural parsing
# =============================================================================

# 1. Core Worker Function: Executes triplex search and sequence segmentation
run_triplex_analysis <- function(fasta_path, output_dir, 
                                   run_min_score = 1,
                                   run_min_len = 6,
                                   run_p_value = 0.05,
                                   run_min_loop = 6,    
                                   run_max_loop = 90) {  
  
  if (!requireNamespace("Biostrings", quietly = TRUE)) library(Biostrings)
  if (!requireNamespace("triplex", quietly = TRUE)) library(triplex)
  
  # Helper: Sanitizes strings for filesystem compatibility
  make_safe_filename <- function(name) {
    safe_name <- gsub("[\\\\/:*?\"<>|]", "_", name)
    return(substr(safe_name, 1, 100))
  }
  
  # Helper: Extracts genomic coordinates and metadata from FASTA headers
  parse_sequence_id <- function(seq_id) {
    parts <- strsplit(seq_id, "\\|")[[1]]
    if (length(parts) >= 3) {
      seq_id_safe <- make_safe_filename(seq_id)
      pos_parts <- strsplit(parts[3], "-")[[1]]
      if (length(pos_parts) == 2) {
        start_pos <- as.numeric(pos_parts[1]); end_pos <- as.numeric(pos_parts[2])
      } else {
        start_pos <- 1; end_pos <- nchar(seq_id)
      }
      return(list(safe_id = seq_id_safe, chromosome = parts[2], start = start_pos, end = end_pos))
    } else {
      return(list(safe_id = make_safe_filename(seq_id), chromosome = "unknown", start = 1, end = nchar(seq_id)))
    }
  }
  
  # Workspace initialization
  output_dir <- normalizePath(output_dir, mustWork = FALSE, winslash = "\\")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  summary_file <- file.path(output_dir, "analysis_summary.csv")
  full_details_file <- file.path(output_dir, "all_triplex_sequences_detailed.csv")
  
  fasta_data <- tryCatch({ readDNAStringSet(fasta_path) }, error = function(e) NULL)
  if (is.null(fasta_data)) return(paste("Read failed:", basename(fasta_path)))
  
  # Initialize CSV headers for persistent storage
  if (!file.exists(summary_file)) {
    write.csv(data.frame(
      SequenceID = character(), SafeID = character(), Chromosome = character(),
      Start = numeric(), End = numeric(), SequenceLength = numeric(),
      TriplexCount = numeric(), BestScore = numeric(), AvgScore = numeric(),
      AnalysisStatus = character(), stringsAsFactors = FALSE
    ), summary_file, row.names = FALSE)
  }
  
  if (!file.exists(full_details_file)) {
    write.csv(data.frame(
      SequenceID = character(), TriplexID = numeric(), Sequence = character(),
      Stem1 = character(), Loop = character(), Stem2 = character(), 
      Score = numeric(), PValue = numeric(), Length = numeric(),
      Type = character(), Strand = character(), stringsAsFactors = FALSE
    ), full_details_file, row.names = FALSE)
  }
  
  buffer_size <- 2000
  summary_buffer <- list()
  details_buffer <- list()
  total_sequences <- length(fasta_data)
  
  for (i in 1:total_sequences) {
    seq_id <- names(fasta_data)[i]
    seq_str <- as.character(fasta_data[[i]])
    seq_info <- parse_sequence_id(seq_id)
    
    # Executes the triplex.search algorithm with user-defined constraints
    result <- tryCatch({
      triplex.search(DNAString(seq_str), 
                       min_score = run_min_score, 
                       min_len   = run_min_len, 
                       p_value   = run_p_value,
                       min_loop  = run_min_loop, 
                       max_loop  = run_max_loop) 
    }, error = function(e) NULL)
    
    if (!is.null(result) && length(result) > 0) {
      # 1. Aggregates summary statistics
      summary_buffer[[length(summary_buffer)+1]] <- data.frame(
        SequenceID = seq_id, SafeID = seq_info$safe_id, Chromosome = seq_info$chromosome,
        Start = seq_info$start, End = seq_info$end, SequenceLength = nchar(seq_str),
        TriplexCount = length(result), BestScore = max(score(result)), AvgScore = mean(score(result)),
        AnalysisStatus = "Success", stringsAsFactors = FALSE
      )
      
      # 2. Structural parsing: Segments triplexes into Stem1, Loop, and Stem2
      full_seqs <- as.character(result)
      n_hits <- length(result)
      stem1_vec <- rep("", n_hits)
      loop_vec  <- rep("", n_hits)
      stem2_vec <- rep("", n_hits)
      
      for (k in 1:n_hits) {
        tryCatch({
          t_start <- start(result)[k]
          l_start_abs <- lstart(result)[k]
          l_end_abs   <- lend(result)[k]
          
          rel_l_start <- l_start_abs - t_start + 1
          rel_l_end   <- l_end_abs - t_start + 1
          full_len    <- nchar(full_seqs[k])
          
          if (rel_l_start > 1) stem1_vec[k] <- substr(full_seqs[k], 1, rel_l_start - 1)
          if (rel_l_end >= rel_l_start) loop_vec[k] <- substr(full_seqs[k], rel_l_start, rel_l_end)
          if (rel_l_end < full_len) stem2_vec[k] <- substr(full_seqs[k], rel_l_end + 1, full_len)
          
        }, error = function(e) {
          stem1_vec[k] <<- "Error"
        })
      }
      
      details_buffer[[length(details_buffer)+1]] <- data.frame(
        SequenceID = seq_id, TriplexID = 1:n_hits, Sequence = full_seqs,
        Stem1 = stem1_vec, Loop = loop_vec, Stem2 = stem2_vec,
        Score = score(result), PValue = pvalue(result), 
        Length = width(result), Type = as.character(type(result)), 
        Strand = as.character(strand(result)), stringsAsFactors = FALSE
      )
    } else {
      # Handles sequences with zero triplex findings
      summary_buffer[[length(summary_buffer)+1]] <- data.frame(
        SequenceID = seq_id, SafeID = seq_info$safe_id, Chromosome = seq_info$chromosome,
        Start = seq_info$start, End = seq_info$end, SequenceLength = nchar(seq_str),
        TriplexCount = 0, BestScore = NA, AvgScore = NA, AnalysisStatus = "No_Triplex", stringsAsFactors = FALSE
      )
    }
    
    # Periodically flushes memory-resident buffers to disk to ensure data persistence
    if (i %% buffer_size == 0 || i == total_sequences) {
      if (length(summary_buffer) > 0) {
        write.table(do.call(rbind, summary_buffer), summary_file, append = TRUE, sep = ",", 
                    row.names = FALSE, col.names = FALSE, qmethod = "double")
        summary_buffer <- list()
      }
      if (length(details_buffer) > 0) {
        write.table(do.call(rbind, details_buffer), full_details_file, append = TRUE, sep = ",", 
                    row.names = FALSE, col.names = FALSE, qmethod = "double")
        details_buffer <- list()
      }
      gc()
    }
  }
  return(paste("Completed:", basename(fasta_path)))
}

# 2. Parallel Processing Manager: Orchestrates batch tasks with resume functionality
batch_process_parallel <- function(input_dir, output_base_dir, cores_to_use, params) {
  library(parallel)
  if (!dir.exists(output_base_dir)) dir.create(output_base_dir, recursive = TRUE)
  
  all_fasta_files <- list.files(input_dir, pattern = "\\.(fasta|fa|fna)$", full.names = TRUE, ignore.case = TRUE)
  if (length(all_fasta_files) == 0) stop("No FASTA files found")
  
  cat("======================================================\n")
  cat("Total files detected:", length(all_fasta_files), "\n")
  
  # Checkpoints logic: Skips folders that already contain valid analysis output
  files_to_process <- c()
  for (f in all_fasta_files) {
    file_name <- tools::file_path_sans_ext(basename(f))
    expected_out_dir <- file.path(output_base_dir, paste0("analysis_", file_name))
    
    if (!dir.exists(expected_out_dir) || length(list.files(expected_out_dir, pattern = "\\.csv$")) == 0) {
      files_to_process <- c(files_to_process, f)
    }
  }
  
  n_total <- length(all_fasta_files)
  n_todo <- length(files_to_process)
  n_done <- n_total - n_todo
  
  cat("Completed (Skipped):", n_done, "\n")
  cat("Pending (To process):", n_todo, "\n")
  
  if (n_todo == 0) {
    cat("All tasks already completed. Pipeline terminated. \n")
    return(NULL)
  }
  
  cat("------------------------------------------------------\n")
  cat("Launching cluster with", cores_to_use, "workers...\n")
  cat("======================================================\n")
  
  cl <- makeCluster(cores_to_use)
  clusterEvalQ(cl, { library(Biostrings); library(triplex); gc() })
  clusterExport(cl, c("run_triplex_analysis"), envir = environment())
  clusterExport(cl, list("params", "output_base_dir"), envir = environment())
  
  run_wrapper <- function(file_path) {
    gc()
    file_name <- tools::file_path_sans_ext(basename(file_path))
    out_dir <- file.path(output_base_dir, paste0("analysis_", file_name))
    tryCatch({
      run_triplex_analysis(file_path, output_dir = out_dir, 
                           run_min_score = params$min_score,
                           run_min_len   = params$min_len,
                           run_p_value   = params$p_value,
                           run_min_loop  = params$min_loop,
                           run_max_loop  = params$max_loop)
    }, error = function(e) paste("Error:", e$message))
  }
  
  results <- parLapply(cl, files_to_process, run_wrapper)
  stopCluster(cl)
  gc()
  return(results)
}

# 3. Entry Point: Defines parameters and executes the batch pipeline
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript triplex_parallel_search_pipeline.R <input_dir> <output_base_dir>")
  }

  input_dir <- args[1]
  output_base_dir <- args[2]
  
  # User-defined search constraints maintained for reproducibility
  run_params <- list(
    min_score = 1,      
    min_len   = 6,      
    p_value   = 0.05,      
    min_loop  = 6,      
    max_loop  = 90      
  )
  
  results <- batch_process_parallel(input_dir, output_base_dir, cores_to_use = 6, params = run_params)
  print(results)
}

main()