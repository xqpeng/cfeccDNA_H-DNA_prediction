#!/user/bin/env python3
# -*- coding: utf-8 -*-
import os
import glob
import pandas as pd
import multiprocessing as mp
import time
import argparse

# ================= Configuration Area =================
# Number of CPU cores for parallel processing
NUM_CORES = 16

# Input filename (detailed triplex sequence records)
INPUT_FILENAME = "all_triplex_sequences_detailed.csv"

# Output filename (filtered records with the best score per sequence)
OUTPUT_FILENAME = "best_triplex_sequences_detailed.csv"
# ======================================================

def process_single_folder(folder_path):
    """
    Worker function: Processes a single sample folder.
    1. Reads the detailed CSV file.
    2. Groups data by SequenceID and selects the entry with the highest Score.
    3. Saves the filtered results to a new CSV file.
    """
    try:
        input_file = os.path.join(folder_path, INPUT_FILENAME)
        output_file = os.path.join(folder_path, OUTPUT_FILENAME)

        # 1. Check if the input file exists
        if not os.path.exists(input_file):
            return f"⚠️ Skipped (File missing): {os.path.basename(folder_path)}"

        # 2. Read CSV with robust encoding handling
        # Pandas handles quotechar='"' by default for quoted data
        try:
            df = pd.read_csv(input_file, encoding='utf-8')
        except UnicodeDecodeError:
            # Fallback to GBK encoding if UTF-8 fails
            df = pd.read_csv(input_file, encoding='gbk')

        if df.empty:
            return f"⚠️ Skipped (File empty): {os.path.basename(folder_path)}"

        # 3. Data Cleaning and Filtering
        # Ensure 'Score' column is treated as numeric for comparison
        if 'Score' not in df.columns:
            return f"❌ Error (Missing Score column): {os.path.basename(folder_path)}"

        df['Score'] = pd.to_numeric(df['Score'], errors='coerce')

        # Core Logic:
        # A. Sort by Score in descending order (highest score first)
        # B. Drop duplicates based on SequenceID, keeping only the first entry
        df_sorted = df.sort_values(by='Score', ascending=False)
        df_best = df_sorted.drop_duplicates(subset=['SequenceID'], keep='first')

        # 4. Save the filtered results
        df_best.to_csv(output_file, index=False, encoding='utf-8')

        # Statistics for logging
        original_count = len(df)
        filtered_count = len(df_best)

        return f"✅ Success: {os.path.basename(folder_path)} ({original_count} -> {filtered_count})"

    except Exception as e:
        return f"❌ Failed to process {os.path.basename(folder_path)}: {str(e)}"


def main():
    parser = argparse.ArgumentParser(description="Filter best score triplex sequences.")
    parser.add_argument("root_dir", help="Root directory containing the batch analysis results")
    args = parser.parse_args()

    ROOT_DIR = args.root_dir

    start_time = time.time()
    print(f"🚀 Starting Best Score Selection for Triplex sequences...")
    print(f"Root Directory: {ROOT_DIR}")

    # 1. Retrieve all subfolders
    # Assumes each sample is located in a first-level subdirectory under ROOT_DIR
    if not os.path.exists(ROOT_DIR):
        print(f"❌ Error: Directory not found: {ROOT_DIR}")
        return

    subfolders = [f.path for f in os.scandir(ROOT_DIR) if f.is_dir()]

    if not subfolders:
        print("❌ No subfolders found. Please check the root directory path.")
        return

    print(f"Sample folders found: {len(subfolders)}")
    print(f"Utilizing CPU cores: {NUM_CORES}")
    print("-" * 50)

    # 2. Parallel execution using Multiprocessing Pool
    with mp.Pool(NUM_CORES) as pool:
        results = pool.map(process_single_folder, subfolders)

    # 3. Generate summary report
    success_count = 0
    error_count = 0

    print("\nProcessing Logs:")
    for res in results:
        # Print only errors or warnings to maintain a clean console output
        if "❌" in res or "⚠️" in res:
            print(res)
            error_count += 1
        else:
            success_count += 1

    end_time = time.time()
    print("-" * 50)
    print(f"🎉 Task Completed!")
    print(f"Successfully processed: {success_count} folders")
    print(f"Errors/Skipped: {error_count} folders")
    print(f"Output files generated: {OUTPUT_FILENAME}")
    print(f"Total elapsed time: {end_time - start_time:.2f} seconds")


if __name__ == "__main__":
    # Required for Windows compatibility with multiprocessing
    mp.freeze_support()
    main()