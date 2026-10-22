# Purpose: Generate feature-matched negative samples with chromosome-level physical isolation for H-DNA benchmarking.

import os
import random
import bisect
import argparse
import pandas as pd
from pyfaidx import Fasta
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

CHROMS_TRAIN = [f'chr{i}' for i in range(1, 16)]
CHROMS_VAL = ['chr16', 'chr17', 'chr18', 'chr19']
CHROMS_TEST = ['chr20', 'chr21', 'chr22']


def calculate_metrics(seq):
    """Calculate sequence length, GC ratio, and purine (A+G) ratio."""
    seq = seq.upper()
    length = len(seq)
    if length == 0:
        return 0, 0.0, 0.0
    gc = (seq.count('G') + seq.count('C')) / length
    purine = (seq.count('A') + seq.count('G')) / length
    return length, gc, purine


def parse_fasta(fasta_path, desc=""):
    """Parse positive FASTA sequences and compute length and nucleotide compositions."""
    print(f"[INFO] Loading {desc} positive sequences from: {fasta_path}")
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Input file not found: {fasta_path}")

    records = []
    with open(fasta_path, 'r', encoding='utf-8') as f:
        seq_id, seq = "", ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id and seq:
                    length, gc, purine = calculate_metrics(seq)
                    records.append({
                        'sample_id': seq_id, 'sequence': seq, 'length': length,
                        'gc_content': gc, 'purine_content': purine,
                        'chrom': 'N/A', 'start': -1, 'end': -1, 'label': 1
                    })
                seq_id = line[1:]
                seq = ""
            else:
                seq += line
        if seq_id and seq:
            length, gc, purine = calculate_metrics(seq)
            records.append({
                'sample_id': seq_id, 'sequence': seq, 'length': length,
                'gc_content': gc, 'purine_content': purine,
                'chrom': 'N/A', 'start': -1, 'end': -1, 'label': 1
            })

    df = pd.DataFrame(records)
    print(f"[INFO] Successfully loaded {len(df)} records for {desc}.")
    return df


def check_overlap(used_intervals, chrom, start, end):
    """Verify non-overlapping intervals within chromosome coordinates."""
    if chrom not in used_intervals:
        return False
    intervals = used_intervals[chrom]
    idx = bisect.bisect_right([iv[1] for iv in intervals], start)
    if idx < len(intervals):
        if intervals[idx][0] < end:
            return True
    return False


def worker_generate_negatives(args):
    """Worker task: sample background sequences matching length, GC, and purine content."""
    task_data, target_chroms, hg38_path, tolerance = args
    genome = Fasta(hg38_path)

    chrom_sizes = {k: len(genome[k]) for k in target_chroms if k in genome}
    chrom_names = list(chrom_sizes.keys())
    chrom_weights = list(chrom_sizes.values())

    negatives = []
    local_used_intervals = {k: [] for k in target_chroms}

    for row in task_data:
        target_len = row['length']
        target_gc = row['gc_content']
        target_purine = row['purine_content']

        found_seq, final_chrom, final_start, final_end = None, None, -1, -1

        # Phase 1: strict feature matching
        for _ in range(200):
            chrom = random.choices(chrom_names, weights=chrom_weights, k=1)[0]
            max_start = chrom_sizes[chrom] - target_len
            if max_start <= 0:
                continue

            start = random.randint(0, max_start)
            end = start + target_len

            if check_overlap(local_used_intervals, chrom, start, end):
                continue

            try:
                seq = genome[chrom][start:end].seq.upper()
                if 'N' in seq:
                    continue
                _, cand_gc, cand_purine = calculate_metrics(seq)
                if abs(cand_gc - target_gc) <= tolerance and abs(cand_purine - target_purine) <= tolerance:
                    found_seq = seq
                    final_chrom, final_start, final_end = chrom, start, end
                    break
            except Exception:
                continue

        # Phase 2: background sampling fallback if strict criteria fail
        if not found_seq:
            for _ in range(100):
                chrom = random.choices(chrom_names, weights=chrom_weights, k=1)[0]
                max_start = chrom_sizes[chrom] - target_len
                if max_start <= 0:
                    continue
                start = random.randint(0, max_start)
                end = start + target_len
                if check_overlap(local_used_intervals, chrom, start, end):
                    continue
                try:
                    seq = genome[chrom][start:end].seq.upper()
                    if 'N' not in seq:
                        found_seq = seq
                        final_chrom, final_start, final_end = chrom, start, end
                        break
                except Exception:
                    continue

        if found_seq:
            bisect.insort(local_used_intervals[final_chrom], (final_start, final_end))
            _, final_gc, final_purine = calculate_metrics(found_seq)
            negatives.append({
                'sample_id': f"neg_{final_chrom}_{final_start}",
                'sequence': found_seq,
                'length': target_len,
                'gc_content': final_gc,
                'purine_content': final_purine,
                'chrom': final_chrom,
                'start': final_start,
                'end': final_end,
                'label': 0
            })

    return negatives


def resolve_global_overlaps(neg_df):
    """Deduplicate overlapping candidate intervals across parallel workers."""
    if neg_df.empty:
        return neg_df
    neg_df = neg_df.sort_values(by=['chrom', 'start']).reset_index(drop=True)
    drop_indices = set()
    for chrom in neg_df['chrom'].unique():
        chrom_data = neg_df[neg_df['chrom'] == chrom]
        prev_end = -1
        for idx, row in chrom_data.iterrows():
            if row['start'] < prev_end:
                drop_indices.add(idx)
            else:
                prev_end = row['end']

    if drop_indices:
        neg_df = neg_df.drop(index=list(drop_indices)).reset_index(drop=True)
    return neg_df


def generate_negatives_parallel(pos_df, target_chroms, desc_name, num_workers, hg38_path, tolerance, output_dir):
    """Execute multi-processing pipeline to sample negative genomic regions."""
    print(f"\n[INFO] Generating {desc_name} negative instances from {target_chroms[0]} to {target_chroms[-1]}...")
    pos_dicts = pos_df.to_dict('records')

    chunk_size = 500
    chunks = [
        [pos_dicts[i:i + chunk_size], target_chroms, hg38_path, tolerance]
        for i in range(0, len(pos_dicts), chunk_size)
    ]

    all_negatives = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for res in tqdm(executor.map(worker_generate_negatives, chunks), total=len(chunks), desc=f"{desc_name} Workers"):
            all_negatives.extend(res)

    neg_df = pd.DataFrame(all_negatives)
    neg_df = resolve_global_overlaps(neg_df)

    pure_neg_path = os.path.join(output_dir, f"pure_negative_{desc_name.lower()}.csv")
    neg_df.to_csv(pure_neg_path, index=False)
    print(f"[INFO] Saved {len(neg_df)} pure negatives to {pure_neg_path}")

    return neg_df


def mix_and_save_dataset(pos_df, neg_df, filename, output_dir, random_state=42):
    """Align positive and negative subsets 1:1, shuffle, and export to CSV."""
    min_len = min(len(pos_df), len(neg_df))
    pos_aligned = pos_df.iloc[:min_len].copy()
    neg_aligned = neg_df.iloc[:min_len].copy()

    mixed_df = pd.concat([pos_aligned, neg_aligned], ignore_index=True)
    mixed_df = mixed_df.sample(frac=1, random_state=random_state).reset_index(drop=True)

    out_path = os.path.join(output_dir, filename)
    mixed_df.to_csv(out_path, index=False)
    print(f"[INFO] Exported {filename} (Total: {len(mixed_df)}, Positives: {min_len}, Negatives: {min_len})")


def main():
    parser = argparse.ArgumentParser(
        description="Dataset constructor with chromosome-level isolation and feature-matched background sampling."
    )
    parser.add_argument("--pos_train", type=str, required=True,
                        help="Path to positive FASTA for training (e.g., chr1-15)")
    parser.add_argument("--pos_val", type=str, required=True,
                        help="Path to positive FASTA for validation (e.g., chr16-19)")
    parser.add_argument("--pos_test", type=str, required=True,
                        help="Path to positive FASTA for testing (e.g., chr20-22)")
    parser.add_argument("--reference_genome", type=str, required=True,
                        help="Path to human reference genome FASTA (e.g., hg38.fa)")
    parser.add_argument("--output_dir", type=str, default="./output_dataset",
                        help="Directory to save final datasets and intermediate files")
    parser.add_argument("--num_workers", type=int, default=8,
                        help="Number of parallel worker processes (default: 8)")
    parser.add_argument("--tolerance", type=float, default=0.05,
                        help="Tolerance threshold for GC and purine matching (default: 0.05)")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Parse pre-split positive sequences
    pos_train = parse_fasta(args.pos_train, "Train")
    pos_val = parse_fasta(args.pos_val, "Validation")
    pos_test = parse_fasta(args.pos_test, "Test")

    # 2. Parallel background sampling with physical chromosome partitioning
    neg_train = generate_negatives_parallel(
        pos_train, CHROMS_TRAIN, "Train", args.num_workers, args.reference_genome, args.tolerance, args.output_dir
    )
    neg_val = generate_negatives_parallel(
        pos_val, CHROMS_VAL, "Val", args.num_workers, args.reference_genome, args.tolerance, args.output_dir
    )
    neg_test = generate_negatives_parallel(
        pos_test, CHROMS_TEST, "Test", args.num_workers, args.reference_genome, args.tolerance, args.output_dir
    )

    # 3. 1:1 balance alignment, shuffle, and serialization
    print("\n[INFO] Assembling balanced datasets...")
    mix_and_save_dataset(pos_train, neg_train, "final_dataset_train.csv", args.output_dir)
    mix_and_save_dataset(pos_val, neg_val, "final_dataset_val.csv", args.output_dir)
    mix_and_save_dataset(pos_test, neg_test, "final_dataset_test.csv", args.output_dir)

    print("\n[INFO] Dataset pipeline execution completed successfully.")


if __name__ == "__main__":
    main()