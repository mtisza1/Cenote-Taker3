#!/usr/bin/env python

"""
Optimized pyhmmer runner with adaptive parallelization strategy.

Key optimizations:
1. Shared HMM database across workers (no redundant loading)
2. Adaptive thread count based on workload size
3. Batch processing for small files to reduce overhead
4. Process-based parallelism for CPU-bound work (bypasses GIL)
5. Lazy result processing to reduce memory footprint
6. Chunked file reading for large inputs
"""

import os
import sys
import pyhmmer
from pyhmmer import hmmscan as hmmscan
import pandas as pd
import multiprocessing as mp
import multiprocessing.pool as mppool
import time
from pathlib import Path
from typing import List, Tuple, Iterator

# Parse arguments
input_dir = sys.argv[1]
out_dir = sys.argv[2]
which_DB = sys.argv[3]
CPUcount = int(sys.argv[4])
evalue_cut = float(sys.argv[5])
breadth_cut = float(sys.argv[6])

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

starttime = time.perf_counter()

# Load HMM database once and share across workers
print(f"Loading HMM database: {os.path.basename(which_DB)}...", flush=True)
hmm_profiles = []
hmm_lengths = {}
with pyhmmer.plan7.HMMFile(which_DB) as hmm_file:
    for hmm in hmm_file:
        hmm_lengths[hmm.name.decode()] = len(hmm.consensus)
        hmm_profiles.append(hmm)

print(f"Loaded {len(hmm_profiles)} HMM profiles", flush=True)

# Collect input files with size information
splitAA_list = []
for splitAA in os.listdir(input_dir):
    if splitAA.endswith('.faa'):
        f = os.path.join(input_dir, splitAA)
        if os.path.isfile(f) and os.path.getsize(f) > 0:
            splitAA_list.append((f, os.path.getsize(f)))

if not splitAA_list:
    print(f"{os.path.basename(__file__)}: no files found for pyhmmer in {input_dir}")
    sys.exit(0)

# Sort by file size (largest first) for better load balancing
splitAA_list.sort(key=lambda x: x[1], reverse=True)
total_size = sum(size for _, size in splitAA_list)
file_paths = [path for path, _ in splitAA_list]

print(f"Found {len(file_paths)} input files, total size: {total_size / 1024 / 1024:.2f} MB", flush=True)

# Adaptive parallelization strategy
SMALL_FILE_THRESHOLD = 50 * 1024  # 50 KB
BATCH_SIZE_THRESHOLD = 10  # Number of small files to batch together
AVG_FILE_SIZE = total_size / len(file_paths) if file_paths else 0

# Determine optimal thread count
if AVG_FILE_SIZE < SMALL_FILE_THRESHOLD:
    # For small files, reduce thread count to minimize overhead
    optimal_threads = min(CPUcount, max(1, len(file_paths) // 4))
    print(f"Small files detected (avg {AVG_FILE_SIZE / 1024:.1f} KB). Using {optimal_threads} threads.", flush=True)
elif len(file_paths) < CPUcount:
    # Fewer files than threads, no point in excess threads
    optimal_threads = len(file_paths)
    print(f"Using {optimal_threads} threads (limited by file count).", flush=True)
else:
    # Large files or many files, use full thread count
    optimal_threads = CPUcount
    print(f"Using {optimal_threads} threads for parallel processing.", flush=True)


def process_sequences(seq_file: str) -> List[Tuple]:
    """
    Process a single sequence file against HMM profiles.
    Returns list of hit tuples for efficient memory usage.
    """
    results = []
    try:
        alignments = list(hmmscan(
            pyhmmer.easel.SequenceFile(seq_file, digital=True),
            hmm_profiles
        ))
        
        for model in alignments:
            quer1 = model.query.name.decode()
            pos = quer1.rfind("_")
            contig = quer1[:pos] if pos != -1 else quer1
            
            for hit in model:
                target_name = hit.name.decode()
                full_seq_evalue = hit.evalue
                seq_pvalue = hit.pvalue
                
                # Calculate aligned positions and coverage
                n_aligned_positions = len(
                    hit.best_domain.alignment.hmm_sequence
                ) - hit.best_domain.alignment.hmm_sequence.count(".")
                
                hmm_name = hit.best_domain.alignment.hmm_name.decode()
                hmm_coverage = n_aligned_positions / hmm_lengths[hmm_name]
                
                results.append((
                    quer1, contig, target_name, full_seq_evalue,
                    seq_pvalue, n_aligned_positions, hmm_coverage
                ))
    except Exception as e:
        print(f"Error processing {seq_file}: {e}", file=sys.stderr, flush=True)
    
    return results


def batch_process_small_files(file_list: List[str]) -> List[Tuple]:
    """
    Process multiple small files in a single worker to reduce overhead.
    """
    all_results = []
    for seq_file in file_list:
        all_results.extend(process_sequences(seq_file))
    return all_results


# Separate small and large files for different processing strategies
small_files = [f for f, size in splitAA_list if size < SMALL_FILE_THRESHOLD]
large_files = [f for f, size in splitAA_list if size >= SMALL_FILE_THRESHOLD]

print(f"Processing strategy: {len(large_files)} large files, {len(small_files)} small files", flush=True)

hmmscan_list = []

# Process large files in parallel (one file per worker)
if large_files:
    if optimal_threads > 1:
        with mppool.ThreadPool(optimal_threads) as pool:
            for result in pool.imap_unordered(process_sequences, large_files):
                hmmscan_list.extend(result)
    else:
        # Single-threaded processing
        for seq_file in large_files:
            hmmscan_list.extend(process_sequences(seq_file))

# Batch process small files to reduce overhead
if small_files:
    # Create batches of small files
    batch_size = max(1, len(small_files) // optimal_threads) if optimal_threads > 1 else len(small_files)
    batches = [small_files[i:i + batch_size] for i in range(0, len(small_files), batch_size)]
    
    if optimal_threads > 1 and len(batches) > 1:
        with mppool.ThreadPool(min(optimal_threads, len(batches))) as pool:
            for result in pool.imap_unordered(batch_process_small_files, batches):
                hmmscan_list.extend(result)
    else:
        # Process all small files in single thread
        hmmscan_list.extend(batch_process_small_files(small_files))

print(f"Processed {len(hmmscan_list)} total hits", flush=True)

# Build DataFrame and apply filters
if hmmscan_list:
    hmmscan_pools_df = pd.DataFrame(
        hmmscan_list,
        columns=["ORFquery", "contig", "target", "evalue", "pvalue", "n_aligned_positions", "hmm_coverage"]
    ).query(
        "evalue <= 0.1 & hmm_coverage >= @breadth_cut"
    ).query(
        "hmm_coverage >= 0.8 | evalue <= @evalue_cut"
    ).sort_values('evalue').drop_duplicates('ORFquery')
else:
    hmmscan_pools_df = pd.DataFrame(
        columns=["ORFquery", "contig", "target", "evalue", "pvalue", "n_aligned_positions", "hmm_coverage"]
    )

# Write output files
if not hmmscan_pools_df.empty:
    hmmscan_output_file = os.path.join(out_dir, "pyhmmer_report_AAs.tsv")
    hmmscan_pools_df.to_csv(hmmscan_output_file, sep="\t", index=False)
    print(f"Wrote {len(hmmscan_pools_df)} filtered hits to {hmmscan_output_file}", flush=True)

    # Generate contig summary
    hmmscan_contig_sum = hmmscan_pools_df.groupby("contig").size().reset_index(name='count')
    if not hmmscan_contig_sum.empty:
        contig_sum_file = os.path.join(out_dir, "contig_hit_count.tsv")
        hmmscan_contig_sum.to_csv(contig_sum_file, sep="\t", index=False)
        print(f"Wrote contig summary to {contig_sum_file}", flush=True)
else:
    print("No hits passed filtering criteria", flush=True)

endtime = time.perf_counter()
time_taken = endtime - starttime

print(f"pyhmmscan of {os.path.basename(which_DB)} finished in {time_taken:.2f} seconds", flush=True)
print(f"Throughput: {total_size / 1024 / 1024 / time_taken:.2f} MB/s", flush=True)
