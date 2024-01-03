#!/usr/bin/env python

import subprocess
import os
import sys
import pyhmmer
from pyhmmer import hmmscan as hmmscan
import pandas as pd
import multiprocessing.pool
import math
import re
import time

input_dir = sys.argv[1]

out_dir = sys.argv[2]

which_DB = sys.argv[3]

CPUcount = sys.argv[4]

evalue_cut = sys.argv[5]

evalue_cut = float(evalue_cut)

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)


starttime = time.perf_counter()

def hmmscanner(seqs):
    scanout = list(hmmscan(pyhmmer.easel.SequenceFile(seqs, digital=True), pyhmmer.plan7.HMMFile(which_DB)))
    return scanout



splitAA_list = []
for splitAA in os.listdir(input_dir):
    if splitAA.endswith('.faa'):
        f = os.path.join(input_dir, splitAA)

        if os.path.isfile(f) and os.path.getsize(f) > 0:
            splitAA_list.append(f)

if not splitAA_list:
    print("no files found for pyhmmer in " + str(input_dir))
    exit

## get lengths of each hmm
## based on code from https://github.com/althonos/pyhmmer/issues/27#issuecomment-1713131288
hmm_lengths = {}
with pyhmmer.plan7.HMMFile(which_DB) as hmm_file:
    for hmm in hmm_file:
        hmm_lengths[hmm.name.decode()] = len(hmm.consensus)

hmmscan_list = []
with multiprocessing.pool.ThreadPool(int(CPUcount)) as pool:
    for alignments in pool.map(hmmscanner, splitAA_list):
        for model in alignments:
            quer1 = model.query_name.decode()
            pos = quer1.rfind("_")
            contig = quer1[:pos]
            for hit in model:
                target_name = hit.name.decode()
                target_acc = hit.accession
                full_seq_evalue = hit.evalue
                seq_pvalue = hit.pvalue      
                n_aligned_positions = len(
                    hit.best_domain.alignment.hmm_sequence
                ) - hit.best_domain.alignment.hmm_sequence.count(".")
                hmm_coverage = (
                    n_aligned_positions / hmm_lengths[hit.best_domain.alignment.hmm_name.decode()]
                )
                
                hmmscan_list.append([quer1, contig, target_name, full_seq_evalue, seq_pvalue, 
                                     n_aligned_positions, hmm_coverage])

hmmscan_pools_df = pd.DataFrame(hmmscan_list, columns=["ORFquery", "contig", "target", "evalue", 
                                                       "pvalue", "n_aligned_positions", "hmm_coverage"])\
    .query("evalue <= 0.1").query("hmm_coverage >= 0.8 | evalue <= @evalue_cut")\
    .sort_values('evalue').drop_duplicates('ORFquery')

if not hmmscan_pools_df.empty:
    hmmscan_output_file = os.path.join(out_dir, "pyhmmer_report_AAs.tsv")

    hmmscan_pools_df.to_csv(hmmscan_output_file,
                            sep = "\t", index = False)


hmmscan_contig_sum = hmmscan_pools_df.groupby("contig").size().reset_index(name='count')

if not hmmscan_contig_sum.empty:
    contig_sum_file = os.path.join(out_dir, "contig_hit_count.tsv")

    hmmscan_contig_sum.to_csv(contig_sum_file,
                            sep = "\t", index = False)

endtime = time.perf_counter()

time_taken = endtime - starttime

print(f"pyhmmscan of {os.path.basename(which_DB)} finished in " + "%.2f" % time_taken + " seconds")