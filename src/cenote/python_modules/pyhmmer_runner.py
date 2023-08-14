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


if not os.path.isdir(out_dir):
    os.makedirs(out_dir)




#print("pyhmmscan pool part")

starttime = time.time()

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
                hmmscan_list.append([quer1, contig, target_name, full_seq_evalue, seq_pvalue])

hmmscan_pools_df = pd.DataFrame(hmmscan_list, columns=["ORFquery", "contig", "target", "evalue", "pvalue"])\
    .sort_values('evalue').drop_duplicates('ORFquery').query("evalue <= 1e-8")

if not hmmscan_pools_df.empty:
    hmmscan_output_file = os.path.join(out_dir, "pyhmmer_report_AAs.tsv")

    hmmscan_pools_df.to_csv(hmmscan_output_file,
                            sep = "\t", index = False)


hmmscan_contig_sum = hmmscan_pools_df.groupby("contig").size().reset_index(name='count')

if not hmmscan_contig_sum.empty:
    contig_sum_file = os.path.join(out_dir, "contig_hit_count.tsv")

    hmmscan_contig_sum.to_csv(contig_sum_file,
                            sep = "\t", index = False)

endtime = time.time()

time_taken = endtime - starttime

print("pyhmmscan part took: " + "%.2f" % time_taken + " seconds")