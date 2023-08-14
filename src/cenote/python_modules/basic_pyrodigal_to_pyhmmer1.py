#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import subprocess
import os
import sys
import pyhmmer
from pyhmmer import hmmscan as hmmscan
import pyrodigal
import pandas as pd
import multiprocessing.pool
import math
import re
import time

input_file = sys.argv[1]

out_dir = sys.argv[2]

hmm_db_virion = sys.argv[3]


#out_dir = os.path.join(str(os.chdir), out_dir)

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)


pyrodigal_output_file = re.sub(".fna$", ".faa", input_file)


print("pyrodigal part")

starttime = time.time()

if os.path.isfile(pyrodigal_output_file):
    os.remove(pyrodigal_output_file)

orf_finder = pyrodigal.OrfFinder(meta = True)


for record in SeqIO.parse(input_file, "fasta"):
    for i, pred in enumerate(orf_finder.find_genes(bytes(record.seq))):
        print(f">{record.id}_{i+1}", file = open(pyrodigal_output_file, "a"))
        print(pred.translate(), file = open(pyrodigal_output_file, "a"))

endtime = time.time()

time_taken = endtime - starttime

print("pyrodigal part took: " + str(time_taken))

print("seqkit split part")

starttime = time.time()

mycpus = os.cpu_count()


split_AA_outfile = re.sub(".faa$", "_hmmscan.faa", pyrodigal_output_file)

split_AA_outdir = os.path.join(out_dir, "temp")

if not os.path.isdir(split_AA_outdir):
    os.makedirs(split_AA_outdir)

subprocess.run(['seqkit', 'split', '-j', str(mycpus), '-p', str(mycpus), 
                '-O', str(split_AA_outdir), '-o', 
                str(split_AA_outfile), str(pyrodigal_output_file)])

endtime = time.time()

time_taken = endtime - starttime

print("seqkit part took: " + str(time_taken))

print("pyhmmscan pool part")

starttime = time.time()

def hmmscanner(seqs):
    scanout = list(hmmscan(pyhmmer.easel.SequenceFile(seqs, digital=True), pyhmmer.plan7.HMMFile(hmm_db_virion)))
    return scanout



splitAA_list = []
for splitAA in os.listdir(split_AA_outdir):
    f = os.path.join(split_AA_outdir, splitAA)

    if os.path.isfile(f):
        splitAA_list.append(f)

hmmscan_list = []
with multiprocessing.pool.ThreadPool() as pool:
    for alignments in pool.map(hmmscanner, splitAA_list):
        for model in alignments:
            quer1 = model.query_name.decode()
            for hit in model:
                target_name = hit.name.decode()
                target_acc = hit.accession
                full_seq_evalue = hit.evalue
                seq_pvalue = hit.pvalue        
                hmmscan_list.append([quer1, target_name, full_seq_evalue, seq_pvalue])

hmmscan_pools_df = pd.DataFrame(hmmscan_list, columns=["query", "target", "evalue", "pvalue"])\
    .sort_values('pvalue').drop_duplicates('query')

hmmscan_output_file = re.sub(".fna$", "_hmmscan.filtered_out.tsv", input_file)

hmmscan_output_file = os.path.join(out_dir, hmmscan_output_file)

hmmscan_pools_df.to_csv(hmmscan_output_file,
                        sep = "\t", index = False)

endtime = time.time()

time_taken = endtime - starttime

print("pyhmmscan part took: " + str(time_taken))