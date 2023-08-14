#!/usr/bin/env python

import os
import sys
import pandas as pd
import math
import re
import numpy as np



### process hallmark taxonomy for ORF calling

mmseqs2_tax_table = sys.argv[1]

repeat_table = sys.argv[2]

out_dir = sys.argv[3]

tax_df = pd.read_csv(mmseqs2_tax_table, header = None, sep = "\t",
                     names = ["query","target","pident","alnlen","evalue","theader","taxlineage"])\
                     .sort_values('evalue').drop_duplicates('query').query("evalue <= 1e-3")

tax_df['ORFcaller'] = np.where(tax_df['taxlineage']
                               .str.contains(
                                   "Caudoviricetes|Crassvirales|Malgrandaviricetes|Tubulavirales|Leviviricetes|Duplopiviricetes|Kalamavirales|Vinavirales|Autolykiviridae",
                                   case = False), 'phanotate', 'prodigal')

tax_df["pos"] = tax_df["query"].str.rfind("_")

tax_df["contig"] = tax_df.apply(lambda x: x["query"][0:x["pos"]], axis = 1)

##override phanotate call due to maximum length
length_df = pd.read_csv(repeat_table, sep = "\t")[['contig', 'out_length_contig']]

merge_tax_df = length_df.merge(tax_df, on = "contig", how = "left")

merge_tax_df['ORFcaller'] = np.where(merge_tax_df['ORFcaller'].isnull(), 
                                     'prodigal', merge_tax_df['ORFcaller'])


### Forcing long contigs (over 300 kb) into prodigal as phanotate does not handle these well
merge_tax_df['ORFcaller'] = np.where(merge_tax_df['out_length_contig'] >= 300000, 
                                     'prodigal', merge_tax_df['ORFcaller'])

merge_tax_df['Note'] = np.where(merge_tax_df['out_length_contig'] >= 300000, 
                                     'over phanotate length limit', 'NA')


### I don't really need to make this file, just for dev
tax_label_file = os.path.join(out_dir, "orf_caller_each_seq.tsv")

merge_tax_df.to_csv(tax_label_file, sep = "\t", index = False)

ORFcaller_majority = merge_tax_df.groupby("contig")['ORFcaller'].agg(pd.Series.mode).to_frame()


##save lists to files
prodigal_seqs = ORFcaller_majority.query("ORFcaller == 'prodigal'")

if not prodigal_seqs.empty:
    prodigal_file = os.path.join(out_dir, "prodigal_seqs1.txt")
    if os.path.isfile(prodigal_file):
        os.remove(prodigal_file)

    for i in prodigal_seqs.index:

        print(i, file = open(prodigal_file, "a"))


phanotate_seqs = ORFcaller_majority.query("ORFcaller == 'phanotate'")

if not phanotate_seqs.empty:
    phanotate_file = os.path.join(out_dir, "phanotate_seqs1.txt")
    if os.path.isfile(phanotate_file):
        os.remove(phanotate_file)

    for i in phanotate_seqs.index:
        print(i, file = open(phanotate_file, "a"))




