#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pandas as pd
import re
import os
import random
import string
import sys

final_contig_file = sys.argv[1]

final_tax_file = sys.argv[2]

repeat_file = sys.argv[3]

temp_dir = sys.argv[4]

out_dir = sys.argv[5]

#phanotate seqs
phanotate_list_file = os.path.join(temp_dir, "hallmark_tax", "phanotate_seqs1.txt")

#prodigal gcode table
prodigal_list_file = os.path.join(temp_dir, "reORF", "prod_split", "contig_gcodes1.txt")

# make out dir

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

# load dataframes

# tax
tax_call_df = pd.read_csv(final_tax_file, sep = "\t")

# phanotate and prodigal
if os.path.isfile(phanotate_list_file) and os.path.getsize(phanotate_list_file) > 0:
    phan_df = pd.read_csv(phanotate_list_file, header = None, names = ['contig'])
    phan_df['gcode'] = 11
else:
    phan_df = pd.DataFrame()

if os.path.isfile(prodigal_list_file) and os.path.getsize(prodigal_list_file) > 0:
    prod_df = pd.read_csv(prodigal_list_file, header = None, sep = "\t", names = ['contig', 'gcode'])
else:
    prod_df = pd.DataFrame()

## combine phanotate and prodigal table
gcode_list = []

for df in phan_df, prod_df:
    if not df.empty:
        gcode_list.append(df)


try:
    gcode_df = pd.concat(gcode_list, ignore_index=True)
except:
    print("nope")

# repeat and length
repeat_df = pd.read_csv(repeat_file, sep = "\t")

#### loop each virus
for seq_record in SeqIO.parse(final_contig_file, "fasta"):

    #make a random alphanumeric ID 5 characters in length
    randID = ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))

    #chunked/pruned seqs need different handling
    if "@" in seq_record.id:

        nameq = seq_record.id.split("@")[0]
        chunkq = seq_record.id.split("@")[1]
        try:
            organism = tax_call_df.query("contig == @nameq & chunk_name == @chunkq")['taxon'].agg(pd.Series.mode)[0]
        except:
            organism = "unclassified virus"

        try:
            lineage = tax_call_df.query("contig == @nameq & chunk_name == @chunkq")['taxonomy_hierarchy'].agg(pd.Series.mode)[0]
        except:
            lineage = "no lineage"

        gcode = gcode_df.query("contig == @nameq")['gcode'].agg(pd.Series.mode)[0]

        topology = "linear"
    else:
        try:
            organism = tax_call_df.query("contig == @seq_record.id")['taxon'].agg(pd.Series.mode)[0]
        except:
            organism = "unclassified virus"

        try:
            lineage = tax_call_df.query("contig == @seq_record.id")['taxonomy_hierarchy'].agg(pd.Series.mode)[0]
        except:
            lineage = "no lineage"

        gcode = gcode_df.query("contig == @seq_record.id")['gcode'].agg(pd.Series.mode)[0]

        top_str = repeat_df.query("contig == @seq_record.id")['dtr_seq'].agg(pd.Series.mode)

        if not top_str.empty:
            topology = "circular"
        else:
            topology = "linear"

    # here's the whole header string
    header = f">{str(seq_record.id)} [organism={organism} sp. ct{randID}] [gcode={gcode}] [topology={topology}] [note: taxonomic lineage {lineage}]"
    #header = ">" + str(seq_record.id) + " [organism=" + organism + " sp. ct" + randID + "] [gcode=" + str(gcode) + "] [topology=" + topology + "] [note: taxonomic lineage " + lineage + "]" 
    
    seq_output_file = os.path.join(out_dir, str(seq_record.id) + ".fsa")

    print(f"{header}\n{seq_record.seq}", file = open(seq_output_file, "a"))