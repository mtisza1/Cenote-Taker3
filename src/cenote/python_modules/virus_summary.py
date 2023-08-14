#!/usr/bin/env python

import os
import sys
import pandas as pd
import math
import re
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

# length/repeat file
length_table = sys.argv[1]
# ct name/original name file
name_table = sys.argv[2]
# gene to contig file
gene_to_contig_table = sys.argv[3]
# taxonomy file
tax_table = sys.argv[4]

sequin_dir = sys.argv[5]

run_title = sys.argv[6]

try:
    main_annot_df = pd.read_csv(gene_to_contig_table, sep = "\t")

    length_df = pd.read_csv(length_table, sep = "\t")

    name_df = pd.read_csv(name_table, sep = "\t", names=['contig', 'input_name'])

    tax_df = pd.read_csv(tax_table, sep = "\t")

except:
    print("couldn't load files for summary")


## merge all files
merge_df = pd.merge(main_annot_df, length_df, on = ["contig", "dtr_seq"], how = "left")

merge_df = pd.merge(merge_df, name_df, on = "contig", how = "left")

merge_df = pd.merge(merge_df, tax_df, on = ["contig", "chunk_name"], how = "left")


merge_df['taxon'] = merge_df['taxon'].fillna("unclassified virus")


## get descriptions from fastas
finalseq_list = []
for fsa in os.listdir(sequin_dir):
    if fsa.endswith('.fsa'):
        f = os.path.join(sequin_dir, fsa)

        if os.path.isfile(f) and os.path.getsize(f) > 0:
            finalseq_list.append(f)

if not finalseq_list:
    print("no files found for seqIO parse " + str(sequin_dir))
    exit


desc_list = []
for seq_file in finalseq_list:
    seq_record = SeqIO.read(seq_file, "fasta")
    try:
        if "@" in seq_record.id:
            contig = seq_record.id.split("@")[0]
            chunkq = seq_record.id.split("@")[1]
        else:
            contig = seq_record.id
            chunkq = None
        fields = re.findall(r'\[.*?\]', seq_record.description)
        organism = re.search(r'\[organism=(.*?)\]', fields[0]).group(1)
        gcode = re.search(r'\[gcode=(.*?)\]', fields[1]).group(1)
        desc_list.append([contig, chunkq, organism, gcode])
    except:
        print("except")

desc_df = pd.DataFrame(desc_list, columns=["contig", "chunk_name", "organism", "genetic_code"])

org_info_df = pd.merge(merge_df, desc_df, on = ["contig", "chunk_name"], how = "left")


## make summary of virus seqs
grouped_df = org_info_df.groupby(['contig', 'contig_length', 'dtr_seq', 'chunk_name', 'chunk_length',
                     'itr_seq', 'input_name', 'taxon', 'taxonomy_hierarchy', 'taxon_level',
                     'avg_hallmark_AAI_to_ref', 'organism', 'genetic_code'], dropna = False)

summary_list = []
for name, group in grouped_df:
    gene_count = group['gene_name'].nunique()
    hallmark_count = group.query("Evidence_source == 'hallmark_hmm'")['gene_name'].nunique()
    hallmark_list = '|'.join(
        list(group.query("Evidence_source == 'hallmark_hmm'")['evidence_description'])
        ).replace("-", " ")
    if name[2]:
        end_type = "DTR"
    elif name[5]:
        end_type = "ITR"
    else:
        end_type = "None"
        
    if gene_count >= 1:
        summary_list.append([name[0], name[6], name[11], name[1], end_type, gene_count, hallmark_count, hallmark_list, name[8]])

summary_df = pd.DataFrame(summary_list, columns=['contig', 'input_name', 'organism', 
                                                 'contig_length', 'end_feature', 'gene_count', 'hallmark_count', 'hallmark_genes', 'taxonomy_hierarchy'])

parentpath = Path(sequin_dir).parents[0]

summary_out = os.path.join(parentpath, f"{run_title}_virus_summary.tsv")

summary_df.to_csv(summary_out, sep = "\t", index = False)
