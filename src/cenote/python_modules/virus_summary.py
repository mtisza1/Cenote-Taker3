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


# ct name/original name file
name_table = sys.argv[1]
# gene to contig file
gene_to_contig_table = sys.argv[2]
# taxonomy file
tax_table = sys.argv[3]
# directory for sequin related files
sequin_dir = sys.argv[4]
# run title
run_title = sys.argv[5]
# prodigal gcodes file
prod_gcodes = sys.argv[6]
# phanotate file
phan_file = sys.argv[7]

try:
    main_annot_df = pd.read_csv(gene_to_contig_table, sep = "\t")

    name_df = pd.read_csv(name_table, sep = "\t", names=['contig', 'input_name'])

    tax_df = pd.read_csv(tax_table, sep = "\t")

except:
    print("couldn't load files for summary")

# phanotate and prodigal
if os.path.isfile(phan_file) and os.path.getsize(phan_file) > 0:
    phan_df = pd.read_csv(phan_file, header = None, names = ['contig'])
    phan_df['gcode'] = 11
    phan_df['ORFcaller'] = 'phanotate'
else:
    phan_df = pd.DataFrame()

if os.path.isfile(prod_gcodes) and os.path.getsize(prod_gcodes) > 0:
    prod_df = pd.read_csv(prod_gcodes, header = None, sep = "\t", names = ['contig', 'gcode'])
    prod_df['ORFcaller'] = 'prodigal'
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


## merge all files
merge_df = pd.merge(main_annot_df, name_df, on = "contig", how = "left")
merge_df = pd.merge(merge_df, tax_df, on = ["contig", "chunk_name"], how = "left")
merge_df = pd.merge(merge_df, gcode_df, on = "contig", how = "left")

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
        #gcode = re.search(r'\[gcode=(.*?)\]', fields[1]).group(1)
        desc_list.append([contig, chunkq, organism])
    except:
        print("except")

desc_df = pd.DataFrame(desc_list, columns=["contig", "chunk_name", "organism"])

org_info_df = pd.merge(merge_df, desc_df, on = ["contig", "chunk_name"], how = "left")

org_info_df['chunk_length'] = np.where(org_info_df['chunk_length'].isnull, 
                                       org_info_df['contig_length'], 
                                       org_info_df['chunk_length'])

## make summary of virus seqs
grouped_df = org_info_df.groupby(['contig', 'chunk_length', 'dtr_seq', 'chunk_name', 
                     'input_name', 'taxon', 'taxonomy_hierarchy', 'taxon_level',
                     'avg_hallmark_AAI_to_ref', 'organism', 'gcode', 'ORFcaller'], dropna = False)

summary_list = []
for name, group in grouped_df:
    if "Chunk" in str(name[3]):
        outname = "@".join([name[0], name[3]])
    else:
        outname = name[0]
    gene_count = group['gene_name'].nunique()
    vir_hall_count = group.query("Evidence_source == 'hallmark_hmm'")['gene_name'].nunique()
    vir_hall_list = '|'.join(
        list(group.query("Evidence_source == 'hallmark_hmm'")['evidence_description'])
        ).replace("-", " ")
    rep_hall_count = group.query("Evidence_source == 'rep_hall_hmm'")['gene_name'].nunique()
    rep_hall_list = '|'.join(
        list(group.query("Evidence_source == 'rep_hall_hmm'")['evidence_description'])
        ).replace("-", " ")
    if all(c in "ATCG" for c in str(name[2])):
        end_type = "DTR"
    else:
        end_type = "None"
        
    if gene_count >= 1:
        summary_list.append([outname, name[4], name[9], name[1], end_type, gene_count, vir_hall_count, rep_hall_count, 
                             vir_hall_list, rep_hall_list, name[6], name[11]])

summary_df = pd.DataFrame(summary_list, columns=['contig', 'input_name', 'organism', 'virus_seq_length', 
                                                 'end_feature', 'gene_count', 'virion_hallmark_count', 'rep_hallmark_count',
                                                 'virion_hallmark_genes', 'rep_hallmark_genes', 'taxonomy_hierarchy', 'ORF_caller'])

parentpath = Path(sequin_dir).parents[0]

summary_out = os.path.join(parentpath, f"{run_title}_virus_summary.tsv")

summary_df.to_csv(summary_out, sep = "\t", index = False)

## pruning summary
# where chunk names are not null:
prune_sum_df = org_info_df[['contig', 'contig_length', 'chunk_length', 'chunk_name', 'chunk_start', 'chunk_stop']]\
    .query("chunk_name == chunk_name").drop_duplicates()

if not prune_sum_df.empty:
    prune_out = os.path.join(parentpath, f"{run_title}_prune_summary.tsv")

    prune_sum_df.to_csv(prune_out, sep = "\t", index = False)