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
# main run directory
main_dir = sys.argv[4]
# run title
run_title = sys.argv[5]
# prodigal gcodes file
prod_gcodes = sys.argv[6]
# phanotate file
phan_file = sys.argv[7]
#ORFcaller arg
caller_arg = sys.argv[8]
# contig-to-organism table
org_table = sys.argv[9]
# samtools coverage table
samcov_table = sys.argv[10]

try:
    main_annot_df = pd.read_csv(gene_to_contig_table, sep = "\t")

    name_df = pd.read_csv(name_table, sep = "\t", names=['contig', 'input_name'])

    tax_df = pd.read_csv(tax_table, sep = "\t")

except:
    print(f"{os.path.basename(__file__)}: couldn't load files for summary")
    sys.exit()

# phanotate and prodigal
if os.path.isfile(phan_file) and os.path.getsize(phan_file) > 0:
    phan_df = pd.read_csv(phan_file, header = None, names = ['contig'])
    phan_df['gcode'] = 11
    phan_df['ORFcaller'] = 'phanotate'
else:
    phan_df = pd.DataFrame()

if os.path.isfile(prod_gcodes) and os.path.getsize(prod_gcodes) > 0:
    prod_df = pd.read_csv(prod_gcodes, header = None, sep = "\t", names = ['contig', 'gcode'])
    if caller_arg == 'prodigal':
        prod_df['ORFcaller'] = 'prodigal'
    else:
        prod_df['ORFcaller'] = 'prodigal-gv'
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
    print(f"{os.path.basename(__file__)}: couldn't load gcode table")

## samtools coverage table for coverage info
try:
    sam_df = pd.read_csv(samcov_table, sep = "\t", header = 0)[['#rname', 'coverage']]

    if (sam_df['#rname'].str.contains('@')).any():
        sam_df[['contig', 'chunk_name']] = sam_df['#rname'].str.split('@', n=1, expand=True)

    else:
        sam_df['contig'] = sam_df['#rname']
        sam_df['chunk_name'] = "NaN"

    sam_df = sam_df[['contig', 'chunk_name', 'coverage']]

    sam_df['chunk_name'] = sam_df['chunk_name'].infer_objects(copy=False).fillna("NaN")
except Exception as e:
    sam_df = pd.DataFrame()



## merge all files
merge_df = pd.merge(main_annot_df, name_df, on = "contig", how = "left")

tax_df['chunk_name'] = tax_df['chunk_name'].infer_objects(copy=False).fillna("NaN")
merge_df['chunk_name'] = merge_df['chunk_name'].infer_objects(copy=False).fillna("NaN")

merge_df = pd.merge(merge_df, tax_df, on = ["contig", "chunk_name"], how = "left")
merge_df = pd.merge(merge_df, gcode_df, on = "contig", how = "left")

if not sam_df.empty:
    merge_df = pd.merge(merge_df, sam_df, on = ["contig", "chunk_name"], how = "left")
else:
    merge_df['coverage'] = "NaN"

merge_df['taxon'] = merge_df['taxon'].infer_objects(copy=False).fillna("unclassified virus")


## get descriptions from fastas

try:
    desc_df = pd.read_csv(org_table, sep = "\t")
except:
    desc_df = pd.DataFrame()
    print(f"{os.path.basename(__file__)}: couldn't load contig-to-organism table")



## ensure merge gets same data type
desc_df['chunk_name'] = desc_df['chunk_name'].infer_objects(copy=False).fillna("NaN")
merge_df['chunk_name'] = merge_df['chunk_name'].infer_objects(copy=False).fillna("NaN")


org_info_df = pd.merge(merge_df, desc_df, on = ["contig", "chunk_name"], how = "left")

org_info_df['chunk_length'] = np.where(org_info_df['chunk_length'].isnull(), 
                                       org_info_df['contig_length'], 
                                       org_info_df['chunk_length'])

## make summary of virus seqs
grouped_df = org_info_df.groupby(['contig', 'chunk_length', 'chunk_name', 
                     'input_name', 'taxon', 'taxonomy_hierarchy', 'taxon_level',
                     'avg_hallmark_AAI_to_ref', 'organism', 'gcode', 'ORFcaller', 'coverage'], dropna = False)

summary_list = []
for name, group in grouped_df:
    if "Chunk" in str(name[2]):
        outname = "@".join([name[0], name[2]])
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
    rdrp_hall_count = group.query("Evidence_source == 'rdrp_hall_hmm'")['gene_name'].nunique()
    rdrp_hall_list = '|'.join(
        list(group.query("Evidence_source == 'rdrp_hall_hmm'")['evidence_description'])
        ).replace("-", " ")
    if not group.dtr_seq.replace('', np.nan).isna().all():
        
        dtr_seqf = group['dtr_seq'].mode()[0]
    else:
        dtr_seqf = "None"
    
    if all(c in "ATCG" for c in dtr_seqf):
        end_type = "DTR"
    else:
        end_type = "None"
        
    if gene_count >= 1:
        summary_list.append([outname, name[3], name[8], name[1], end_type, gene_count, vir_hall_count, rep_hall_count, 
                             rdrp_hall_count, vir_hall_list, rep_hall_list, rdrp_hall_list, name[5], name[10], 
                             name[9], name[11]])

summary_df = pd.DataFrame(summary_list, columns=['contig', 'input_name', 'organism', 'virus_seq_length', 
                                                 'end_feature', 'gene_count', 'virion_hallmark_count', 'rep_hallmark_count',
                                                 'RDRP_hallmark_count', 'virion_hallmark_genes', 'rep_hallmark_genes', 
                                                 'RDRP_hallmark_genes', 'taxonomy_hierarchy', 'ORF_caller',
                                                 'gcode', 'avg_read_depth'])


summary_df['virus_seq_length'] = summary_df['virus_seq_length'].infer_objects(copy=False).fillna(0)

summary_df = summary_df.astype(dtype= {"virus_seq_length": "int64"})

summary_out = os.path.join(main_dir, f"{run_title}_virus_summary.tsv")

summary_df.to_csv(summary_out, sep = "\t", index = False)

## pruning summary
# where chunk names are not null:
prune_sum_df = org_info_df[['contig', 'contig_length', 'chunk_length', 'chunk_name', 'chunk_start', 'chunk_stop']]\
    .query("chunk_name == chunk_name").drop_duplicates()

prune_sum_df = prune_sum_df[prune_sum_df['contig_length'].notna()]

prune_sum_df = prune_sum_df.astype(dtype= {"contig_length": "int64", "chunk_length": "int64", 
                                           "chunk_start": "int64", "chunk_stop": "int64"})

if not prune_sum_df.empty:
    prune_out = os.path.join(main_dir, f"{run_title}_prune_summary.tsv")

    prune_sum_df.to_csv(prune_out, sep = "\t", index = False)