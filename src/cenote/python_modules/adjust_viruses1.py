#!/usr/bin/env python

import os
import sys
import pandas as pd
import glob
import numpy as np
import itertools

virus_coord_dir = sys.argv[1]

gene_annotation_file = sys.argv[2]

out_dir = sys.argv[3]

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

try:

    coord_files = glob.glob(os.path.join(virus_coord_dir, "*.viruses.tsv"))

    df_from_each_coord = (pd.read_csv(coord, sep = "\t", header = None,
                                    names = ["contig", "chunk_start", "chunk_stop", 
                                             "chunk_name", "hallmark_count"])
                                    for coord in coord_files)
    chunk_coord_df = pd.concat(df_from_each_coord, ignore_index=True)

except:
    print("no chunk coords")

adjusted_gene_df = pd.DataFrame()

gene_contig_df = pd.read_csv(gene_annotation_file, sep = "\t")

grouped_df = gene_contig_df.groupby('contig')

for name, group in grouped_df:
    contig_length1 = group['contig_length'].agg(pd.Series.mode)

    dtr_status = group['dtr_seq'].agg(pd.Series.mode)

    if len(dtr_status) == 0 and int(contig_length1.iloc[0]) >=10000:
        ## lift genes from each chunk

        ## @name is how to call the name variable in query
        thisgroup_chunk_df = chunk_coord_df.query("contig == @name")

        for index, row in thisgroup_chunk_df.iterrows():
            chunk_start = row.loc['chunk_start'] - 1
            if chunk_start < 0:
                chunk_start = 0
            chunk_stop = row.loc['chunk_stop']
            chunk_name = row.loc['chunk_name']
            
            thischunk_gene_df = group.copy().query("gene_start >= @chunk_start & gene_stop <= @chunk_stop")

            thischunk_gene_df['gene_start'] = thischunk_gene_df['gene_start'] - chunk_start
            thischunk_gene_df['gene_stop'] = thischunk_gene_df['gene_stop'] - chunk_start
            thischunk_gene_df['chunk_name'] = chunk_name
            thischunk_gene_df['chunk_length'] = chunk_stop - chunk_start
            thischunk_gene_df['chunk_start'] = chunk_start
            thischunk_gene_df['chunk_stop'] = chunk_stop

            adjusted_gene_df = pd.concat([adjusted_gene_df, thischunk_gene_df], ignore_index=True)

    else:
        #dtrs and seqs under 10kb
        ## for these just return the gene coordinates as is
        ## make chunk column

        group['chunk_name'] = "NA"
        group['chunk_length'] = "NA"
        group['chunk_start'] = 0
        group['chunk_stop'] = group['contig_length']

        adjusted_gene_df = pd.concat([adjusted_gene_df, group], ignore_index=True)

adjusted_gene_file = os.path.join(out_dir, "contig_gene_annotation_summary.pruned.tsv")

adjusted_gene_df.to_csv(adjusted_gene_file, sep = "\t", index = False)

bed_df = adjusted_gene_df[['contig', 'chunk_start', 'chunk_stop', 'chunk_name']].drop_duplicates()

bed_df['new_chunk'] = np.where(bed_df['chunk_name'] == "NA", bed_df['contig'], bed_df['contig'] + "@" + bed_df['chunk_name'])

bed_df = bed_df.drop('chunk_name', axis = 1)

bed_file = os.path.join(out_dir, "prune_coords.bed")

bed_df.to_csv(bed_file, sep = "\t", index = False, header = False)
