#!/usr/bin/env python

import os
import sys
import pandas as pd
import math
import re
import numpy as np
from pathlib import Path

gene_contig_file = sys.argv[1]

tRNA_table = sys.argv[2]

phrogs_table = sys.argv[3]

out_dir = sys.argv[4]

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

# load gene/contig table
gene_contig_df = pd.read_csv(gene_contig_file, sep = "\t")


# load and format tRNA table from tRNAscan-SE
# check if there are any tRNAs predicted
if  os.path.isfile(tRNA_table) and os.path.getsize(tRNA_table) > 0:
    tRNA_df = pd.read_csv(tRNA_table, index_col=False, sep = "\t", 
                          names = ['con_chunk', 'tRNA_num', 'gene_c1', 'gene_c2', 
                                   'evidence_description', 'tRNA_codon', 'other1', 
                                   'other2', 'tRNA_score'])

    tRNA_df['at_pos'] = tRNA_df['con_chunk'].str.find("@")

    tRNA_df['contig'] = tRNA_df.apply(
        lambda x: x["con_chunk"][0:x["at_pos"]], axis = 1)

    tRNA_df['chunk_name'] = tRNA_df.apply(
        lambda x: x["con_chunk"][x["at_pos"]+1:-1], axis = 1)

    tRNA_df['gene_start'] = np.where(tRNA_df['gene_c1'] < tRNA_df['gene_c2'], tRNA_df['gene_c1'], tRNA_df['gene_c2'])

    tRNA_df['gene_stop'] = np.where(tRNA_df['gene_c1'] > tRNA_df['gene_c2'], tRNA_df['gene_c1'], tRNA_df['gene_c2'])

    tRNA_df['gene_orient'] = np.where(tRNA_df['gene_c1'] < tRNA_df['gene_c2'], "+", "-")

    tRNA_df['evidence_acession'] = "tRNAscan-SE score: " + tRNA_df['tRNA_score'].astype(str)

    tRNA_df['Evidence_source'] = "tRNAscan-SE"

    tRNA_df['gene_name'] = "tRNA-" + tRNA_df['evidence_description'].astype(str)

    tRNA_df = tRNA_df[['contig', 'chunk_name', 'gene_start', 'gene_stop', 'gene_name', 
                    'gene_orient', 'evidence_description', 'evidence_acession', 'Evidence_source']]

    # concat tables
    more_feature_df = pd.concat([gene_contig_df, tRNA_df], ignore_index=True)
else:
    more_feature_df = gene_contig_df


## check for phrogs table and parse
if os.path.isfile(phrogs_table) and os.path.getsize(phrogs_table) > 0:
    phrogs_pyh_df = pd.read_csv(phrogs_table, sep = "\t")[['ORFquery', 'target']]

    phrogs_pyh_df["gene_name"] = phrogs_pyh_df["ORFquery"]

    phrogs_pyh_df["slash_pos"] = phrogs_pyh_df["target"].str.find("/")
    phrogs_pyh_df["fdash_pos"] = phrogs_pyh_df["target"].str.find("-")


    phrogs_pyh_df["evidence_acession"] = phrogs_pyh_df.apply(
        lambda x: x["target"][x["slash_pos"]+1:x["fdash_pos"]], axis = 1)

    phrogs_pyh_df["evidence_description"] = phrogs_pyh_df.apply(lambda x: x["target"][x["fdash_pos"]+1:], axis = 1)

    phrogs_pyh_df = phrogs_pyh_df[['gene_name', 'evidence_acession', 'evidence_description']]

    phrogs_pyh_df['Evidence_source'] = 'phrogs_hmm'

else:
    print("nope")
    phrogs_pyh_df = pd.DataFrame()


## combine gene+contig table and phrogs table, replacing gene annotations
if not phrogs_pyh_df.empty:
    merged_df = pd.merge(more_feature_df, phrogs_pyh_df, on = "gene_name", how = "left",
                         suffixes = ("", "_y"))

    merged_df['evidence_acession'] = np.where(merged_df['evidence_acession_y']
                               .notnull(), merged_df['evidence_acession_y'], 
                               merged_df['evidence_acession'])

    merged_df['evidence_description'] = np.where(merged_df['evidence_description_y']
                               .notnull(), merged_df['evidence_description_y'], 
                               merged_df['evidence_description'])
    
    merged_df['Evidence_source'] = np.where(merged_df['Evidence_source_y']
                               .notnull(), merged_df['Evidence_source_y'], 
                               merged_df['Evidence_source'])
    
    merged_df = merged_df.drop(['evidence_acession_y', 'evidence_description_y', 
                                'Evidence_source_y'], axis = 1)
    
else:
    merged_df = gene_contig_df


#merged_df = merged_df.astype({'gene_start': 'int32', 'gene_stop': 'int32', 'gene_stop': 'int32', 
#                              'contig_length': 'int32', 'chunk_length': 'int32', 
#                              'chunk_start': 'int32', 'chunk_stop': 'int32'}).dtypes

## save genes to contigs annotation file to main run directory
parentpath = Path(out_dir).parents[0]

final_annotation_out = os.path.join(parentpath, "final_genes_to_contigs_annotation_summary.tsv")

merged_df.to_csv(final_annotation_out, sep = "\t", index = False)


## need a chunk value for all viruses to properly group
merged_df['chunk_name'] = merged_df['chunk_name'].fillna("nochunk")

chunk_grouped_df = merged_df.groupby(['contig', 'chunk_name'], dropna = False)


## coords in the right order for tbl files
def tbl_first_second(gstart, gstop, gorient):
    if gorient == "+":
        return gstart, gstop, gorient
    else:
        return gstop, gstart, gorient

#### loop each virus 
for name, seq_group in chunk_grouped_df:

    # header line for tbl
    # first for chunked/pruned
    if not name[1] == "nochunk":
        tbl_output_file = os.path.join(out_dir, name[0] + "@" + name[1] + ".tbl")

        print(f">Feature {name[0]}@{name[1]} Table1", file = open(tbl_output_file, "a"))
        
    # here for not chunked/pruned
    else:
        tbl_output_file = os.path.join(out_dir, name[0] + ".tbl")

        print(f">Feature {name[0]} Table1", file = open(tbl_output_file, "a"))

    # now each row is a feature that needs to be parsed and printed correctly
    for index, row in seq_group.iterrows():
        trna_number = 1

        first_c = tbl_first_second(row['gene_start'], row['gene_stop'], row['gene_orient'])[0]
        second_c = tbl_first_second(row['gene_start'], row['gene_stop'], row['gene_orient'])[1]

        if row['Evidence_source'] in {'hallmark_hmm', 'common_virus_hmm', 'phrogs_hmm'}: #my hmms
            typeq = "CDS"
            tagstr = ("protein_id" + "\tlcl|" + row['gene_name'])
            productstr = re.sub("-", " ", row['evidence_description'])
            inferencestr = ("inference\tprotein motif " + str(row['evidence_acession']))
        
        elif row['Evidence_source'] == "mmseqs_cdd":  #mmseqs_cdd
            typeq = "CDS"            
            tagstr = ("protein_id" + "\tlcl|" + row['gene_name'])
            productstr = re.sub("\..*", "", row['evidence_description'])
            inferencestr = ("inference\tprotein motif CDD:" + str(row['evidence_acession']))

        elif row['Evidence_source'] == "tRNAscan-SE": ## tRNAs
            typeq = "gene"
            tagstr = ("gene\t" + str(name[0]) + "tRNA" + str(trna_number))
            productstr = row['gene_name']
            inferencestr = ("inference\t" + row['evidence_acession'])
            trna_number =+ 1

        elif pd.isnull(row['Evidence_source']): #hypos
            typeq = "CDS"
            tagstr = ("protein_id" + "\tlcl|" + row['gene_name'])
            productstr = row['evidence_description']
            inferencestr = ("note\tno search hits")
        else:
            raise Exception("this shouldn't happen")
            typeq = "help"
            tagstr = "help"
            productstr = "help"
            inferencestr = "help"        
        
        print(f"{first_c}\t{second_c}\t{typeq}\n\t\t\t{tagstr}\n\t\t\tproduct\t{productstr}\n\t\t\t{inferencestr}", 
              file = open(tbl_output_file, "a"))

        