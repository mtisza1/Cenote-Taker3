#!/usr/bin/env python

import os
import re
import sys
import pandas as pd
import glob
import numpy as np

import itertools
from itertools import tee
import csv

import statistics
from statistics import mean
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math
import collections
import bisect

from prune_virus_coords1 import prune_chunks

repeat_table = sys.argv[1]

phan_tab_directory = sys.argv[2]

prod_tab_directory = sys.argv[3]

virion_pyhmmer_table = sys.argv[4]

comm_pyhmmer_table = sys.argv[5]

rep_pyhmmer_table = sys.argv[6]

rdrp_pyhmmer_table = sys.argv[7]

mmseqs_CDD_table = sys.argv[8]

viral_cdds_list = sys.argv[9]

out_dir = sys.argv[10]

hall_type = sys.argv[11]

PROPHAGE = sys.argv[12]

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

fig_out_dir = os.path.join(out_dir, "prune_figures")

if not os.path.isdir(fig_out_dir):
    os.makedirs(fig_out_dir)

## look for phanotate gene tables
try:

    phan_files = glob.glob(os.path.join(phan_tab_directory, "*.bed"))

    df_from_each_phan = (pd.read_csv(phan, sep = "\t", header = None,
                                    names = ["contig", "gene_bstart", "gene_bstop", "gene_name", 
                                            "gene_score", "gene_orient"])
                                    for phan in phan_files)
    phan_gene_df = pd.concat(df_from_each_phan, ignore_index=True)

    phan_gene_df['gene_start'] = (phan_gene_df['gene_bstart'] + 1).astype('int')
    phan_gene_df['gene_stop'] = phan_gene_df['gene_bstop'].astype('int')


    phan_gene_df = phan_gene_df.drop(['gene_score', 'gene_bstart', 'gene_bstop'], axis = 1)

except:
    print("no phanotate tables")
    phan_gene_df = pd.DataFrame()

## look for prodigal gene tables
try:
    prod_files = glob.glob(os.path.join(prod_tab_directory, "*.gff"))

    df_from_each_prod = (pd.read_csv(prod, sep = "\t", header = None, comment='#',
                                    names = ["contig", "gene_caller", "feature_type", "gene_start", 
                                            "gene_stop", "score", "gene_orient", "frame", "attribute"])
                                    for prod in prod_files)
    prod_gene_df = pd.concat(df_from_each_prod, ignore_index=True)

    prod_gene_df = prod_gene_df[["contig", "gene_start", "gene_stop", "attribute", "gene_orient"]]

    prod_gene_df["semi_pos"] = prod_gene_df["attribute"].str.find(";")
    prod_gene_df["gene_IDstr"] = prod_gene_df.apply(
        lambda x: x["attribute"][0:x["semi_pos"]], axis = 1)
    
    prod_gene_df["under_pos"] = prod_gene_df["gene_IDstr"].str.rfind("_")
    prod_gene_df["gene_num"] = prod_gene_df.apply(
        lambda x: x["gene_IDstr"][x["under_pos"]:len(x["gene_IDstr"])], axis = 1)
    
    prod_gene_df["gene_name"] = prod_gene_df["contig"] + prod_gene_df["gene_num"]

    prod_gene_df = prod_gene_df[["contig", "gene_start", "gene_stop", "gene_name", "gene_orient"]]

except:
    print("no prodigal tables")
    prod_gene_df = pd.DataFrame()

## combine gene tables
if not phan_gene_df.empty and not prod_gene_df.empty:
    print("both")
    both_list = [phan_gene_df, prod_gene_df]
    just_gene_df = pd.concat(both_list, ignore_index=True)
elif not phan_gene_df.empty:
    print("phan")
    just_gene_df = phan_gene_df
elif not prod_gene_df.empty:
    print("prod")
    just_gene_df = prod_gene_df

## get table with lengths and repeats for each hallmark contig
try:
    length_df = pd.read_csv(repeat_table, sep = "\t")[['contig', 'out_length_contig', 'dtr_seq']]

    length_df = length_df.rename(columns={"out_length_contig": "contig_length"})
except:
    print("nope")
    exit

## combine gene and contig tables
try:
    basal_df = just_gene_df.merge(length_df, on = "contig", how = "left")
except:
    print("nope")
    exit

## load and parse table for first pyhmmer search (hallmarks)\

def parse_pyhmmer_table(tab_file, evidence, categ):
    try:
        ppyh_df = pd.read_csv(tab_file, sep = "\t")[['ORFquery', 'target']]


        ppyh_df["gene_name"] = ppyh_df["ORFquery"]

        ppyh_df["slash_pos"] = ppyh_df["target"].str.find("/")
        ppyh_df["fdash_pos"] = ppyh_df["target"].str.find("-")


        ppyh_df["evidence_acession"] = ppyh_df.apply(
            lambda x: x["target"][x["slash_pos"]+1:x["fdash_pos"]], axis = 1)

        ppyh_df["evidence_description"] = ppyh_df.apply(lambda x: x["target"][x["fdash_pos"]+1:], 
                                                                    axis = 1)

        ppyh_df = ppyh_df[['gene_name', 'evidence_acession', 'evidence_description']]

        ppyh_df['Evidence_source'] = str(evidence)

        ppyh_df['vscore_category'] = str(categ)

    except:
        print("nope")
        ppyh_df = pd.DataFrame()
    return ppyh_df

virion_ppyh_df = parse_pyhmmer_table(virion_pyhmmer_table, "hallmark_hmm", "common_virus")
comm_pyh_df = parse_pyhmmer_table(comm_pyhmmer_table, "common_virus_hmm", "common_virus")
rep_pyh_df = parse_pyhmmer_table(rep_pyhmmer_table, "rep_hall_hmm", "common_virus")
rdrp_pyh_df = parse_pyhmmer_table(rdrp_pyhmmer_table, "rdrp_hall_hmm", "common_virus")


## load and parse table for mmseqs CDD search
try:
    cdd_df = pd.read_csv(mmseqs_CDD_table, sep = "\t")[['query', 'target', 'description']]

    cdd_df["gene_name"] = cdd_df["query"]

    cdd_df = cdd_df[['gene_name', 'target', 'description']]

    cdd_df = cdd_df.rename(columns={"target": "evidence_acession", "description" : "evidence_description"})

    cdd_df['Evidence_source'] = 'mmseqs_cdd'
except:
    print("no CDD")


## load file with list of additional common virus genes
try:
    virlist_df = pd.read_csv(viral_cdds_list, header = None, names = ["evidence_acession"])

    virlist_df['vscore_category'] = 'common_virus'
except:
    print("no virlist")

## combine mmseqs CDD search table and common virus gene list
try:
    comb_cdd_df = cdd_df.merge(virlist_df, on = "evidence_acession", how = "left")

    comb_cdd_df['vscore_category'] = np.where(comb_cdd_df['evidence_acession']
                                .str.contains(
                                    "PHA0",
                                    case = False), 'common_virus', comb_cdd_df['vscore_category'])

    comb_cdd_df['vscore_category'] = np.where(comb_cdd_df['vscore_category'].isna(), 'nonviral_gene',
                                            comb_cdd_df['vscore_category'])
except:
    print("no CDD mmseqs table")


## combine pyhmmer and mmseqs tables with contig/gene table for all gene annotations
gene_ann_list = []

if not virion_ppyh_df.empty:
    gene_ann_list.append(virion_ppyh_df)

if not comm_pyh_df.empty:
    gene_ann_list.append(comm_pyh_df)

if not rep_pyh_df.empty:
    gene_ann_list.append(rep_pyh_df)

if not rdrp_pyh_df.empty:
    gene_ann_list.append(rdrp_pyh_df)

if not comb_cdd_df.empty:
    gene_ann_list.append(comb_cdd_df)

try:
    gene_ann_df = pd.concat(gene_ann_list, ignore_index=True)

    contig_gene_df = basal_df.merge(gene_ann_df, on = "gene_name", how = "left")

    contig_gene_df['vscore_category'] = np.where(contig_gene_df['vscore_category'].isna(), 
                                                 'hypothetical_protein',
                                                 contig_gene_df['vscore_category'])
    
    contig_gene_df['evidence_description'] = np.where(contig_gene_df['evidence_description'].isna(), 
                                                      'hypothetical protein',
                                                      contig_gene_df['evidence_description'])

except:
    print("nope")

## save annotation table to file
contig_gene_outfile = os.path.join(out_dir, "contig_gene_annotation_summary.tsv")

contig_gene_df.to_csv(contig_gene_outfile, sep = "\t", index = False)

## save hallmark genes in bed format
hallmark_df = contig_gene_df.query("Evidence_source == 'hallmark_hmm'")
hallmark_df = hallmark_df[['contig', 'gene_start', 'gene_stop']]

hallmark_gene_outfile = os.path.join(out_dir, "contig_gene_annotation_summary.hallmarks.bed")

hallmark_df.to_csv(hallmark_gene_outfile, sep = "\t", index = False, header = False)


#try:
grouped_df = contig_gene_df.query("contig_length >= 10000")\
    .query("dtr_seq.isnull()").groupby('contig')



try:
    if PROPHAGE == "True":
        for name, group in grouped_df:

            prune_chunks(name, group, fig_out_dir, hall_type)
    else:
        print("--prune_prophage set to False, not pruning")
    
except:
    print("No non-DTR virus contigs >= 10,000 nt. So pruning will not happen")

    os.mkdir(os.path.join(out_dir, "prune_figures"))

