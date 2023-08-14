#!/usr/bin/env python

import os
import sys
import pandas as pd
import math
import re
import numpy as np

mmseqs2_tax_table = sys.argv[1]

annotation_table = sys.argv[2]

out_dir =  sys.argv[3]

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

## load, parse, combine tables
tax_df = pd.read_csv(mmseqs2_tax_table, header = None, sep = "\t",
                     names = ["query","target","pident","alnlen","evalue","theader","taxlineage"])\
                     .sort_values('evalue').drop_duplicates('query').query("evalue <= 1e-3")

gene_contig_df = pd.read_csv(annotation_table, sep = "\t")[['contig', 'chunk_name', 'gene_start',
                                                            'gene_stop', 'gene_name']]

chunk_wise_df = gene_contig_df.merge(tax_df, left_on = "gene_name", right_on = "query", how = "inner")

## break out all the taxa of interest
chunk_wise_df['genus'] = np.where(chunk_wise_df['taxlineage'].str.contains(";g_"), 
                                     chunk_wise_df['taxlineage'].apply(lambda st: st[st.find(";g_")+3:].split(";")[0]), 
                                     "NA")

chunk_wise_df['family'] = np.where(chunk_wise_df['taxlineage'].str.contains(";f_"), 
                                     chunk_wise_df['taxlineage'].apply(lambda st: st[st.find(";f_")+3:].split(";")[0]), 
                                     "NA")

chunk_wise_df['order'] = np.where(chunk_wise_df['taxlineage'].str.contains(";o_"), 
                                     chunk_wise_df['taxlineage'].apply(lambda st: st[st.find(";o_")+3:].split(";")[0]), 
                                     "NA")

chunk_wise_df['taxclass'] = np.where(chunk_wise_df['taxlineage'].str.contains(";c_"), 
                                     chunk_wise_df['taxlineage'].apply(lambda st: st[st.find(";c_")+3:].split(";")[0]), 
                                     "NA")

chunk_wise_df['phylum'] = np.where(chunk_wise_df['taxlineage'].str.contains(";p_"), 
                                     chunk_wise_df['taxlineage'].apply(lambda st: st[st.find(";p_")+3:].split(";")[0]), 
                                     "NA")

## group by chunk
group_chunk_df = chunk_wise_df.groupby(['contig', 'chunk_name'], dropna = False)

## function to decide taxonomy
def taxon_decider(name, group, taxonomy_list):

    # get info for each taxonomy level
    group_count = len(group.index)

    # modal/most common genus name
    g_mode = group['genus'].agg(pd.Series.mode)

    g_mode_occur = len(group.query("genus == @g_mode[0]").index)

    g_mean_AAI = group.query("genus == @g_mode[0]")['pident'].mean()

    #
    f_mode = group['family'].agg(pd.Series.mode)

    f_mode_occur = len(group.query("family == @f_mode[0]").index)

    f_mean_AAI = group.query("family == @f_mode[0]")['pident'].mean()    
    #
    o_mode = group['order'].agg(pd.Series.mode)

    o_mode_occur = len(group.query("order == @o_mode[0]").index)

    o_mean_AAI = group.query("order == @o_mode[0]")['pident'].mean()  
    #
    c_mode = group['taxclass'].agg(pd.Series.mode)

    c_mode_occur = len(group.query("taxclass == @c_mode[0]").index)

    c_mean_AAI = group.query("taxclass == @c_mode[0]")['pident'].mean()  
    #
    p_mode = group['phylum'].agg(pd.Series.mode)

    p_mode_occur = len(group.query("phylum == @p_mode[0]").index)

    p_mean_AAI = group.query("phylum == @p_mode[0]")['pident'].mean()  

    ## calls genus if genus is not "NA", 75% of hallmarks agree, and mean AAI o alignments is >= 80%
    if g_mode[0] != "NA" and (g_mode_occur / group_count) >= 0.75 and g_mean_AAI >= 80:
        taxonomy = g_mode[0]
        level = "genus"
        mean_AAI = g_mean_AAI

    ## calls family if family is not "NA", 75% of hallmarks agree, and mean AAI o alignments is >= 35%
    elif f_mode[0] != "NA" and (f_mode_occur  / group_count) >= 0.75 and f_mean_AAI >= 35:
        taxonomy = f_mode[0]
        level = "family"
        mean_AAI = f_mean_AAI

    ## calls order if order is not "NA", 75% of hallmarks agree, and mean AAI o alignments is >= 20%
    elif o_mode[0] != "NA" and (o_mode_occur  / group_count) >= 0.75 and o_mean_AAI >= 20:
        taxonomy = o_mode[0]
        level = "order"
        mean_AAI = o_mean_AAI

    ## calls class if class is not "NA", 75% of hallmarks agree
    elif c_mode[0] != "NA" and (c_mode_occur  / group_count) >= 0.75:
        taxonomy = c_mode[0]
        level = "class"
        mean_AAI = c_mean_AAI

    ## calls phylum if phylum is not "NA", 75% of hallmarks agree
    elif p_mode[0] != "NA" and (p_mode_occur  / group_count) >= 0.75:
        taxonomy = p_mode[0]
        level = "phylum"
        mean_AAI = p_mean_AAI

    ## anything else is unclassified
    else:
        taxonomy = "unclassified virus"
        level = "NA"
        mean_AAI = "NA"

        
    try:
        lineage_to_tax_S = np.where(group['taxlineage'].str.contains(taxonomy), 
                                     group['taxlineage'].apply(lambda st: st[:st.find(taxonomy)+len(taxonomy)]), 
                                     taxonomy)
    
        lineage_to_tax = lineage_to_tax_S[0]
    except:
        lineage_to_tax = "NA"
    
    return taxonomy_list.append([name[0], name[1], taxonomy, lineage_to_tax, level, mean_AAI])


taxonomy_list_ct = []

for name_ct, group_ct in group_chunk_df:
    taxon_decider(name_ct, group_ct, taxonomy_list_ct)


taxonomy_call_df = pd.DataFrame(taxonomy_list_ct, columns = ['contig', 'chunk_name', 'taxon', 'taxonomy_hierarchy', 'taxon_level', 'avg_hallmark_AAI_to_ref'])

## save tax call table to file
tax_call_outfile = os.path.join(out_dir, "virus_taxonomy_summary.tsv")

taxonomy_call_df.to_csv(tax_call_outfile,
                        sep = "\t", index = False)