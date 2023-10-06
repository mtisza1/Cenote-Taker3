#!/usr/bin/env python

import pandas as pd
import os
import sys

virion_count_dir = sys.argv[1]

virion_count_file = os.path.join(virion_count_dir, "contig_hit_count.tsv")

virion_tsv_file = os.path.join(virion_count_dir, "pyhmmer_report_AAs.tsv")


rep_count_dir = sys.argv[2]

rep_count_file = os.path.join(rep_count_dir, "contig_hit_count.tsv")

rep_tsv_file = os.path.join(rep_count_dir, "pyhmmer_report_AAs.tsv")


hallmark_count = sys.argv[3]

hallmark_count = int(hallmark_count)

hallmark_type = sys.argv[4]

out_dir = sys.argv[5]

names_file = sys.argv[6]

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
    
## load and parse hit count files
try:
    virion_dt = pd.read_csv(virion_count_file, sep = "\t", names=['contig', 'virion_hit_count'], skiprows = 1)
except:
    print("no virion hmm tsv")
    virion_dt  = pd.DataFrame(columns = ['contig', 'virion_hit_count'])

try:
    rep_dt = pd.read_csv(rep_count_file, sep = "\t", names=['contig', 'rep_hit_count'], skiprows = 1)
except:
    print("no rep hmm tsv")
    rep_dt  = pd.DataFrame(columns = ['contig', 'rep_hit_count'])

names_dt = pd.read_csv(names_file, sep = "\t", names=['contig', 'original_name'])[['contig']]

merge_dt = pd.merge(virion_dt, rep_dt, on = 'contig', how = 'outer')

merge_dt = pd.merge(merge_dt, names_dt, on = 'contig', how = 'outer')

merge_dt = merge_dt.fillna(0)

merge_dt['total_hit_count'] = merge_dt['virion_hit_count'] + merge_dt['rep_hit_count']

## depending on settings, call contigs with minimum hallmark genes
if hallmark_type == 'virion':
    contigs_w_min_hall = merge_dt.query("virion_hit_count >= @hallmark_count")['contig']
else:
    contigs_w_min_hall = merge_dt.query("total_hit_count >= @hallmark_count")['contig']

if not contigs_w_min_hall.empty:

    hallmark_contigs = len(contigs_w_min_hall.index)
    print(f"{hallmark_contigs} contigs with at least {hallmark_count} {hallmark_type} hallmark genes found")

    contig_keep_file = os.path.join(out_dir, "contigs_to_keep.txt")

    contigs_w_min_hall.to_csv(contig_keep_file, sep = "\t", index = False, header = False)
    
    halls_orig_contigs_file = os.path.join(out_dir, "hallmarks_per_orig_contigs.tsv")

    merge_dt.to_csv(halls_orig_contigs_file, sep = "\t", index = False)   
    
else:
    print(f"no contigs with at least {hallmark_count} {hallmark_type} hallmark genes found")


try:
    virion_hits_dt = pd.read_csv(virion_tsv_file, sep = "\t")
except:
    virion_hits_dt = pd.DataFrame()

try:
    rep_hits_dt = pd.read_csv(rep_tsv_file, sep = "\t")
except:
    rep_hits_dt = pd.DataFrame()

df_list = []
for hit_t in [virion_hits_dt, rep_hits_dt]:
    if not hit_t.empty:
        df_list.append(hit_t)

try:
    full_hits_dt = pd.concat(df_list, ignore_index=True)
except:
    print("couldn't make full hits list")

hm_contigs_hits_dt = merge_dt = pd.merge(contigs_w_min_hall, full_hits_dt, on = 'contig', how = 'left')

hm_contigs_hits_l = hm_contigs_hits_dt['ORFquery']

if not hm_contigs_hits_l.empty:

    hallmark_list_file = os.path.join(out_dir, "hallmarks_for_keepcontigs1.txt")

    hm_contigs_hits_l.to_csv(hallmark_list_file, sep = "\t", index = False, header = False)
else:
    print("no hallmark genes list")
