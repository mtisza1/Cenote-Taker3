#!/usr/bin/env python
import pandas as pd
import sys
import os

mmseqs_search_table = sys.argv[1]

cdd_annotation_table = sys.argv[2]

out_dir = sys.argv[3]

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

mmseqs2_CDD_search_df = pd.read_csv(mmseqs_search_table,
                                    sep = "\t", header = None,
                                    names = ["query", "target", "sequence_identity", "align_length",
                                             "evalue", "bitscore"])\
    .sort_values('evalue').drop_duplicates('query').query("evalue <= 1e-3")

cdd_annotation_df = pd.read_csv(cdd_annotation_table,sep = "\t", header = None,
                                names = ["cdd_num", "cdd_accession", "shortname", "description",
                                         "other_num"])

mmseqs2_CDD_annotation_df = mmseqs2_CDD_search_df.merge(cdd_annotation_df, 
                                                        left_on = "target", 
                                                        right_on = "cdd_accession",
                                                        how = "left")\
                                                        .sort_values('query')

mmseqs_cdd_sum_file = os.path.join(out_dir, "summary_no2_AAs_vs_CDD.besthit.tsv")

mmseqs2_CDD_annotation_df.to_csv(mmseqs_cdd_sum_file, sep = "\t", index = False)