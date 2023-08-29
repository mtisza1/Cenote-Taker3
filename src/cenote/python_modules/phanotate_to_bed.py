#!/usr/bin/env python

import os
import sys
import pandas as pd
import numpy as np

f = sys.argv[1]

outf = sys.argv[2]

if os.path.isfile(f) and os.path.getsize(f) > 0:
    phan_df = pd.read_csv(f, header = None, index_col=False, sep = "\t", comment='#',
                          names = ["first_coord","second_coord","orient","contig","score"])
    
    phan_df = phan_df.astype({'first_coord': 'int32', 'second_coord': 'int32'})

    phan_df['start'] = np.where(phan_df['second_coord'] > phan_df['first_coord'], 
                                phan_df['first_coord'] - 1,
                                phan_df['second_coord'] - 1)
    
    phan_df['stop'] = np.where(phan_df['second_coord'] > phan_df['first_coord'], 
                                phan_df['second_coord'],
                                phan_df['first_coord'])
    
    phan_df['gene_number'] = phan_df.reset_index().index
    phan_df['gene_name'] = str(phan_df['contig']) + "_" + str(phan_df['gene_number'])

    bed_df = phan_df[['contig', 'start', 'stop', 'gene_name', 'score', 'orient']]

    bed_df = bed_df.astype({'start': 'int32', 'stop': 'int32'})
    
    bed_df.to_csv(outf, sep = "\t", index = False, header = False)
else:
    print("", file = open(outf, "a"))