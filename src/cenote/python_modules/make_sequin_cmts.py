#!/usr/bin/env python

import subprocess
import os
import sys
import pandas as pd


samtools_cov = sys.argv[1]

sequin_dir = sys.argv[2]

assembler = sys.argv[3]

seqtech = sys.argv[4]

# make list of fsa files
fsa_list = []
for fsa in os.listdir(sequin_dir):
    if fsa.endswith('.fsa'):
        f = os.path.join(sequin_dir, fsa)

        if os.path.isfile(f) and os.path.getsize(f) > 0:
            fsa_list.append(f)

if not fsa_list:
    print(f"{os.path.basename(__file__)}: no .fsa files found in " + str(sequin_dir))
    sys.exit()

# makes empty dataframe if samtools coverage table wasn't created
try:
    sam_df = pd.read_csv(samtools_cov, sep = "\t", header = 0)
except:
    sam_df = pd.DataFrame()


# matches coverage to contig by name
for fsa_file in fsa_list:
    with open(fsa_file) as fsaf:
        contig = fsaf.readline().strip('\n').strip('>').split(' ')[0]

        try:
            coverage: float = float(sam_df['coverage'][sam_df['#rname'] == contig].iloc[0,])

        except:
            coverage = 0
        
        out_cmt = os.path.splitext(fsa_file)[0]+'.cmt'
        
        try:
            # cmt file 
            print(f'StructuredCommentPrefix\t##Genome-Assembly-Data-START##\n'
                f'Assembly Method\t{assembler}\n'
                f'Genome Coverage\t{coverage:.2f}x\n'
                f'Sequencing Technology\t{seqtech}\n'
                f'Annotation Pipeline\tCenote-Taker 3\n'
                f'URL\thttps://github.com/mtisza1/Cenote-Taker3',
                file = open(out_cmt, "a"))
        except Exception as e:
            print(e)