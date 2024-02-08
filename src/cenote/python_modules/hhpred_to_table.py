#!/usr/bin/env python
import os
import sys
from io import StringIO
import pandas as pd

input_dir = sys.argv[1]

out_dir = sys.argv[2]


hhrout_list = []
for hhrout in os.listdir(input_dir):
    if hhrout.endswith('.hhr'):
        f = os.path.join(input_dir, hhrout)

        if os.path.isfile(f) and os.path.getsize(f) > 0:
            hhrout_list.append(f)

if not hhrout_list:
    print(f"{os.path.basename(__file__)}: no hhpred result files found in " + str(input_dir))
    sys.exit()

rows = []

for hhresult_file in hhrout_list:
    with open(hhresult_file, 'r') as hrh:
        for line in hrh:
            if line.startswith('Query         '):
                gene_name = line.strip('Query         ').split(' ')[0]
            
            #with this logic, only files with results (lines starting with >) append the df
            if line.startswith('>'):
                full_desc = line.rstrip('\n').strip('>')

                ## CDD models
                if line.startswith('>cd') or line.startswith('>sd'):
                    accession = full_desc.split(' ')[0]

                    try:
                        #annotation = full_desc.split('; ')[1].split('. ')[0]
                        annotation = full_desc.split(' ')[1].split(';')[0]
                    except:
                        annotation = 'hypothetical protein'

                ## PFAM models
                elif line.startswith('>PF'):
                    accession = full_desc.split(' ; ')[0]

                    try:
                        annotation = full_desc.split('; ')[2]
                    except:
                        annotation = 'hypothetical protein'

                ## PDB models
                else:
                    accession = full_desc.split(' ')[0]

                    try:
                        annotationpd = full_desc.split(' ')[1:]
                        annotationpd = ' '.join(annotationpd)
                        if ']' in annotationpd:
                            annotation = annotationpd.split(']')[1].split(';')[0]
                        else:
                            annotation = annotationpd.split(';')[0]
                    except:
                        annotation = 'hypothetical protein'

                full_stats = next(hrh, '').strip()
                probability = full_stats.split('  ')[0].split('=')[1]
                evalue = full_stats.split('  ')[1].split('=')[1]
                Aligned_cols = full_stats.split('  ')[3].split('=')[1]
                row = [gene_name, accession, annotation, probability, evalue, Aligned_cols]

                rows.append(row)

hhpred_df = pd.DataFrame(rows, columns=['gene_name', 'evidence_acession', 'evidence_description', 
                                        'evidence_prob', 'evalue', 'Aligned_cols'])

hhpred_df = hhpred_df.astype(dtype= {'evidence_prob': 'float', 'evalue': 'float', 'Aligned_cols': 'int32'})

hhpred_df = hhpred_df.query("evalue <= 1e-5")\
    .sort_values('evidence_prob', ascending = False).drop_duplicates('gene_name')

if not hhpred_df.empty:
    hhpred_output_file = os.path.join(out_dir, "hhpred_report_AAs.tsv")

    hhpred_df.to_csv(hhpred_output_file,
                            sep = "\t", index = False)