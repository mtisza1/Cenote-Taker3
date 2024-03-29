#!/usr/bin/env python
import os
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

in_fna = sys.argv[1]
min_length = sys.argv[2]
vir_sum = sys.argv[3]
prune_sum = sys.argv[4]
genes_sum = sys.argv[5]
ann_mode = sys.argv[6]
prune = sys.argv[7]

###

num_seqs = sum( 1 for _ in enumerate(SeqIO.parse(in_fna, "fasta")))

gene_df = pd.read_csv(genes_sum, sep = '\t', header = 0)

gene_df['chunk_name'] = gene_df['chunk_name'].infer_objects(copy=False).fillna("NaN")

gene_df['contig_chunk'] = gene_df['contig'] + "@" + gene_df['chunk_name']

num_virus_contigs = len(gene_df['contig_chunk'].drop_duplicates())

all_genes = len(gene_df['gene_name'].drop_duplicates())
nonhypo_genes = len(gene_df.query('evidence_description != "hypothetical protein"')['gene_name'].drop_duplicates())

try:
      prune_df = pd.read_csv(prune_sum, sep = '\t', header = 0)

      over_10kb = len(prune_df.query('contig_length >= 10_000')['contig'].drop_duplicates() )

      pruned_seqs = len(prune_df.query('contig_length >= 10_000 & chunk_length < contig_length')['contig'].drop_duplicates())
except:
      over_10kb = 0
      pruned_seqs = 0


lc = '\033[91m'
pc = '\033[35m'
endc = '\033[0m'
print('')
if ann_mode == 'True':
      print(f'>>> {pc}{num_seqs}{lc} contigs were annotated.{endc}')
      print(f'>>> {lc}In all, {pc}{nonhypo_genes / all_genes:.0%}{lc} of '
            f'virus genes were annotated with functional information.{endc}')
else:
      print(f'>>> {pc}{num_seqs}{lc} contigs over {pc}{min_length}{lc} nt were '
            f'searched and {pc}{num_virus_contigs}{lc} viruses were detected and annotated.{endc}')
      print(f'>>> {lc}In all, {pc}{nonhypo_genes / all_genes:.0%}{lc} of '
            f'virus genes were annotated with functional information.{endc}')
      if prune == 'True':
            print(f'>>> {pc}{over_10kb}{lc} contig(s) over 10kb went through pruning module '
                  f'and {pc}{pruned_seqs}{lc} were shortened by pruning.{endc}')

print('')