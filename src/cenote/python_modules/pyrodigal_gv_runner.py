#!/usr/bin/env python

import subprocess
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pyrodigal
import pyrodigal_gv
import multiprocessing.pool
import math
import re
import time

input_fasta = sys.argv[1]

out_dir = sys.argv[2]

CPUcount = sys.argv[3]

orf_caller = sys.argv[4]

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

starttime = time.time()

records = SeqIO.parse(input_fasta, "fasta")

if orf_caller == "prodigal":
    orf_finder = pyrodigal.GeneFinder(meta = True)
else:
    orf_finder = pyrodigal_gv.ViralGeneFinder(meta = True)

def _find_genes(record):
    genes = orf_finder.find_genes(str(record.seq))
    return (record.id, genes)

with multiprocessing.pool.ThreadPool(int(CPUcount)) as pool:
    with open(os.path.join(out_dir, "pyrodigal_gv_AAs.prod.faa"), "w") as dst:
        with open(os.path.join(out_dir, "pyrodigal_gv_AAs.prod.gff"), "w") as dgff:
            for record_id, genes in pool.imap(_find_genes, records):
                genes.write_translations(dst, sequence_id=record_id)
                genes.write_gff(dgff, sequence_id=record_id)

endtime = time.time()

time_taken = endtime - starttime

print("pyrodigal part took: " + "%.2f" % time_taken + " seconds")