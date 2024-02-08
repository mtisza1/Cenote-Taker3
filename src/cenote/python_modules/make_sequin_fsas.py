#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pandas as pd
import re
import os
import random
import string
import sys

final_contig_file = sys.argv[1]

final_tax_file = sys.argv[2]

repeat_file = sys.argv[3]

temp_dir = sys.argv[4]

out_dir = sys.argv[5]

phanotate_list_file = sys.argv[6]

prodigal_list_file = sys.argv[7]

ISO_SOURCE = sys.argv[8]

COLLECT_DATE = sys.argv[9]

META_TYPE = sys.argv[10]

SRR = sys.argv[11]

SRX = sys.argv[12]

BIOSAMP = sys.argv[13]

PRJ = sys.argv[14]

MOL_TYPE = sys.argv[15]

DATA_SOURCE = sys.argv[16]

GENBANK = sys.argv[17]

# make out dir
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

# load dataframes

# tax
tax_call_df = pd.read_csv(final_tax_file, sep = "\t")

# phanotate and prodigal
if os.path.isfile(phanotate_list_file) and os.path.getsize(phanotate_list_file) > 0:
    phan_df = pd.read_csv(phanotate_list_file, header = None, names = ['contig'])
    phan_df['gcode'] = 11
else:
    phan_df = pd.DataFrame()


if os.path.isfile(prodigal_list_file) and os.path.getsize(prodigal_list_file) > 0:
    prod_df = pd.read_csv(prodigal_list_file, header = None, sep = "\t", names = ['contig', 'gcode'])
else:
    prod_df = pd.DataFrame()


## combine phanotate and prodigal table
gcode_list = []

for df in phan_df, prod_df:
    if not df.empty:
        gcode_list.append(df)


try:
    gcode_df = pd.concat(gcode_list, ignore_index=True)
except:
    print(f"{os.path.basename(__file__)}: genetic code table not found")
    gcode_df = pd.DataFrame()

# repeat and length
repeat_df = pd.read_csv(repeat_file, sep = "\t")

#### loop each virus

desc_list = []
for seq_record in SeqIO.parse(final_contig_file, "fasta"):

    #make a random alphanumeric ID 5 characters in length
    randID = ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))

    #chunked/pruned seqs need different handling
    if "@" in seq_record.id:

        nameq = seq_record.id.split("@")[0]
        chunkq = seq_record.id.split("@")[1]
        try:
            organism = tax_call_df.query("contig == @nameq & chunk_name == @chunkq")['taxon'].mode()[0]
        except:
            organism = "unclassified virus"

        try:
            lineage = tax_call_df.query("contig == @nameq & chunk_name == @chunkq")['taxonomy_hierarchy'].mode()[0]
        except:
            lineage = "no lineage"
        if not gcode_df.empty:
            gcode = gcode_df.query("contig == @nameq")['gcode'].mode()[0]
        else:
            gcode = 1

        topology = "linear"
    else:
        try:
            organism = tax_call_df.query("contig == @seq_record.id")['taxon'].mode()[0]
        except:
            organism = "unclassified virus"

        try:
            lineage = tax_call_df.query("contig == @seq_record.id")['taxonomy_hierarchy'].mode()[0]
        except:
            lineage = "no lineage"

        if not gcode_df.empty:
            gcode = gcode_df.query("contig == @seq_record.id")['gcode'].mode()[0]
        else:
            gcode = 1
            
        top_str = repeat_df.query("contig == @seq_record.id")['dtr_seq'].mode()

        if not top_str.empty:
            topology = "circular"
        else:
            topology = "linear"
    
    if GENBANK == 'True':
        # here's the whole header string
        header = (f">{str(seq_record.id)} [organism={organism} sp. ct{randID}] [gcode={gcode}] "
            f"[topology={topology}] [note: taxonomic lineage {lineage}] [isolation_source={ISO_SOURCE}] "
            f"[collection_date={COLLECT_DATE}] [metagenome_source={META_TYPE}] [SRA={SRR}] "
            f"[note=genome binned from sequencing reads available in {SRX}] [Biosample={BIOSAMP}] "
            f"[Bioproject={PRJ}] [moltype={MOL_TYPE}] [isolation_source={DATA_SOURCE}]")
        
        seq_output_file = os.path.join(out_dir, str(seq_record.id) + ".fsa")

        print(f"{header}\n{seq_record.seq}", file = open(seq_output_file, "a"))

    ## add outtable with contig, chunk, organism name (full)
    try:
        if "@" in seq_record.id:
            contig = seq_record.id.split("@")[0]
            chunkq = seq_record.id.split("@")[1]
        else:
            contig = seq_record.id
            chunkq = None

        fullorg = f'{organism} sp. ct{randID}'

        desc_list.append([contig, chunkq, fullorg])
    except:
        print(f"{os.path.basename(__file__)}: seq record info parse failed.")


desc_df = pd.DataFrame(desc_list, columns=["contig", "chunk_name", "organism"])

if not desc_df.empty:
    org_output_file = os.path.join(temp_dir, "contig_to_organism.tsv")

    desc_df.to_csv(org_output_file, sep = "\t", index = False)
