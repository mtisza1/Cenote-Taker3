#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import re
import sys
import pandas as pd

fasta_file = sys.argv[1]

orig_hallmark_tab = sys.argv[2]

length_circ = sys.argv[3]

length_lin = sys.argv[4]

hall_circ = sys.argv[5]

hall_lin = sys.argv[6]

out_dir = sys.argv[7]

wrap = sys.argv[8]

maxlength = sys.argv[9]


if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

output_allf = os.path.join(out_dir, "trimmed_TRs_hallmark_contigs.fasta")

if os.path.isfile(output_allf):
    os.remove(output_allf)


## The following functions are copied wholesale from CheckV as I find them to be fast and correct
## link https://bitbucket.org/berkeleylab/checkv/src/master/checkv/modules/complete_genomes.py

def fetch_dtr(fullseq, min_length=20):
    startseq = fullseq[0:min_length]
    # find index positions of all matches of startseq in fullseq
    # only keep matches occuring in 2nd half of string
    matches = [
        m.start() for m in re.finditer("(?={0})".format(re.escape(startseq)), fullseq)
    ]
    matches = [_ for _ in matches if _ >= len(fullseq) / 2]
    for matchpos in matches:
        # determine if the match extends to the contig end
        endseq = fullseq[matchpos:]
        if fullseq[0 : len(endseq)] == endseq:
            return endseq
    return ""

def reverse_complement(seq):
    if sys.version_info > (3, 0):
        trans = str.maketrans("ACTG", "TGAC")
    else:
        trans = string.maketrans("ACTG", "TGAC")
    return seq[::-1].translate(trans)

def fetch_itr(seq, min_len=20, max_len=1000):
    rev = reverse_complement(seq)
    # see if minimal substring occurs at end
    if seq[:min_len] == rev[:min_len]:
        # extend to maximum substring, up to <max_len>
        i = min_len + 1
        while seq[:i] == rev[:i] and i <= max_len:
            i += 1
        return seq[: i - 1]
    # no match
    else:
        return ""




terminal_r_list = []
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    dtr_seq = fetch_dtr(str(seq_record.seq))

    if not dtr_seq or len(seq_record) > float(maxlength):
        dtr_seq = "NA"

    itr_seq = fetch_itr(str(seq_record.seq))

    if not itr_seq:
        itr_seq = "NA"

    if not dtr_seq == "NA" and wrap.lower() == "true":
        print(f">{seq_record.id}", file = open(output_allf, "a"))
        print(seq_record.seq[:-len(dtr_seq)], file = open(output_allf, "a"))

        terminal_r_list.append([seq_record.id, len(seq_record.seq), len(seq_record.seq[:-len(dtr_seq)]), dtr_seq, itr_seq])
    else:
        print(f">{seq_record.id}", file = open(output_allf, "a"))
        print(seq_record.seq, file = open(output_allf, "a"))

        terminal_r_list.append([seq_record.id, len(seq_record.seq), len(seq_record.seq), dtr_seq, itr_seq])

terminal_df = pd.DataFrame(terminal_r_list, columns=["contig", "in_length_contig", "out_length_contig", 
                                                     "dtr_seq", "itr_seq"])

all_table = os.path.join(out_dir, "hallmark_contigs_terminal_repeat_summary.tsv")

terminal_df.to_csv(all_table,
                        sep = "\t", index = False)


hall_dt = pd.read_csv(orig_hallmark_tab, sep = "\t", header = 0)

full_dt = pd.merge(hall_dt, terminal_df, on = 'contig')


keep_dt = full_dt[(full_dt['dtr_seq'].notna() 
        & (full_dt['total_hit_count'] >= int(hall_circ) ) 
        & (full_dt['in_length_contig'] >= int(length_circ) )) 
        | ((full_dt['total_hit_count'] >= int(hall_lin) )
        & (full_dt['in_length_contig'] >= int(length_lin) ))][["contig", "in_length_contig", "out_length_contig", 
                                                            "dtr_seq", "itr_seq"]]
thresh_table = os.path.join(out_dir, "threshold_contigs_terminal_repeat_summary.tsv")

keep_dt.to_csv(thresh_table, sep = "\t", index = False)

threshold_contigs = keep_dt['contig']

thresh_list = os.path.join(out_dir, "contigs_over_threshold.txt")
threshold_contigs.to_csv(thresh_list, sep = "\t", index = False, header = False)


