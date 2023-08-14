#!/usr/bin/env python

#!# this is purely chatgpt that I'm slapping in as a placeholder 230720

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import subprocess
import os
from pyhmmer import hmmscan as hmmscan
import pyrodigal

def run_prodigal(input_file):
    # Run Prodigal to predict ORFs
    #cmd = f"prodigal -i {input_file} -a {output_file}"
    #subprocess.run(cmd, shell=True)
    orf_finder = pyrodigal.OrfFinder(meta = True)
    genes = orf_finder.find_genes(input_file)

def parse_prodigal_output(output_file, fasta_file):
    # Parse Prodigal output to extract ORF sequences
    orf_records = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        orf_file = f"{output_file}.faa"
        seq_length = len(seq_record.seq)
        cmd = f"grep -A 1 '{seq_record.id}' {output_file} > {orf_file}"
        subprocess.run(cmd, shell=True)
        with open(orf_file) as f:
            for line in f:
                if line.startswith(">"):
                    # Extract ORF sequence and create a SeqRecord
                    orf_id = line.strip().lstrip('>')
                    orf_seq = next(f).strip()
                    if len(orf_seq) > 0:
                        orf_record = SeqRecord(Seq(orf_seq), id=orf_id, name=orf_id, description="")
                        orf_records.append(orf_record)
        os.remove(orf_file)
    return orf_records

def main():
    input_file = "input.fasta"   # Replace with the path to your input multi-record FASTA file
    output_file = "prodigal_output.gff"   # Replace with the desired output GFF file

    # Step 1: Run Prodigal to predict ORFs
    run_prodigal(input_file, output_file)

    # Step 2: Parse Prodigal output and create ORF records
    orf_records = parse_prodigal_output(output_file, input_file)

    # Step 3: Save the predicted ORFs to a FASTA file
    orf_fasta_file = "predicted_orfs.fasta"
    SeqIO.write(orf_records, orf_fasta_file, "fasta")
    print(f"Predicted ORFs saved to {orf_fasta_file}")

if __name__ == "__main__":
    main()
