#!/usr/bin/env python

import sys
import os
import argparse
import re

#Author: Jessica Gomez-Garrido, CNAG.
#Contact email: jessica.gomez@cnag.eu
#Date:20250612

def parse_paf_and_get_longest_alignments(paf_file, match, sex_chroms):
    """
    Parses a PAF file and returns the longest alignment per target sequence.

    Args:
        paf_file (str): Path to the PAF file.

    Returns:
        dict: Dictionary where keys are target names and values are tuples:
              (query_name, alignment_length, line_data).
    """
    if not os.path.isfile(paf_file):
        sys.exit(f"Error: File not found: {paf_file}")

    longest_alignments={}

    try:
        with open(paf_file, "r") as f:
            for line_num, line in enumerate(f, start=1):
                hit = line.rstrip().split('\t')
                if len(hit) < 12:
                  sys.exit(f"Error: Invalid PAF format on line {line_num}: Not enough fields.")
 
                try:
                    target_start = int(hit[7])
                    target_end = int(hit[8])           
                    query_name = hit[0]
                    target_name = hit[5]
                    mqual = int(hit[11])
                    strand = hit[4]
                    
                    if match:
                        if match not in query_name or "unloc" in line or target_name in sex_chroms or query_name in sex_chroms:
                            continue
                    
                except ValueError:
                     sys.exit(f"Error: Non-integer target coordinate on line {line_num}.")

                alignment_length = target_end - target_start 
                if alignment_length < 0:
                    sys.exit(f"Error: Negative alignment length on line {line_num}.")
                if mqual != 60:
                    continue
              
                if (query_name not in longest_alignments or alignment_length > longest_alignments[query_name][1]):
                    longest_alignments[query_name] = (target_name, alignment_length, strand)

    except Exception as e:
        sys.exit(f"Error reading file: {e}")

    return longest_alignments

def process_fasta(fasta_path, alignments, output_fasta_path):

    query_dict = {}
    correspondence_data = []

    try:
        with open(fasta_path, 'r') as f:
            current_seq = ''
            current_name = ''
            for line in f:
                line = line.strip()
                if line.startswith('>'):  # Header line
                    final_seq = current_seq
                    if current_name:
                        if current_name in alignments and '-' in alignments[current_name]:
                            final_seq = revcomp(current_seq) 
                        query_dict[current_name] = format_sequence(final_seq)
                    current_name = line[1:]  # Remove '>'
                    current_seq = ''
                else:
                    current_seq += line 
            
            if current_name:
                final_seq = current_seq
                if current_name in alignments and '-' in alignments[current_name]:
                    final_seq = revcomp(current_seq) 
                query_dict[current_name] = format_sequence(final_seq)

    except Exception as e:
        sys.exit(f"Error reading FASTA file: {e}")

    for i in alignments:
        print (i, alignments[i])
    for keys in query_dict:
        if keys in alignments:
            args.fasta_output.write(">" + alignments[keys][0] + '\n')
        elif "unloc" in keys:
            name = keys.split("_unloc")
            if name[0] in alignments:
                args.fasta_output.write(">" + alignments[name[0]][0] + "_unloc" + name[1] + '\n')
            else:
                args.fasta_output.write(">" + keys+ '\n')
        else: 
            args.fasta_output.write(">" + keys+ '\n')
        args.fasta_output.write(query_dict[keys] + "\n")

def revcomp(sequence):
    complement = {
        'A': 'T', 'T': 'A',
        'C': 'G', 'G': 'C',
        'a': 't', 't': 'a',
        'c': 'g', 'g': 'c',
        'N': 'N', 'n': 'n'  # N stands for any nucleotide
    }

    # Create the reverse complement
    rev_comp = ''.join(complement.get(base, base) for base in reversed(sequence))

    return rev_comp

def format_sequence(seq, line_length=60):
    """Split a sequence into lines of specified length (default: 60)."""
    return '\n'.join(seq[i:i+line_length] for i in range(0, len(seq), line_length))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--paf-file",  required=True, help="PAF file from minimap2. Names in the query assembly will be replaced by corresponding sequence names in the target assembly.")
    parser.add_argument("-f", "--fasta-file", required=False, help="fasta query file")
    parser.add_argument("-x", "--sex_chroms", required=False, nargs="+", help="list of chromosomes that you do not want to rename")
    parser.add_argument("-k", "--lookup-table", required=True, help="Output file that will contain the correspondance between the two assemblies")
    parser.add_argument("-o", "--fasta-output", required=False, type=argparse.FileType('w'), default=sys.stdout,
                        help="Output fasta file that will contain the renamed sequences, reverse complemented if aligned in the negative strand. Default: stdout")
    parser.add_argument("-m", "--match_seqs", required=False, default="SUPER", help="Match only sequences with this prefix. Default: SUPER" )
    args = parser.parse_args()
    
    paf_path = args.paf_file
    fasta_path, output_fasta_path = "", ""
    if args.fasta_file:
        fasta_path = args.fasta_file
    if args.fasta_output:
        output_fasta_path = args.fasta_output
    lookup_path = args.lookup_table
    match=args.match_seqs
    sex_chroms = []
    if args.sex_chroms:
        sex_chroms = args.sex_chroms
    alignments = parse_paf_and_get_longest_alignments(paf_path, match, sex_chroms)


    if fasta_path:
        process_fasta(fasta_path, alignments, args.fasta_output)

    if args.fasta_output is not sys.stdout:
        args.fasta_output.close()

    with open(lookup_path, 'w') as f:
        for target, (query, length, strand) in alignments.items():
            try:
                f.write(f"{target}\t{query}\t{strand}\n")
            except Exception as e:
                sys.exit(f"Error writing correspondence file: {e}")



