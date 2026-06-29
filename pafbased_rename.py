#!/usr/bin/env python3

import sys
import os
import argparse
import re

#Author: Jessica Gomez-Garrido, CNAG.
#Contact email: jessica.gomez@cnag.eu
#Date:20250612

def parse_paf_and_choose_orientation(paf_file, match, sex_chroms):
    """
    Parse a PAF and, for each query sequence, decide its dominant target and the
    orientation that maximizes forward (positive) colinearity with the target
    assembly.

    Orientation is decided by the TOTAL high-confidence aligned length per strand,
    not by the single longest block. So if a query scaffold is a join of two
    target superscaffolds, the whole scaffold is flipped to keep the *larger*
    amount of aligned sequence forward-colinear (the smaller joined piece may end
    up inverse, which is acceptable). The dominant target name is the one with the
    most total aligned length.

    Args:
        paf_file (str): Path to the PAF file.

    Returns:
        dict: query_name -> (dominant_target_name, total_aligned_length, strand)
              where strand is '+' or '-' (the orientation to apply to the query).
    """
    if not os.path.isfile(paf_file):
        sys.exit(f"Error: File not found: {paf_file}")

    target_len = {}    # query -> {target_name: total aligned length}
    strand_len = {}    # query -> {'+': total fwd length, '-': total rev length}

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
                if strand not in ('+', '-'):
                    continue

                target_len.setdefault(query_name, {})
                target_len[query_name][target_name] = \
                    target_len[query_name].get(target_name, 0) + alignment_length
                strand_len.setdefault(query_name, {'+': 0, '-': 0})
                strand_len[query_name][strand] += alignment_length

    except Exception as e:
        sys.exit(f"Error reading file: {e}")

    alignments = {}
    for query_name in target_len:
        dominant_target = max(target_len[query_name], key=target_len[query_name].get)
        total_len = sum(target_len[query_name].values())
        orientation = '-' if strand_len[query_name]['-'] > strand_len[query_name]['+'] else '+'
        alignments[query_name] = (dominant_target, total_len, orientation)

    return alignments

def process_fasta(fasta_path, alignments, output_fasta_path, keep_names=False):

    query_dict = {}

    def reverse_if_needed(name, seq):
        # orientation is the 3rd element of the tuple ('+' or '-')
        if name in alignments and alignments[name][2] == '-':
            return revcomp(seq)
        return seq

    try:
        with open(fasta_path, 'r') as f:
            current_seq = ''
            current_name = ''
            for line in f:
                line = line.strip()
                if line.startswith('>'):  # Header line
                    if current_name:
                        query_dict[current_name] = format_sequence(reverse_if_needed(current_name, current_seq))
                    current_name = line[1:].split()[0]  # Remove '>' and any description
                    current_seq = ''
                else:
                    current_seq += line

            if current_name:
                query_dict[current_name] = format_sequence(reverse_if_needed(current_name, current_seq))

    except Exception as e:
        sys.exit(f"Error reading FASTA file: {e}")

    for i in alignments:
        print(i, alignments[i], file=sys.stderr)
    for keys in query_dict:
        if keep_names:
            out_name = keys
        elif keys in alignments:
            out_name = alignments[keys][0]
        elif "unloc" in keys:
            name = keys.split("_unloc")
            if name[0] in alignments:
                out_name = alignments[name[0]][0] + "_unloc" + name[1]
            else:
                out_name = keys
        else:
            out_name = keys
        args.fasta_output.write(">" + out_name + '\n')
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
    parser.add_argument("-K", "--keep-names", action="store_true", help="Keep the query sequence names instead of renaming to the target; still reorient (reverse-complement) to match the target assembly.")
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
    alignments = parse_paf_and_choose_orientation(paf_path, match, sex_chroms)


    if fasta_path:
        process_fasta(fasta_path, alignments, args.fasta_output, keep_names=args.keep_names)

    if args.fasta_output is not sys.stdout:
        args.fasta_output.close()

    with open(lookup_path, 'w') as f:
        for target, (query, length, strand) in alignments.items():
            try:
                f.write(f"{target}\t{query}\t{strand}\n")
            except Exception as e:
                sys.exit(f"Error writing correspondence file: {e}")



