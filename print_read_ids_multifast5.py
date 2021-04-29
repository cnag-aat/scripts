#!/usr/bin/env python
from ont_fast5_api.fast5_interface import get_fast5_file
from argparse import ArgumentParser

parser = ArgumentParser("")
parser.add_argument('-i', '--input_file', required=True,
                        help="MultiRead fast5 file")
args = parser.parse_args()
fast5_filepath = args.input_file # This can be a single- or multi-read file
with get_fast5_file(fast5_filepath, mode="r") as f5:
    for read in f5.get_reads():
        #raw_data = read.get_raw_data()
        print(read.read_id + "\t" + fast5_filepath)

