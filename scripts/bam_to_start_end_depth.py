#!/usr/bin/env python3


import pysam
from collections import defaultdict
import argparse


# args
parser = argparse.ArgumentParser(description="Compute read depth on alignment positions (first, last, or both) for each position in BAM file.")
parser.add_argument("input_bam", help="Input BAM file")
parser.add_argument("output_prefix", help="Output prefix for BED files")
parser.add_argument("--pos", choices=["first", "last", "both"], default="both", 
                    help="Compute depth on first alignment position, last, or both (default: last)")
args = parser.parse_args()

bamfile = pysam.AlignmentFile(args.input_bam, "rb")

# init depth counters
first_depth = defaultdict(lambda: defaultdict(int))
last_depth = defaultdict(lambda: defaultdict(int))

# extract chosen alignment position and compute depth
for read in bamfile.fetch():  # default: only mapped reads (cf. 'until_eof' parameters)
    chr_name = bamfile.get_reference_name(read.reference_id)
    
    if args.pos in ["last", "both"]:
        last_pos = read.reference_end
        last_depth[chr_name][last_pos] += 1
    
    if args.pos in ["first", "both"]:
        first_pos = read.reference_start + 1
        first_depth[chr_name][first_pos] += 1

# write results on bed file(s)
def write_bed(file_path, depth_dict):
    with open(file_path, "w") as bedfile:
        for chr_name in bamfile.references:
            chr_length = bamfile.get_reference_length(chr_name)
            for pos in range(1, chr_length + 1):
                depth = depth_dict[chr_name].get(pos, 0)
                bedfile.write(f"{chr_name}\t{pos}\t{depth}\n")

if args.pos in ["last", "both"]:
    write_bed(f"{args.output_prefix}_last_aln_pos_depth.bed", last_depth)

if args.pos in ["first", "both"]:
    write_bed(f"{args.output_prefix}_first_aln_pos_depth.bed", first_depth)

bamfile.close()
