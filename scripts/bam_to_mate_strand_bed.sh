#!/bin/bash

# This script processes a BAM file to generate strand-specific and mate-specific
# coverage files in a simple 3-column BED format (chromosome, position, coverage).
# The script is compatible only with paired-end data and assumes the input BAM
# file is sorted by coordinates. No index is used for input BAM.
# The coverage/depth can be performed at 2 resolution: read or last aligned read position.
#
# The script generates two main sets of outputs:
# 1. Properly mapped reads, separated by strand (forward or reverse).
# 2. Reads mapped without a mate, separated by strand and mate identity (first or second).
#
# Note on file naming for strand:
#   The strand designation (forward or reverse) is chosen arbitrarily and does not depend on
#   the library type. However, since the same nomenclature is consistently used across files, 
#   this does not affect the analysis (in fact this denomination would correspond to a 'fr' 
#   library and it would be the opposite for an 'rf' library).
#
#   This script is **not compatible** with 'rr' or 'ff' libraries.
#
# Arguments:
# -i : Path to the input BAM file.
# -o : Directory to save the output files.
# -l : Depth count level: consider the entire reads ['read'] or simply the last aligned position of the reads ['position'] (last in the context of the last cycle of the read, i.e. the closest to the other mate for fr or rf libraries when inner distance is positive). (default: read)"
# -d : Optional; if set to "false", keeps intermediate BAM files (default: true, meaning BAM files are deleted).
#
# Usage:
#   bash process_bam.sh -i <input_bam> -o <output_dir> [option(s)]
#
# Requirements:
# - `samtools` for BAM manipulation.
# - `bedtools` for creating coverage files in BED format. (if '-l read')
# - `pysam` on `python3` for creating coverage files in BED format. (if '-l position' - calling `./script/bam_to_start_end_depth.py`)

set -e

# help
usage() {
    echo "Usage: $0 -i <input_bam> -o <output_dir> [-d true|false]"
    echo "Options:"
    echo "  -i    Path to the input BAM file"
    echo "  -o    Path to the output directory. Must not exist."
    echo "  -l    Depth count level: consider the entire reads ['read'] or simply the last aligned position of the reads ['position'] (last in the context of the last cycle of the read, i.e. the closest to the other mate for fr or rf libraries when inner distance is positive). (default: read)"
    echo "  -d    Delete intermediate BAM files (default: true)"
    echo "  -s    Path of scripts directories, must contains 'bam_to_start_end_depth.py' (used if counting at position reslution instead of read resolution) (default: ./scripts/)."
    echo "NOTE: - compatible only about paired-end data"
    echo "      - input bam must be sorted by coordinates"
    echo "      - no need bam index"
    exit 1
}

# args
delete_bam_files="true"
count_level="read"
scripts_path="./scripts/"
while getopts "i:o:l:d:s:" opt; do
    case $opt in
        i) input_bam="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        l) count_level="$OPTARG" ;;
        d) delete_bam_files="$OPTARG" ;;
        s) scripts_path="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "$input_bam" || -z "$output_dir" ]]; then
    echo "Error: -i and -o options are required."
    usage
fi

if [[ "$count_level" != "position" && "$count_level" != "read" ]]; then
    echo "Error: Invalid value for -l option. Must be 'position' or 'read' (your value: '${count_level}')."
    usage
fi

mkdir "${output_dir}"

#################################################
#### PROPERLY MAPPED BY STRAND (R1 + R2)
#################################################
echo -ne "\n#### PROPERLY MAPPED BY STRAND"
out_dir="${output_dir}/properly_map_by_strand/"
mkdir "${out_dir}"

### properly mapped & in strand
# first in strand & properly mapped: (-f64 -F16) & (-f2)
echo -ne " . "
samtools view -h -b -f66 -F16 -o "${out_dir}/properly_map_first_in_strand.bam" "${input_bam}"
# second in strand & properly mapped: (-f128 -f16) & (-f2)
echo -ne " . "
samtools view -h -b -f146 -o "${out_dir}/properly_map_second_in_strand.bam" "${input_bam}"
# merge bam: first and second
echo -ne " . "
samtools merge -o "${out_dir}/properly_map_first_and_second_in_strand.bam" "${out_dir}/properly_map_first_in_strand.bam" "${out_dir}/properly_map_second_in_strand.bam"
# bam to bed
if [[ "$count_level" == "read" ]]; then
    bedtools genomecov -ibam "${out_dir}/properly_map_first_and_second_in_strand.bam" -d -split > "${out_dir}/properly_map_first_and_second_in_strand.bed"
elif [[ "$count_level" == "position" ]]; then
    samtools index "${out_dir}/properly_map_first_in_strand.bam"
    ${scripts_path}bam_to_start_end_depth.py --pos last "${out_dir}/properly_map_first_in_strand.bam" "${out_dir}/properly_map_first_in_strand"
    samtools index "${out_dir}/properly_map_second_in_strand.bam"
    ${scripts_path}bam_to_start_end_depth.py --pos first "${out_dir}/properly_map_second_in_strand.bam" "${out_dir}/properly_map_second_in_strand"
    # merge bed: first and second
    paste "${out_dir}/properly_map_first_in_strand_last_aln_pos_depth.bed" "${out_dir}/properly_map_second_in_strand_first_aln_pos_depth.bed" | \
      while read -r chr1 pos1 cov1 chr2 pos2 cov2; do \
        echo -e "$chr1\t$pos1\t$((cov1 + cov2))"; \
      done > "${out_dir}/properly_map_first_and_second_in_strand.bed"
fi

### properly mapped & in reverse strand
# first in reverse & properly mapped: (-f64 -f16) & (-f2)
echo -ne " . "
samtools view -h -b -f82 -o "${out_dir}/properly_map_first_in_reverse.bam" "${input_bam}"
# second in reverse & properly mapped: (-f128 -F16) & (-f2)
echo -ne " . "
samtools view -h -b -f130 -F16 -o "${out_dir}/properly_map_second_in_reverse.bam" "${input_bam}"
# merge first and second
echo -ne " . "
samtools merge -o "${out_dir}/properly_map_first_and_second_in_reverse.bam" "${out_dir}/properly_map_first_in_reverse.bam" "${out_dir}/properly_map_second_in_reverse.bam"
# bam to bed
if [[ "$count_level" == "read" ]]; then
    bedtools genomecov -ibam "${out_dir}/properly_map_first_and_second_in_reverse.bam" -d -split > "${out_dir}/properly_map_first_and_second_in_reverse.bed"
elif [[ "$count_level" == "position" ]]; then
    samtools index "${out_dir}/properly_map_first_in_reverse.bam"
    ${scripts_path}bam_to_start_end_depth.py --pos first "${out_dir}/properly_map_first_in_reverse.bam" "${out_dir}/properly_map_first_in_reverse"
    samtools index "${out_dir}/properly_map_second_in_reverse.bam"
    ${scripts_path}bam_to_start_end_depth.py --pos last "${out_dir}/properly_map_second_in_reverse.bam" "${out_dir}/properly_map_second_in_reverse"
    # merge bed: first and second
    paste "${out_dir}/properly_map_first_in_reverse_first_aln_pos_depth.bed" "${out_dir}/properly_map_second_in_reverse_last_aln_pos_depth.bed" | \
      while read -r chr1 pos1 cov1 chr2 pos2 cov2; do \
        echo -e "$chr1\t$pos1\t$((cov1 + cov2))"; \
      done > "${out_dir}/properly_map_first_and_second_in_reverse.bed"
fi

if [ "$delete_bam_files" == "true" ]; then
    rm "${out_dir}"/*.bam
fi

#################################################
#### MAPPED WITHOUT MATE BY READS AND STRAND
#################################################
echo -ne "\n\n#### MAPPED WITHOUT MATE BY READS AND STRAND"
out_dir="${output_dir}/single_mate_map_by_mate_and_strand/"
mkdir "${out_dir}"

# first in strand & only first map: (-f64 -F16) & (-F4 -f64 -f8)
echo -ne " . "
samtools view -h -b -f72 -F20 -o "${out_dir}/only_first_map_in_strand.bam" "${input_bam}"
if [[ "$count_level" == "read" ]]; then
    bedtools genomecov -ibam "${out_dir}/only_first_map_in_strand.bam" -d -split > "${out_dir}/only_first_map_in_strand.bed"
elif [[ "$count_level" == "position" ]]; then
    samtools index "${out_dir}/only_first_map_in_strand.bam"
    ${scripts_path}bam_to_start_end_depth.py --pos last "${out_dir}/only_first_map_in_strand.bam" "${out_dir}/only_first_map_in_strand"
    mv "${out_dir}/only_first_map_in_strand_last_aln_pos_depth.bed" "${out_dir}/only_first_map_in_strand.bed"
fi
# second in strand & only second map: (-f128 -f16) & (-F4 -f128 -f8)
echo -ne " . "
samtools view -h -b -f152 -F4 -o "${out_dir}/only_second_map_in_strand.bam" "${input_bam}"
if [[ "$count_level" == "read" ]]; then
    bedtools genomecov -ibam "${out_dir}/only_second_map_in_strand.bam" -d -split > "${out_dir}/only_second_map_in_strand.bed"
elif [[ "$count_level" == "position" ]]; then
    samtools index "${out_dir}/only_second_map_in_strand.bam"
    ${scripts_path}bam_to_start_end_depth.py --pos first "${out_dir}/only_second_map_in_strand.bam" "${out_dir}/only_second_map_in_strand"
    mv "${out_dir}/only_second_map_in_strand_first_aln_pos_depth.bed" "${out_dir}/only_second_map_in_strand.bed"
fi
# first in reverse & only first map: (-f64 -f16) & (-F4 -f64 -f8)
echo -ne " . "
samtools view -h -b -f88 -F4 -o "${out_dir}/only_first_map_in_reverse.bam" "${input_bam}"
if [[ "$count_level" == "read" ]]; then
    bedtools genomecov -ibam "${out_dir}/only_first_map_in_reverse.bam" -d -split > "${out_dir}/only_first_map_in_reverse.bed"
elif [[ "$count_level" == "position" ]]; then
    samtools index "${out_dir}/only_first_map_in_reverse.bam"
    ${scripts_path}bam_to_start_end_depth.py --pos first "${out_dir}/only_first_map_in_reverse.bam" "${out_dir}/only_first_map_in_reverse"
    mv "${out_dir}/only_first_map_in_reverse_first_aln_pos_depth.bed" "${out_dir}/only_first_map_in_reverse.bed"
fi
# second in reverse & only second map: (-f128 -F16) & (-F4 -f128 -f8)
echo -ne " . "
samtools view -h -b -f136 -F20 -o "${out_dir}/only_second_map_in_reverse.bam" "${input_bam}"
if [[ "$count_level" == "read" ]]; then
    bedtools genomecov -ibam "${out_dir}/only_second_map_in_reverse.bam" -d -split > "${out_dir}/only_second_map_in_reverse.bed"
elif [[ "$count_level" == "position" ]]; then
    samtools index "${out_dir}/only_second_map_in_reverse.bam"
    ${scripts_path}bam_to_start_end_depth.py --pos last "${out_dir}/only_second_map_in_reverse.bam" "${out_dir}/only_second_map_in_reverse"
    mv "${out_dir}/only_second_map_in_reverse_last_aln_pos_depth.bed" "${out_dir}/only_second_map_in_reverse.bed"
fi

if [ "$delete_bam_files" == "true" ]; then
    rm "${out_dir}"/*.bam
fi

echo -ne "\n\n"