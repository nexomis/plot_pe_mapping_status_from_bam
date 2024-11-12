#!/bin/bash

set -e

# This script automates the process of:
#  - Converting a BAM file into specific BED files: It calls 'bam_to_mate_strand_bed.sh'
#    script to convert the input BAM file into multiple BED files, representing
#    depth data for each strand. (recquires 'samtools' and 'bedtools')
#  - Generating depth plots: After converting the BAM file, it runs 
#    'plot_multiple_bed_depth.r' script to generate depth plots from the 
#    resulting BED files, and saves the plots as a PDF. (recquires 'ggplot2' on r)
# Usage:
# ./run_bam_to_bed_to_graph.sh -i <input_bam> -o <output_dir> -c <figure_basename_pdf> [options]

# Default values
delete_bam_files="true"
plot_title="Depth"
output_cov_pdf=""
input_bam=""
output_dir=""

# help
usage() {
    echo "Usage: $0 -i <input_bam> -o <output_dir> -c <output_bn_cov_pdf> [-d <delete_bam_files>] [-p <plot_title>]"
    echo ""
    echo "This script processes a BAM file to generate specific BED files (by calling of 'bam_to_mate_strand_bed.sh') and corresponding depth plots (by calling of 'plot_multiple_bed_depth.r')."
    echo ""
    echo "Arguments:"
    echo "  -i <input_bam>           Recquired. Path to the input BAM file. Must be sorted by coordinates."
    echo "  -o <output_dir>          Recquired. Path to the output directory for BED files and plots."
    echo "  -l <count_level>         Optional. Depth count level: consider the entire reads ['read'] or simply the last aligned position of the reads ['position'] (last in the context of the last cycle of the read, i.e. the closest to the other mate for fr or rf libraries when inner distance is positive). (default: read)"
    echo "  -d <delete_bam_files>    Optional. true|false Whether to delete BAM files after processing (default: true)."
    echo "  -p <plot_title>          Optional. Title for the depth plots (default: 'Depth')."
    echo "  -c <output_bn_cov_pdf>   Recquired. Base name to the output PDF file for depth plots."
    echo "  -s <scripts_dir_path>    Optional. Path of scripts directories. Must contains 'bam_to_mate_strand_bed.sh', 'plot_multiple_bed_depth.r' and 'bam_to_start_end_depth.py' (default: ./scripts/)."
    exit 1
}


# parse args
count_level="read"
delete_bam_files="true"
plot_title="Depth"
scripts_dir_path="./scripts/"
while getopts "i:o:l:d:p:c:s:" opt; do
    case $opt in
        i) input_bam="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        l) count_level="$OPTARG" ;;
        d) delete_bam_files="$OPTARG" ;;
        p) plot_title="$OPTARG" ;;
        c) output_bn_cov_pdf="$OPTARG" ;;
        s) scripts_dir_path="$OPTARG" ;;
        *) usage ;;
    esac
done


# check args
if [ -z "$input_bam" ] || [ -z "$output_dir" ] || [ -z "$output_bn_cov_pdf" ]; then
    echo "Error: Missing mandatory arguments. '$output_bn_cov_pdf'"
    usage
fi

if [ "$delete_bam_files" != "true" ] && [ "$delete_bam_files" != "false" ]; then
    echo "Error: Invalid value for -d (delete_bam_files). It must be 'true' or 'false'."
    usage
fi

# run 'bam_to_mate_strand_bed.sh'
echo "Generating BED files from BAM..."
${scripts_dir_path}bam_to_mate_strand_bed.sh \
  -i "${input_bam}" \
  -o "${output_dir}" \
  -d "${delete_bam_files}" \
  -s "${scripts_dir_path}" \
  -l "${count_level}"


# run 'plot_multiple_bed_depth.r'
echo "Generating depth plots from BED files..."
Rscript ${scripts_dir_path}plot_multiple_bed_depth.r \
  --bed_files="${output_dir}/single_mate_map_by_mate_and_strand/only_first_map_in_reverse.bed,${output_dir}/single_mate_map_by_mate_and_strand/only_second_map_in_reverse.bed,${output_dir}/properly_map_by_strand/properly_map_first_and_second_in_reverse.bed" \
  --bed_files_opposite="${output_dir}/single_mate_map_by_mate_and_strand/only_first_map_in_strand.bed,${output_dir}/single_mate_map_by_mate_and_strand/only_second_map_in_strand.bed,${output_dir}/properly_map_by_strand/properly_map_first_and_second_in_strand.bed" \
  --output_pdf="${output_dir}/${output_bn_cov_pdf}.pdf" \
  --plot_title="${plot_title}" \
  --cov_axes="1,1,2" \
  --plot_legend="only_first_map,only_second_map,properly_map(first+second)"


echo "Process completed."
