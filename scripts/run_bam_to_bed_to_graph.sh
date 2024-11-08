#!/bin/bash

set -e

# This script automates the process of:
#  - Converting a BAM file into specific BED files: It calls 'bam_to_mate_strand_bed.sh'
#    script to convert the input BAM file into multiple BED files, representing
#    coverage data for each strand. (recquires 'samtools' and 'bedtools')
#  - Generating coverage plots: After converting the BAM file, it runs 
#    'plot_multiple_bed_coverage.r' script to generate coverage plots from the 
#    resulting BED files, and saves the plots as a PDF. (recquires 'ggplot2' on r)
# Usage:
# ./run_bam_to_bed_to_graph.sh --input_bam <input_bam> --output_dir <output_dir> [--delete_bam <yes|no>] [--plot_title <title>] [--output_cov_pdf <output_pdf>]

# Default values
delete_bam_files="true"
plot_title="Coverage"
output_cov_pdf=""
input_bam=""
output_dir=""

# help
usage() {
    echo "Usage: $0 -i <input_bam> -o <output_dir> -c <output_bn_cov_pdf> [-d <delete_bam_files>] [-p <plot_title>]"
    echo ""
    echo "This script processes a BAM file to generate specific BED files (by calling of 'bam_to_mate_strand_bed.sh') and corresponding coverage plots (by calling of 'plot_multiple_bed_coverage.r')."
    echo ""
    echo "Arguments:"
    echo "  -i <input_bam>           Path to the input BAM file. Must be sorted by coordinates (mandatory)"
    echo "  -o <output_dir>          Path to the output directory for BED files and plots (mandatory)"
    echo "  -d <delete_bam_files>    true|false Whether to delete BAM files after processing (default: true)"
    echo "  -p <plot_title>          Title for the coverage plots (default: 'Coverage')"
    echo "  -c <output_bn_cov_pdf>   Base name to the output PDF file for coverage plots (mandatory)"
    echo "  -h                       Display this help message"
    exit 1
}


# parse args
while getopts "i:o:d:p:c:" opt; do
    case $opt in
        i) input_bam="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        d) delete_bam_files="$OPTARG" ;;
        p) plot_title="$OPTARG" ;;
        c) output_bn_cov_pdf="$OPTARG" ;;
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
./scripts/bam_to_mate_strand_bed.sh \
  -i "${input_bam}" \
  -o "${output_dir}" \
  -d "${delete_bam_files}"


# run 'plot_multiple_bed_coverage.r'
echo "Generating coverage plots from BED files..."
Rscript ./scripts/plot_multiple_bed_coverage.r \
  --bed_files="${output_dir}/single_mate_map_by_mate_and_strand/only_first_map_in_reverse.bed,${output_dir}/single_mate_map_by_mate_and_strand/only_second_map_in_reverse.bed,${output_dir}/properly_map_by_strand/properly_map_first_and_second_in_reverse.bed" \
  --bed_files_opposite="${output_dir}/single_mate_map_by_mate_and_strand/only_first_map_in_strand.bed,${output_dir}/single_mate_map_by_mate_and_strand/only_second_map_in_strand.bed,${output_dir}/properly_map_by_strand/properly_map_first_and_second_in_strand.bed" \
  --output_pdf="${output_dir}/${output_bn_cov_pdf}.pdf" \
  --plot_title="${plot_title}" \
  --cov_axes="1,1,2" \
  --plot_legend="only_first_map,only_second_map,properly_map(first+second)"


echo "Process completed."
