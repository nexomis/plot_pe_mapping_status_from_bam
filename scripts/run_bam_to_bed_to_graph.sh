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
min_x_graduation=200
output_depth_pdf=""
input_bam=""
output_dir=""
plotly_out_dir=""
annot_gff_file=""
annot_feat_id_regex='.*Target=([^; ]+).*'
annot_feat_id_catch="\\1"
annot_interest_type="mRNA,CDS,transcript"
annot_name_on_plot=FALSE
relative_heigt_combined_plot="3,1"
max_pb_by_A4_width=""
pos_min=""
pos_max=""
area="FALSE"

# help
usage() {
    echo "Usage: $0 -input_bam <input_bam> -output_dir <output_dir> -output_bn_depth_pdf <output_bn_depth_pdf> [-d <delete_bam_files>] [-p <plot_title>]"
    echo ""
    echo "This script processes a BAM file to generate specific BED files (by calling of 'bam_to_mate_strand_bed.sh') and corresponding depth plots (by calling of 'plot_multiple_bed_depth.r')."
    echo ""
    echo "Arguments:"
    echo ""
    echo "  >global"
    echo "    -input_bam         Recquired. Path to the input BAM file. Must be sorted by coordinates."
    echo "    -output_dir        Recquired. Path to the output directory for BED files and plots."
    echo "    -scripts_dir_path  Optional.  Path of scripts directories if not in PATH. Must contains 'bam_to_mate_strand_bed.sh', 'plot_multiple_bed_depth.r' and 'bam_to_start_end_depth.py' (default: NULL)."
    echo ""
    echo "  >depth calculation"
    echo "    -count_level       Optional.  Depth count level: consider the entire reads ['read'] or simply the last aligned position of the reads ['position'] (last in the context of the last cycle of the read, i.e. the closest to the other mate for fr or rf libraries when inner distance is positive). (default: read)"
    echo "    -delete_bam_files  Optional.  true|false Whether to delete BAM files after processing (default: true)."
    echo ""
    echo "  >depth graph"
    echo "    -plot_title           Optional.  Title for the depth plots (default: 'Depth')."
    echo "    -output_bn_depth_pdf  Recquired. Base name to the output PDF file for depth plots."
    echo "    -min_x_graduation     Optional.  Minimum x-axis graduation in base pairs (default: 200)"
    echo "    -plotly_out_dir       Optional.  Directory to save individual chromosome plots as HTML (interactive Plotly plots) (if provided) (default: NA)"
    echo "    -max_pb_by_A4_width   Optional.  Maximum plot width in pixels for automatic page sizing (default: 5000)"
    echo "    -pos_min              Optional.  First position to consider for plot (if specified, applied for all chromosomes) (default: '0')"
    echo "    -pos_max              Optional.  Last position to consider for plot (if specified, applied to all chromosomes) (default: last position of each chromosome)"
    echo "    -area                 Optional.  Boolean. If 'TRUE', draw depth graphics using 'geom_area()' instead of 'geom_line()'. (default: FALSE)"
    echo ""
    echo "  >annot graph"
    echo "    -annot_gff_file                Optional. Path to input GFF file for annotation plotting (if provided) (default: NA)."
    echo "    -annot_feat_id_regex           Optional. Pattern (sub) to extract feature ID from the GFF attributes (default: '.*Target=([^; ]+).*')."
    echo "    -annot_feat_id_catch           Optional. Replacement (sub) for the feature ID (default: '\\1')."
    echo "    -annot_interest_type           Optional. Vector of feature types to include (default: 'mRNA,CDS,transcript')."
    echo "    -annot_name_on_plot            Optional. Logical; displays feature names on the plot (default: FALSE)"
    echo "    -relative_heigt_combined_plot  Optional. Relative heights for the combined depth and annotation plot (default: '3,1')."
    echo ""
    exit 1
}


# parse args
count_level="read"
delete_bam_files="true"
plot_title="Depth"
scripts_dir_path="./scripts/"
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -input_bam) input_bam="$2"; shift ;;
        -output_dir) output_dir="$2"; shift ;;
        -count_level) count_level="$2"; shift ;;
        -delete_bam_files) delete_bam_files="$2"; shift ;;
        -plot_title) plot_title="$2"; shift ;;
        -output_bn_depth_pdf) output_bn_depth_pdf="$2"; shift ;;
        -scripts_dir_path) scripts_dir_path="$2"; shift ;;
        -plotly_out_dir) plotly_out_dir="$2"; shift ;;
        -annot_gff_file) annot_gff_file="$2"; shift ;;
        -max_pb_by_A4_width) max_pb_by_A4_width="$2"; shift ;;
        -pos_min) pos_min="$2"; shift ;;
        -pos_max) pos_max="$2"; shift ;;
        -area) area="$2"; shift ;;
        -min_x_graduation) min_x_graduation="$2"; shift ;;
        -annot_feat_id_regex) annot_feat_id_regex="$2"; shift ;;
        -annot_feat_id_catch) annot_feat_id_catch="$2"; shift ;;
        -annot_interest_type) annot_interest_type="$2"; shift ;;
        -annot_name_on_plot) annot_name_on_plot="$2"; shift ;;
        -relative_heigt_combined_plot) relative_heigt_combined_plot="$2"; shift ;;
        *) usage ;;
    esac
    shift
done


# check args
if [ -z "$input_bam" ] || [ -z "$output_dir" ] || [ -z "$output_bn_depth_pdf" ]; then
    echo "Error: Missing mandatory arguments. '$output_bn_depth_pdf'"
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
static_plot_parameters='--depth_axes="2,1,1" --plot_legend="properly_map(first+second),only_first_map,only_second_map" --plot_colors="grey50,darkblue,darkred" --plot_alpha="0.4"'

dynamic_plot_parameters=""
if [ -n "$plotly_out_dir" ]; then dynamic_plot_parameters+=" --plotly_out_dir=${plotly_out_dir}"; fi
if [ -n "$max_pb_by_A4_width" ]; then dynamic_plot_parameters+=" --max_pb_by_A4_width=${max_pb_by_A4_width}"; fi
if [ -n "$pos_min" ]; then dynamic_plot_parameters+=" --pos_min=${pos_min}"; fi
if [ -n "$pos_max" ]; then dynamic_plot_parameters+=" --pos_max=${pos_max}"; fi

if [ -n "$annot_gff_file" ]; then
    dynamic_plot_parameters+=" --annot_gff_file=\"${annot_gff_file}\" \
 --annot_feat_id_regex='${annot_feat_id_regex}' \
 --annot_feat_id_catch='${annot_feat_id_catch}' \
 --annot_interest_type=\"${annot_interest_type}\" \
 --annot_name_on_plot=\"${annot_name_on_plot}\" \
 --relative_heigt_combined_plot=\"${relative_heigt_combined_plot}\""
fi

echo "Rscript ${scripts_dir_path}plot_multiple_bed_depth.r \
  --bed_files=\"${output_dir}/properly_map_by_strand/properly_map_first_and_second_in_reverse.bed,${output_dir}/single_mate_map_by_mate_and_strand/only_first_map_in_reverse.bed,${output_dir}/single_mate_map_by_mate_and_strand/only_second_map_in_reverse.bed\" \
  --bed_files_opposite=\"${output_dir}/properly_map_by_strand/properly_map_first_and_second_in_strand.bed,${output_dir}/single_mate_map_by_mate_and_strand/only_first_map_in_strand.bed,${output_dir}/single_mate_map_by_mate_and_strand/only_second_map_in_strand.bed\" \
  --output_pdf_bn=\"${output_dir}/${output_bn_depth_pdf}\" \
  --plot_title=\"${plot_title}\" \
  --min_x_graduation=\"${min_x_graduation}\" \
  --area=\"${area}\" \
  ${static_plot_parameters} \
  ${dynamic_plot_parameters}" | bash

echo "Process completed."
