## PLOT PAIRED-END MAPPING STATUS FROM BAM

This script automates the process of:  
 - Converting a BAM file into specific BED files: It calls `bam_to_mate_strand_bed.sh` script to convert the input BAM file into multiple BED files, representing depth data for each strand. (recquires `samtools` (and `bedtools` if depth count is performed at the resolution of read)). *This script itself calls 'bam_to_start_end_depth.py' if the depth count is performed at the resolution of the last aligned position (recquires `pysam` on `python3`).*
 - Generating depth plots: After converting the BAM file, it runs `plot_multiple_bed_depth.r` script to generate depth plots from the resulting BED files, and saves the plots as a PDF. (recquires `ggplot2` on `r`)  

### Usage:
  `./run_bam_to_bed_to_graph.sh --input_bam <input_bam> --output_dir <output_dir> [option(s)]`  

### Arguments:
```
  -i <input_bam>           Recquired. Path to the input BAM file. Must be sorted by coordinates (mandatory)  
  -o <output_dir>          Recquired. Path to the output directory for BED files and plots (mandatory)  
  -c <output_bn_cov_pdf>   Recquired. Base name to the output PDF file for depth plots (mandatory)  
  -d <delete_bam_files>    Optional. true|false Whether to delete BAM files after processing (default: true)  
  -p <plot_title>          Optional. Title for the depth plots (default: 'Depth')  
  -l <count_level>         Optional. Depth count level: consider the entire reads ['read'] or simply the last aligned position of the reads ['position'] (last in the context of the last cycle of the read, i.e. the closest to the other mate for fr or rf libraries when inner distance is positive) (default: read).  
  -s <scripts_dir_path>    Optional. Path of scripts directories. Must contains 'bam_to_mate_strand_bed.sh', 'plot_multiple_bed_depth.r' and 'bam_to_start_end_depth.py' (default: ./scripts/).  
```    
  
### Modules

#### bam_to_mate_strand_bed.sh

This script processes a BAM file to generate strand-specific and mate-specific depth files in a simple 3-column BED format (chromosome, position, depth).  
The script is compatible only with paired-end data and assumes the input BAM file is sorted by coordinates. No BAM index is used.
  
The script generates two main sets of outputs:  
  - Properly mapped reads, separated by strand (forward or reverse).  
  - Reads mapped without a mate, separated by strand and mate identity (first or second).  

Note about file naming for strand:  
  The strand designation (forward or reverse) is chosen arbitrarily and does not depend on the library type. However, since the same nomenclature is consistently used across files, this does not affect the analysis (in fact this denomination would correspond to a 'fr' library and it would be the opposite for an 'rf' library).  

This script is **not compatible with 'rr' or 'ff' libraries**.  

##### Arguments:
```
  -i <input_bam>         Recquired. Path to the input BAM file.  
  -o <output_dir>        Recquired. Directory to save the output files.  
  -d <delete_bam_files>  Optional. If set to "false", keeps intermediate BAM files (default: true, meaning BAM files are deleted).  
  -l <count_level>       Optional. Depth count level: consider the entire reads ['read'] or simply the last aligned position of the reads ['position'] (last in the context of the last cycle of the read, i.e. the closest to the other mate for fr or rf libraries when inner distance is positive) (default: read).  
  -s <scripts_dir_path>  Optional. Path of scripts directories, must contains 'bam_to_start_end_depth.py' (used if counting at position reslution instead of read resolution) (default: ./scripts/).  
```  
 
##### Usage:
  `bash process_bam.sh -i <input_bam> -o <output_dir> [option(s)]`  

##### Requirements:
 - `samtools` for BAM manipulation.  
 - `bedtools` for creating depth files in BED format (if `-l read`).  
 - `pysam` on `python3` (if `-l position`).  

#### plot_multiple_bed_depth.r

This script generates depth plots from BED files and saves them in a PDF. It takes multiple BED files as input, provided they reference the same genome and contain three columns (chromosome name, position, and depth). Each chromosome's depth data is represented as a line graph, with one plot per chromosome. Optional parameters allow customization of each plot's legend, color, and y-axis position.  
  
If opposite strand BED files are provided, these are displayed as a mirrored image on an inverted y-axis. This setup helps visualize complementary strand data alongside the primary strand.  

##### Arguments:
```
  --bed_files=<list_path>           Required. Comma-separated list of BED files with depth data (no spaces).  
  --output_pdf=<file>               Required. Path to the output PDF file where the plots will be saved.  
  --plot_title=<string>             Optional. Title of the plot (default: 'Depth').  
  --plot_legend=<list_string>       Optional. Comma-separated list of names for the each line plot (default basename of '--bed_files')  
  --plot_colors=<list_string>       Optional. Comma-separated list of color codes or names for each line plot.  
  --bed_files_opposite=<list_path>  Optional. Comma-separated list of BED files for opposite strand.  
  --cov_axes=<int_list>             Optional. Comma-separated list specifying the axis for each BED file (1 for primary, 2 for secondary).  
  --plot_alpha=<float>              Optional. Transparency of plot lines (default: 0.4).  
```  
##### Usage
  `plot_multiple_bed_depth.R --bed_files=<input/bed_files_path_1,input/bed_files_path_2> --output_pdf=<out/dir/depth.pdf> [options]`  

##### Requirements:
 - `ggplot2` (on `R`)  
