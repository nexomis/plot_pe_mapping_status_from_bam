#!/usr/bin/env Rscript

## Required packages
library(ggplot2)

## Define function

#' Plot Coverage from BED Files to PDF
#'
#' This function generates coverage plots from BED files and saves the result in a PDF. 
#' It reads coverage data from each BED file and plots them as line graphs with one graph by chromosme. 
#' You can provide multiple BED files from same reference.
#' If opposite strand files are provided, they are plotted on an inverted y-axis.
#'
#' @param bed_files A character vector containing the paths to the BED files with coverage data. 
#'                  Each file must have three columns: chromosome name, position, and coverage.
#' @param output_pdf The path to the output PDF file where the plots will be saved.
#' @param plot_title The title of the plot. For example, can be the sample name. Default is "Coverage".
#' @param plot_legend A character vector containing names for the plots. If NULL, the function automatically generates names from the 'bed_files' file basenames.
#' @param plot_colors A character vector of color codes or names for each plot. If NULL, colors are automatically assigned.
#' @param bed_files_opposite_strand A character vector containing the paths to the BED files for the opposite strand. If NULL, no opposite strand data is plotted.
#' @param cov_axes A numeric vector specifying which axis (1 for primary, 2 for secondary) each BED file (and opposite strand file) should be plotted on. 
#'                 If NULL, all files are plotted on the primary axis by default.
#' @param plot_alpha The transparency of the plot lines. Default is 0.5.
#'
#' @return This function does not return a value but generates and saves a PDF file with the coverage plots.
#' 
#' @details
#' This function plots coverage from simple BED files as input.
#' Each BED file must contain three columns: the chromosome name, position, and coverage.
#' Others formats (such as those with intervals) are not compatible.
#' Coverage must not be negative (conflict with the representation of the opposite strand)
#' 
#' It accepts multiple BED files through the bed_files parameter, provided all the files
#' contain the same chromosomes (to ensure they refer to the same genomic reference).
#' 
#' The function generates one plot per chromosome, with each BED file represented by a curve.
#' The name, color, and axis (primary or secondary) of each curve are determined
#' by the optional parameters plot_legend, plot_colors, and plot_axis.
#' 
#' Additional BED files can be provided through the bed_files_opposite_strand parameter.
#' These files must also contain the same chromosomes as those in bed_files.
#' The coverage from these files will be plotted as a mirror image on an inverted y-axis.
#'
#' @examples
#' # Example usage
#' plot_multiple_bed_coverage(
#'   bed_files = c("gloabal_pos.bed", "only_R1_map_pos.bed", "only_improperly_pos.bed"),
#'   bed_files_opposite_strand = c("gloabal_neg.bed", "only_R1_map_neg.bed", "only_improperly_neg.bed"),
#'   output_pdf = "coverage_plot_sample_XX.pdf",
#'   cov_axes = c(1, 2, 2),
#'   plot_legend = c("all", "R1", "improperly"),
#'   plot_colors = c("green", "red", "blue")
#' )
#' 
#' @export

plot_multiple_bed_coverage <- function(bed_files, output_pdf,
                                 plot_title = "Coverage",
                                 plot_legend = NULL, plot_colors = NULL,
                                 bed_files_oposite_strand = NULL,
                                 cov_axes = NULL,
                                 plot_alpha = 0.4) {
  
  # missing of missing parameters
  if (is.null(plot_legend)) {
    print("INFO: 'plot_legend' is automatically determined uniquely from the basename of 'bed_files' ('bed_files_oposite_strand' is not considered for this step)")
    plot_legend <- sapply(basename(bed_files), function(x) tools::file_path_sans_ext(x))
    plot_legend <- make.unique(plot_legend)
  }
  
  if (is.null(plot_colors)) {
    plot_colors <- scales::hue_pal()(length(bed_files))
  }
  
  # check conformity of some parameters
  if (length(bed_files) != length(plot_legend) ||
      length(bed_files) != length(plot_colors)) {
    stop("Error: 'bed_files', 'plot_colors', and 'plot_legend' must have the same length")
  }
  
  # check conformity of 'cov_axes'
  if (!is.null(cov_axes)) {
    if (length(cov_axes) != length(bed_files)) {
      stop("Error: 'cov_axes' must have the same length as 'bed_files'")
    }
    if (!all(cov_axes %in% c(1, 2))) {
      stop("Error: 'cov_axes' must contain only values '1' (for primary axis) or '2' (for secondary axis)")
    }
  } else {
    cov_axes <- rep(1, length(bed_files))
  }
  
  # import all bed files and set their color
  bed_list <- lapply(seq_along(bed_files), function(i) {
    data <- read.table(bed_files[i], header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE,
                       col.names = c("chr", "pos", "cov"))
    data$color <- plot_legend[i]
    data
  })
  
  # import opposite strand bed files if provided, set their color and inverse coverage
  if (!is.null(bed_files_oposite_strand)) {
    if (length(bed_files_oposite_strand) != length(bed_files)) {
      stop("Error: 'bed_files_oposite_strand' must have the same length as 'bed_files'")
    }
    
    bed_opposite_list <- lapply(seq_along(bed_files_oposite_strand), function(i) {
      data <- read.table(bed_files_oposite_strand[i], header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE,
                         col.names = c("chr", "pos", "cov"))
      data$cov <- -data$cov
      data$color <- plot_legend[i]
      data
    })
  } else {
    bed_opposite_list <- list()
  }
  
  # check that all input bed have the same chr
  chrs <- unique(bed_list[[1]]$chr)
  for (bed_data in bed_list) {
    if (!all(unique(bed_data$chr) %in% chrs) ||
        !all(chrs %in% unique(bed_data$chr))) {
      stop("Error: All 'bed_files' files must contain exactly the same chr")
    }
  }
  
  if (length(bed_opposite_list) > 0) {
    for (bed_data in bed_opposite_list) {
      if (!all(unique(bed_data$chr) %in% chrs) ||
          !all(chrs %in% unique(bed_data$chr))) {
        stop("Error: All 'bed_files_oposite_strand' must contain exactly the same chr as 'bed_files'")
      }
    }
  }
  
  pdf(output_pdf, width = 10, height = 6)
  
  # draw one plot by chr
  for (chr in chrs) {
    # ggplot structure
    p <- ggplot() +
      labs(title = plot_title,
           subtitle = paste0(chr,
                             "\n(negative values correspond to the coverage of the given file as the opposite strand.)")) +
      geom_hline(yintercept = 0, color = "grey") +
      scale_color_manual(name = NULL,
                         values = setNames(plot_colors, plot_legend)) +
      theme(legend.position = "bottom", legend.direction = "horizontal")
    
    primary_axis_data <- list()
    secondary_axis_data <- list()
    max_cov_primary <- 0
    max_cov_secondary <- 0
    
    # draw one line by input bed file
    for (i in seq_along(bed_files)) {
      chr_data <- bed_list[[i]][bed_list[[i]]$chr == chr, ]
      
      # separate chr_data by axis and update max cov of axis
      if (cov_axes[i] == 1) {
        primary_axis_data[[i]] <- chr_data
        max_cov_primary <- max(max_cov_primary, max(chr_data$cov, na.rm = TRUE))
      } else {
        secondary_axis_data[[i]] <- chr_data
        max_cov_secondary <- max(max_cov_secondary, max(abs(chr_data$cov), na.rm = TRUE))
      }
      
      # if 'bed_opposite_list': separate chr_data by axis and update max cov of axis
      if (length(bed_opposite_list) > 0) {
        chr_opposite_data <- bed_opposite_list[[i]][bed_opposite_list[[i]]$chr == chr, ]
        
        if (cov_axes[i] == 1) {
          primary_axis_data[[length(bed_files) + i]] <- chr_opposite_data
          max_cov_primary <- max(max_cov_primary, max(abs(chr_opposite_data$cov), na.rm = TRUE))
        } else {
          secondary_axis_data[[length(bed_files) + i]] <- chr_opposite_data
          max_cov_secondary <- max(max_cov_secondary, max(abs(chr_opposite_data$cov), na.rm = TRUE))
        }
      }
    }
    
    # draw data on primary_axis
    for (data in primary_axis_data) {
      if (!is.null(data)) {
        p <- p + geom_line(data = data,
                           aes(x = pos, y = cov, color = color),
                           linewidth = 1, alpha = plot_alpha)
      }
    }
    
    # secondary_axis
    if (length(secondary_axis_data) > 0) {
      # compute scalling factor between both axis
      if (max_cov_secondary > 0) {
        secondary_axis_factor <- max_cov_primary / max_cov_secondary
      } else {
        secondary_axis_factor <- 1
      }
      p <- p + scale_y_continuous(name = "Coverage", 
                                  sec.axis = sec_axis(~ . / secondary_axis_factor,
                                                      name = "Coverage\n(in dotted)"))
      # draw data on secondary_axis
      for (data in secondary_axis_data) {
        if (!is.null(data)) {
          p <- p + geom_line(data = data,
                             aes(x = pos, y = cov * secondary_axis_factor, color = color),
                             linewidth = 1, linetype = "dotted", alpha = plot_alpha)
        }
      }
    }
    
    print(p)
  }
  
  dev.off()
}


## Help
print_help <- function() {
  cat("
Usage: plot_multiple_bed_coverage.R --bed_files=<input/bed_files_path_1,input/bed_files_path_2> --output_pdf=<out/dir/coverage.pdf> [options]

Plot Coverage from BED Files to PDF
This script generates coverage plots from BED files and saves the result in a PDF file.
It accepts multiple BED files from the same reference and supports opposite strand files.
Each BED file must have three columns: chromosome name, position, and coverage.

Options:
  --bed_files           Required. Comma-separated list of BED files with coverage data (no spaces).
  --output_pdf          Required. Path to the output PDF file where the plots will be saved.
  --plot_title          Optional. Title of the plot (default: 'Coverage').
  --plot_legend         Optional. Comma-separated list of names for the each line plot (default basename of '--bed_files')
  --plot_colors         Optional. Comma-separated list of color codes or names for each line plot.
  --bed_files_opposite  Optional. Comma-separated list of BED files for opposite strand.
  --cov_axes            Optional. Comma-separated list specifying the axis for each BED file (1 for primary, 2 for secondary).
  --plot_alpha          Optional. Transparency of plot lines (default: 0.4).

Example:
  ./plot_multiple_bed_coverage.R --bed_files=gloabal_pos.bed,only_R1_map_pos.bed --output_pdf=coverage.pdf --plot_legend='All,R1' --plot_colors='green,red'
\n")
}


#### args
args <- commandArgs(trailingOnly = TRUE)

if ("--help" %in% args || length(args) == 0) {
  print_help()
  quit()
}


## extract args
parse_arg <- function(flag, args) {
  val <- args[grep(paste0("^", flag, "="), args)]
  if (length(val) > 0) {
    return(sub(paste0(flag, "="), "", val))
  }
  return(NULL)
}

bed_files <- parse_arg("--bed_files", args)
output_pdf <- parse_arg("--output_pdf", args)
plot_title <- parse_arg("--plot_title", args)
plot_legend <- parse_arg("--plot_legend", args)
plot_colors <- parse_arg("--plot_colors", args)
bed_files_opposite <- parse_arg("--bed_files_opposite", args)
cov_axes <- parse_arg("--cov_axes", args)
plot_alpha <- parse_arg("--plot_alpha", args)

bed_files <- if (!is.null(bed_files)) strsplit(bed_files, ",")[[1]] else NULL
plot_legend <- if (!is.null(plot_legend)) strsplit(plot_legend, ",")[[1]] else NULL
plot_colors <- if (!is.null(plot_colors)) strsplit(plot_colors, ",")[[1]] else NULL
bed_files_opposite <- if (!is.null(bed_files_opposite)) strsplit(bed_files_opposite, ",")[[1]] else NULL
cov_axes <- if (!is.null(cov_axes)) as.integer(strsplit(cov_axes, ",")[[1]]) else NULL
plot_alpha <- if (!is.null(plot_alpha)) as.numeric(plot_alpha) else 0.4


## check args
if (is.null(bed_files) || is.null(output_pdf)) {
  cat("Error: '--bed_files' and '--output_pdf' are required\n")
  quit("no", 1)
}

## Call plot_multiple_bed_coverage with args
plot_multiple_bed_coverage(
  bed_files = bed_files,
  output_pdf = output_pdf,
  plot_title = plot_title,
  plot_legend = plot_legend,
  plot_colors = plot_colors,
  bed_files_oposite_strand = bed_files_opposite,
  cov_axes = cov_axes,
  plot_alpha = plot_alpha
)


