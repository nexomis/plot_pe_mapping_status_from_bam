#!/usr/bin/env Rscript

## Required packages
library(ggplot2)
library(dplyr)
library(gggenomes)
library(patchwork)
library(cowplot)
library(plotly)
library(htmlwidgets)


## Define function

#' Generate Depth Plot from BED Files and Save as PDF/HTML
#'
#' This function creates depth plots from BED files, generating one plot per chromosome and saving it to a PDF.
#' It also includes options to save individual plots as HTML files using Plotly for interactive exploration.
#' Additionally, gene annotations from a GFF file can be overlaid on the plots.
#'
#' @param bed_files Recquired. A character vector of paths to BED files containing depth data (chromosome, position, and depth).
#' @param output_pdf_bn Recquired. Output file path for the combined PDF plot (whithout extension).
#' @param plot_title Optional. Title for the plot. Defaults to "Depth".
#' @param plot_legend Optional. names for each plot. If NULL, names are derived from the file names in `bed_files`.
#' @param plot_colors Optional. colors for each plot. If NULL, colors are automatically assigned.
#' @param bed_files_opposite_strand Optional. Paths to BED files for opposite strand data.
#' @param depth_axes Optional. Numeric vector indicating the axis for each BED file (1 for primary, 2 for secondary).
#'                   Defaults to primary (1) for all.
#' @param plot_alpha Optional. Transparency level of plot lines, default is 0.4.
#' @param annot_gff_file Optional. Path to input GFF file for annotation plotting (if provided).
#' @param annot_feat_id_regex Optional. Pattern (sub) to extract feature ID from the GFF attributes.
#' @param annot_feat_id_catch Optional. Replacement (sub) for the feature ID.
#' @param annot_interest_type Optional. Vector of feature types to include (default: c("mRNA", "CDS", "transcript")).
#' @param annot_name_on_plot Optional. Logical; if TRUE, displays feature names on the plot (default: FALSE).
#' @param relative_heigt_combined_plot Optional. Relative heights for the combined depth and annotation plot (default: c(3, 1)).
#' @param max_pb_by_A4_width Optional. Maximum plot width in pixels for automatic page sizing (default: 5000).
#' @param min_x_graduation Optional. Minimum x-axis graduation in base pairs (default: 200).
#' @param plotly_out_dir Optional. Directory to save individual chromosome plots as HTML (interactive Plotly plots) (if provided).
#'
#' @return This function does not return a value but generates and saves a PDF/HTML files with depth plots.
#' 
#' @details
#' This function plots depth from simple BED files as input.
#' Each BED file must contain three columns: the chromosome name, position, and depth. Others formats (such as those with intervals) are not compatible.
#' Depth must not be negative (conflict with the representation of the opposite strand)
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
#' The depth from these files will be plotted as a mirror image on an inverted y-axis.
#'
#' @examples
#' # Example usage
#' plot_multiple_bed_depth(
#'   bed_files = c("global_pos.bed", "only_R1_map_pos.bed", "only_improperly_pos.bed"),
#'   bed_files_opposite_strand = c("global_neg.bed", "only_R1_map_neg.bed", "only_improperly_neg.bed"),
#'   output_pdf_bn = "depth_plot_sample_XX",
#'   depth_axes = c(1, 2, 2),
#'   plot_legend = c("all", "R1", "improperly"),
#'   plot_colors = c("green", "red", "blue"),
#'   annot_gff_file = "annotations.gff",
#'   plotly_out_dir = "plotly_chromosome_plots/"
#' )
#' 
#' @export
plot_multiple_bed_depth <- function(bed_files,
                                    output_pdf_bn,
                                    plot_title = "Depth",
                                    plot_legend = NULL,
                                    plot_colors = NULL,
                                    bed_files_oposite_strand = NULL,
                                    depth_axes = NULL,
                                    plot_alpha = 0.4,
                                    annot_gff_file = NULL,
                                    annot_feat_id_regex = ".*Target=([^; ]+).*",
                                    annot_feat_id_catch = "\\1",
                                    annot_interest_type = c("mRNA", "CDS", "transcript"),
                                    annot_name_on_plot = FALSE,
                                    relative_heigt_combined_plot = c(3, 1),
                                    max_pb_by_A4_width=5000,
                                    min_x_graduation=200,
                                    plotly_out_dir = NULL) {
  
  # determine value for some missing parameters
  if (is.null(plot_legend)) {
    print("INFO: 'plot_legend' is automatically determined uniquely from the basename of 'bed_files' ('bed_files_oposite_strand' is not considered for this naming step)")
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
  
  # check conformity of 'depth_axes'
  if (!is.null(depth_axes)) {
    if (length(depth_axes) != length(bed_files)) {
      stop("Error: 'depth_axes' must have the same length as 'bed_files'")
    }
    if (!all(depth_axes %in% c(1, 2))) {
      stop("Error: 'depth_axes' must contain only values '1' (for primary axis) or '2' (for secondary axis)")
    }
  } else {
    depth_axes <- rep(1, length(bed_files))
  }
  
  # import all bed files and set their color
  bed_list <- lapply(seq_along(bed_files), function(i) {
    data <- read.table(bed_files[i], header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE,
                       col.names = c("chr", "pos", "depth"))
    data$color <- plot_legend[i]
    data
  })
  
  # import opposite strand bed files if provided, set their color and inverse depth
  if (!is.null(bed_files_oposite_strand)) {
    if (length(bed_files_oposite_strand) != length(bed_files)) {
      stop("Error: 'bed_files_oposite_strand' must have the same length as 'bed_files'")
    }
    
    bed_opposite_list <- lapply(seq_along(bed_files_oposite_strand), function(i) {
      data <- read.table(bed_files_oposite_strand[i], header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE,
                         col.names = c("chr", "pos", "depth"))
      data$depth <- -data$depth
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
  
  # pdf with page width relatif to max chr len (limit of max pb by width of A4 foramt)
  max_chr_len <- max(tapply(bed_list[[1]]$pos, bed_list[[1]]$chr, max))
  base_page_width = 297/25.4  # A4 in inches
  nb_page = max(1, max_chr_len/max_pb_by_A4_width)
  max_width_page = base_page_width*nb_page
  pdf(paste0(output_pdf_bn, ".pdf"), width = max_width_page, height = (210/25.4))
  
  # draw one plot by chr
  for (chr in chrs) {
    chr_len <- max(bed_list[[1]]$pos[bed_list[[1]]$chr == chr])
    
    ## depth plot

    p_depth <- plot_depth_from_multiple_bed_data_and_signle_targeted_chr(chr = chr,
                                                                         bed_list = bed_list,
                                                                         bed_opposite_list = bed_opposite_list,
                                                                         depth_axes = depth_axes,
                                                                         plot_alpha = plot_alpha)
    p_depth <- p_depth +
      labs(title = plot_title,
           subtitle = paste0(chr,
                             "\n(negative values correspond to the depth of the given file as the opposite strand.)")) +
      geom_hline(yintercept = 0, color = "grey") +
      scale_color_manual(name = NULL,
                         values = setNames(plot_colors, plot_legend)) +
      theme_bw() +
      theme(panel.border = element_blank(),
            plot.margin = margin(t = 1, r = 1, b = 0, l = 1, unit = "line"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0),
            axis.title.x = element_blank()) +
      scale_x_continuous(
        position = "top",
        limits = c(0, chr_len),
        breaks = seq(0, round_down_to_lower_power_of_ten(chr_len), by = min(round_down_to_lower_power_of_ten(chr_len)/25, min_x_graduation)),
        labels = scales::comma)
    
    # save depth plot on plotly format
    if (!is.null(plotly_out_dir)) {
      if (dir.exists(plotly_out_dir)) {
        stop(paste0("ERROR : Output directories for saving in plotly format (plotly_out_dir: '", plotly_out_dir,"') must not exist."))
      }
      dir.create(plotly_out_dir)
      saveWidget(ggplotly(p_depth),
                 paste0(plotly_out_dir, "/chr_", chr, ".html"))
    }
    
    # legend parameters are not supported by plotly (e.g 'legend.direction = "horizontal"')
    p_depth <- p_depth +
      theme(legend.direction = "horizontal", 
            legend.position = "top",
            legend.justification = c(0, 0))
    
    ## annot plot
    if (!is.null(annot_gff_file)) {
      p_annot <- plot_annot_from_single_gff_file(gff_file = annot_gff_file,
                                                 chr = chr,
                                                 feat_id_regex = annot_feat_id_regex,
                                                 feat_id_catch = annot_feat_id_catch,
                                                 interest_type = annot_interest_type,
                                                 plot_alpha = plot_alpha,
                                                 name_on_plot = annot_name_on_plot)
      
      p_annot <- p_annot + theme_bw() +
        theme(legend.direction = "horizontal", 
              legend.position = "bottom",
              legend.justification = c(0, 0),
              panel.border = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              panel.spacing = unit(0, "lines"),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              plot.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "line"),
              axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
        labs(y = "Annotation") + 
        scale_x_continuous(
          limits = c(0, chr_len),
          breaks = seq(0, chr_len, by = min(round_down_to_lower_power_of_ten(chr_len)/25, min_x_graduation)),
          labels = scales::comma
        )
      
      ## combine both plot and specify width following chr_len
      combined_plot <- p_depth / p_annot +
        plot_layout(heights = relative_heigt_combined_plot)
    } else {
      combined_plot <- p_depth +
        theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "line"))
    }
    width_ratio <- chr_len / max_chr_len
    aligned_plot <- ggdraw() +
      draw_plot(combined_plot, x = 0, y = 0, width = width_ratio, height = 1)
    
    print(aligned_plot)
  }
  
  dev.off()
  return(combined_plot)
}



#' Round the number one power down of ten lower.
#'
#' Round the number one power down of ten lower.. For exemple, usefull do determine breack_by of axis plots.
#' For exemple: 14390 -> 14000 /15830 -> 16000
#' 
#' @param n Numeric value to round.
#' @return The input number rounded to the nearest power of ten down by one power.
round_down_to_lower_power_of_ten <- function(n) {
  power_of_ten <- 10^(floor(log10(n)) - 1)
  return(round(n / power_of_ten) * power_of_ten)
}




#' Generate Depth Plot for a Single Chromosome from Multiple BED Files
#'
#' Creates a depth plot for a specified chromosome from multiple BED files, with support for both primary and secondary axes.
#'
#' @param p Base ggplot object to add depth data onto.
#' @param chr Target chromosome name to plot.
#' @param bed_list List of data frames, each containing depth data for a BED file.
#' @param bed_opposite_list Optional list of data frames for opposite-strand depth data.
#' @param depth_axes Numeric vector specifying the axis (1 or 2) for each BED file.
#' @param plot_alpha Transparency level of plot lines.
#' @return A ggplot object with the depth data plotted for the specified chromosome.
plot_depth_from_multiple_bed_data_and_signle_targeted_chr <- function(p = ggplot(),
                                                                      chr,
                                                                      bed_list,
                                                                      bed_opposite_list = NULL,
                                                                      depth_axes = NULL,
                                                                      plot_alpha = 0.4) {
  primary_axis_data <- list()
  secondary_axis_data <- list()
  max_depth_primary <- 0
  max_depth_secondary <- 0
  
  # draw one line by input bed file
  for (i in seq_along(bed_list)) {
    chr_data <- bed_list[[i]][bed_list[[i]]$chr == chr, ]
    
    # separate chr_data by axis and update max depth of axis
    if (depth_axes[i] == 1) {
      primary_axis_data[[i]] <- chr_data
      max_depth_primary <- max(max_depth_primary, max(chr_data$depth, na.rm = TRUE))
    } else {
      secondary_axis_data[[i]] <- chr_data
      max_depth_secondary <- max(max_depth_secondary, max(abs(chr_data$depth), na.rm = TRUE))
    }
    
    # if 'bed_opposite_list': separate chr_data by axis and update max depth of axis
    if (length(bed_opposite_list) > 0) {
      chr_opposite_data <- bed_opposite_list[[i]][bed_opposite_list[[i]]$chr == chr, ]
      
      if (depth_axes[i] == 1) {
        primary_axis_data[[length(bed_list) + i]] <- chr_opposite_data
        max_depth_primary <- max(max_depth_primary, max(abs(chr_opposite_data$depth), na.rm = TRUE))
      } else {
        secondary_axis_data[[length(bed_list) + i]] <- chr_opposite_data
        max_depth_secondary <- max(max_depth_secondary, max(abs(chr_opposite_data$depth), na.rm = TRUE))
      }
    }
  }
  
  # draw data on primary_axis
  for (data in primary_axis_data) {
    if (!is.null(data)) {
      p <- p + geom_line(data = data,
                         aes(x = pos, y = depth, color = color),
                         linewidth = 1, alpha = plot_alpha)
    }
  }
  
  # secondary_axis
  if (length(secondary_axis_data) > 0) {
    # compute scalling factor between both axis
    if (max_depth_secondary > 0) {
      secondary_axis_factor <- max(max_depth_primary, 30) / max(max_depth_secondary, 30)
    } else {
      secondary_axis_factor <- 1
    }
    p <- p + scale_y_continuous(name = "Depth", 
                                sec.axis = sec_axis(
                                  transform = ~ . / secondary_axis_factor,
                                  name = "Depth\n(in dotted)",
                                  breaks = seq(if (length(bed_opposite_list) == 0) 0 else -max(round_down_to_lower_power_of_ten(max_depth_secondary), 30),
                                               max(round_down_to_lower_power_of_ten(max_depth_secondary), 30),
                                               by=max(round_down_to_lower_power_of_ten(max_depth_secondary)/10, 5)),
                                  labels = scales::comma),
                                limits = c(if (length(bed_opposite_list) == 0) 0 else -max(round_down_to_lower_power_of_ten(max_depth_primary), 30),
                                           max(max_depth_primary, 30)),
                                breaks = seq(if (length(bed_opposite_list) == 0) 0 else -max(round_down_to_lower_power_of_ten(max_depth_primary), 30),
                                             max(max_depth_primary, 30),
                                             by = max(round_down_to_lower_power_of_ten(max_depth_primary)/10, 5)),
                                labels = scales::comma)
    # draw data on secondary_axis
    for (data in secondary_axis_data) {
      if (!is.null(data)) {
        p <- p + geom_line(data = data,
                           aes(x = pos, y = depth * secondary_axis_factor, color = color),
                           linewidth = 1, linetype = "dotted", alpha = plot_alpha)
      }
    }
  } else {
    p <- p + scale_y_continuous(limits = c(if (length(bed_opposite_list) == 0) 0 else -max(round_down_to_lower_power_of_ten(max_depth_primary), 30),
                                           max(round_down_to_lower_power_of_ten(max_depth_primary), 30)),
                                breaks = seq(-max(round_down_to_lower_power_of_ten(max_depth_primary), 30),
                                             max(round_down_to_lower_power_of_ten(max_depth_primary), 30),
                                             by = max(round_down_to_lower_power_of_ten(max_depth_primary)/10, 5)),
                                labels = scales::comma)
  }
  return(p)
}



#' Generate Gene Annotation Plot for a single chromosome from GFF File
#'
#' This function creates an annotation plot based on features in a GFF file for one signle targeted chromosome.
#'
#' @param gff_file Path to the GFF file containing gene features.
#' @param chr Target chromosome name to plot.
#' @param feat_id_regex Regex to extract feature ID from the attributes column.
#' @param feat_id_catch Replacement pattern for extracted ID.
#' @param interest_type Vector of feature types to include (e.g., "mRNA", "CDS").
#' @param name_on_plot Logical; if TRUE, displays feature names on the plot.
#' @param plot_alpha Transparency level of plot features.
#' @return A ggplot object representing gene annotations.
plot_annot_from_single_gff_file <- function(gff_file,
                                            chr,
                                            feat_id_regex = ".*Target=([^; ]+).*",
                                            feat_id_catch = "\\1",
                                            interest_type = c("mRNA", "CDS", "transcript"),
                                            plot_alpha = 0.4,
                                            name_on_plot = FALSE) {
  gff_data <- read.delim(gff_file, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  
  # rename and filter additional column and row following chromosome and feature type.
  colnames(gff_data)[1:9] <- c("seq_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  gff_clean <- gff_data[, 1:9] %>%
    filter(seq_id == chr & type %in% interest_type)
  
  # extract feat_id ID from the 'attributes' column using regex parameters
  gff_clean$feat_id <- sub(feat_id_regex, feat_id_catch, gff_clean$attributes)
  
  # plot
  p_annot <- gggenomes(genes = gff_clean) + 
    geom_gene(aes(fill = feat_id), position = "pile", alpha = plot_alpha,
              size = 2, shape = c(2, 2), rna_size = 1, stroke = 0.5) + 
    geom_seq_label()
  name_on_plot = FALSE
  if (name_on_plot){
    p_annot <- p_annot + geom_gene_note(aes(label = feat_id)) 
  }
  
  return(p_annot)
}



if (!interactive()) {
  
  ## Help
  print_help <- function() {
    cat("
  Usage: plot_multiple_bed_depth.R --bed_files=<input/bed_files_path_1,input/bed_files_path_2> --output_pdf_bn=<out/dir/depth> [options]
  
  Plot Depth from BED Files to PDF
  This script generates depth plots from BED files, saving them in a PDF. Optional HTML (Plotly) plots per chromosome can also be saved.
  Annotations from a GFF file can be added if provided.
  
  Options:
    --bed_files                    Required. Comma-separated list of BED files with depth data (no spaces).
    --output_pdf_bn                Required. Path to the output PDF file.
    --plot_title                   Optional. Title for the plot (default: 'Depth').
    --plot_legend                  Optional. Names for each line plot (default: file names in '--bed_files').
    --plot_colors                  Optional. Color codes or names for each plot.
    --bed_files_opposite           Optional. Comma-separated list of BED files for opposite strand data.
    --depth_axes                   Optional. Axis for each BED file (1 for primary (line), 2 for secondary (dotted)).
    --plot_alpha                   Optional. Transparency of plot lines (default: 0.4).
    --annot_gff_file               Optional. Path to GFF file for annotation.
    --annot_feat_id_regex          Optional. Regex to extract feature ID from GFF attributes (default: '.*Target=([^; ]+).*').
    --annot_feat_id_catch          Optional. Regex replacement pattern for feature ID (default: '\\1').
    --annot_interest_type          Optional. Feature types to include from GFF file (default: c('mRNA', 'CDS', 'transcript')).
    --annot_name_on_plot           Optional. Display feature names on plot (default: FALSE).
    --relative_heigt_combined_plot Optional. Relative heights of depth and annotation plots (default: c(3, 1)).
    --max_pb_by_A4_width           Optional. Maximum page width in pixels for auto-sizing (default: 5000).
    --min_x_graduation             Optional. Minimum x-axis graduation in bp (default: 200).
    --plotly_out_dir               Optional. Directory to save individual chromosome plots as HTML (interactive).
  
  Example:
    ./plot_multiple_bed_depth.R --bed_files=gloabal_pos.bed,only_R1_map_pos.bed --output_pdf_bn=depth.pdf --plot_legend='All,R1' --plot_colors='green,red' --plotly_out_dir=plotly_chromosome_plots
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
  output_pdf_bn <- parse_arg("--output_pdf_bn", args)
  plot_title <- parse_arg("--plot_title", args)
  plot_legend <- parse_arg("--plot_legend", args)
  plot_colors <- parse_arg("--plot_colors", args)
  bed_files_opposite <- parse_arg("--bed_files_opposite", args)
  depth_axes <- parse_arg("--depth_axes", args)
  plot_alpha <- parse_arg("--plot_alpha", args)
  annot_gff_file <- parse_arg("--annot_gff_file", args)
  annot_feat_id_regex <- parse_arg("--annot_feat_id_regex", args)
  annot_feat_id_catch <- parse_arg("--annot_feat_id_catch", args)
  annot_interest_type <- parse_arg("--annot_interest_type", args)
  annot_name_on_plot <- parse_arg("--annot_name_on_plot", args)
  relative_heigt_combined_plot <- parse_arg("--relative_heigt_combined_plot", args)
  max_pb_by_A4_width <- parse_arg("--max_pb_by_A4_width", args)
  min_x_graduation <- parse_arg("--min_x_graduation", args)
  plotly_out_dir <- parse_arg("--plotly_out_dir", args) 
  
  # format args and default value
  bed_files <- if (!is.null(bed_files)) strsplit(bed_files, ",")[[1]] else NULL
  plot_legend <- if (!is.null(plot_legend)) strsplit(plot_legend, ",")[[1]] else NULL
  plot_colors <- if (!is.null(plot_colors)) strsplit(plot_colors, ",")[[1]] else NULL
  bed_files_opposite <- if (!is.null(bed_files_opposite)) strsplit(bed_files_opposite, ",")[[1]] else NULL
  depth_axes <- if (!is.null(depth_axes)) as.integer(strsplit(depth_axes, ",")[[1]]) else NULL
  plot_alpha <- if (!is.null(plot_alpha)) as.numeric(plot_alpha) else 0.4
  plot_title <- if (!is.null(plot_title)) plot_title else "Depth"
  annot_gff_file <- if (!is.null(annot_gff_file)) annot_gff_file else NULL
  annot_feat_id_regex <- if (!is.null(annot_feat_id_regex)) annot_feat_id_regex else ".*Target=([^; ]+).*"
  annot_feat_id_catch <- if (!is.null(annot_feat_id_catch)) annot_feat_id_catch else "\\1"
  annot_interest_type <- if (!is.null(annot_interest_type)) strsplit(annot_interest_type, ",")[[1]] else c("mRNA", "CDS", "transcript")
  annot_name_on_plot <- if (!is.null(annot_name_on_plot)) as.logical(annot_name_on_plot) else FALSE
  relative_heigt_combined_plot <- if (!is.null(relative_heigt_combined_plot)) as.numeric(strsplit(relative_heigt_combined_plot, ",")[[1]]) else c(3, 1)
  max_pb_by_A4_width <- if (!is.null(max_pb_by_A4_width)) as.numeric(max_pb_by_A4_width) else 5000
  min_x_graduation <- if (!is.null(min_x_graduation)) as.numeric(min_x_graduation) else 200
  plotly_out_dir <- if (!is.null(plotly_out_dir)) plotly_out_dir else NULL
  
  ## check args
  if (is.null(bed_files) || is.null(output_pdf_bn)) {
    cat("Error: '--bed_files' and '--output_pdf_bn' are required (use '=' to asign value to parameter)\n")
    quit("no", 1)
  }
  
  ## Call plot_multiple_bed_depth with args

  plot_multiple_bed_depth(
    bed_files = bed_files,
    output_pdf_bn = output_pdf_bn,
    plot_title = plot_title,
    plot_legend = plot_legend,
    plot_colors = plot_colors,
    bed_files_oposite_strand = bed_files_opposite,
    depth_axes = depth_axes,
    plot_alpha = plot_alpha,
    annot_gff_file = annot_gff_file,
    annot_feat_id_regex = annot_feat_id_regex,
    annot_feat_id_catch = annot_feat_id_catch,
    annot_interest_type = annot_interest_type,
    annot_name_on_plot = annot_name_on_plot,
    relative_heigt_combined_plot = relative_heigt_combined_plot,
    max_pb_by_A4_width = max_pb_by_A4_width,
    min_x_graduation = min_x_graduation,
    plotly_out_dir =plotly_out_dir
  )
}
