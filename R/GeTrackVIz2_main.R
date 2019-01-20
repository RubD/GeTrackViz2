#' @title RNA bedgraph plotter
#' @description This function plots RNA bedgraph data
#' @param rna_bdg bedgraph file for RNA expression data
#' @param name_rna_bdg name of sample
#' @param color_rna_bdg color of track
#' @param mychr chromosome
#' @param start_loc start location on chromomse
#' @param end_loc end location on chromosome
#' @param binsize size of bin
#' @param breaks_wanted number of breaks on x-axis
#' @param y_title title for y-axis
#' @param chrom_col column name for chromosomes
#' @param start_col column name for start
#' @param end_col column name for end
#' @param counts_col column name with counts
#' @param strand_col column name with strand info, default = NA
#' @param min_strand_sign minus strand sign, default = '-'
#' @param plus_strand_sign plus strand sign, default = '+'
#' @param make_minus_strand_negative reverse strand
#' @param forced_y_min hard minimum y-axis value
#' @param forced_y_max hard maximum y-axis value
#' @param axis_line_size line size of x-axis
#' @param marginvec_mm margin around plot, default = c(0,0,0,0)
#' @param use_barplot plot track as barplot, default = F
#' @param alpha_fill alpha for lineplot, default = 1
#' @param reverse_strand boolean: flip strand x-axis direction
#' @param show_labels boolean: show x-axis lables
#' @param print_plot boolean: print individual plot
#' @import ggplot2 data.table
#' @return ggplot object
#' @keywords RNA, bedgraph
#' @examples
#'     plot_RNA_bedgraph(rna_bdg)
#'
#'
#' @export
plot_RNA_bedgraph <- function(rna_bdg, name_rna_bdg, color_rna_bdg, mychr, start_loc, end_loc, binsize = 100,
                              breaks_wanted = 4, y_title = 'RPKM',
                              chrom_col = 'chr', start_col = 'start', end_col = 'end', counts_col = 'RPKM', strand_col = NA,
                              min_strand_sign = '-', plus_strand_sign = '+',
                              make_minus_strand_negative = T,
                              forced_y_min = NULL, forced_y_max = NULL,
                              axis_line_size = 0.2, marginvec_mm = c(0,0,0,0),
                              use_barplot = F, alpha_fill = 0.7,
                              reverse_strand = FALSE,
                              show_labels = F, print_plot = F) {



  # get full range of bins
  full_range <- data.table(chr = mychr, start = seq(start_loc, end_loc, by = binsize))
  full_range[, end := start + (binsize-1)]
  full_range[, id := paste0('id_', 1:nrow(full_range))]
  #print(full_range)


  # if there is a strand column
  if(!is.na(strand_col)) {

    cat('\n stranded data \n')

    rna_bdg <- rna_bdg[, c(chrom_col, start_col, end_col, counts_col, strand_col), with = F]
    colnames(rna_bdg) <- c('chr', 'start', 'end', 'counts', 'strand')

    bdg_subset_plus <- rna_bdg[strand == plus_strand_sign]
    bdg_subset_plus <- bdg_subset_plus[chr == mychr & start >= start_loc & end <= end_loc]
    bdg_subset_plus[, name := paste0('name_', 1:nrow(bdg_subset_plus))]
    bdg_subset_plus <- bdg_subset_plus[, .(chr, start, end, name, counts)]


    bdg_subset_minus <- rna_bdg[strand == min_strand_sign]
    bdg_subset_minus <- bdg_subset_minus[chr == mychr & start >= start_loc & end <= end_loc]
    bdg_subset_minus[, name := paste0('name_', 1:nrow(bdg_subset_minus))]
    bdg_subset_minus <- bdg_subset_minus[, .(chr, start, end, name, counts)]



    # calculate overlap
    overlap_plus <- bedOverlap(bdg_subset_plus, full_range, minOverlap = 2)
    overlap_plus <- overlap_plus[, .(chr, start, end, name, counts)]; setkey(overlap_plus, chr, start, end)
    overlap_plus <- overlap_plus[, sum(counts), by = c('chr', 'start', 'end', 'name')]

    overlap_minus <- bedOverlap(bdg_subset_minus, full_range, minOverlap = 2)
    overlap_minus <- overlap_minus[, .(chr, start, end, name, counts)]; setkey(overlap_minus, chr, start, end)
    overlap_minus <- overlap_minus[, sum(counts), by = c('chr', 'start', 'end', 'name')]



    # PLUS #
    # add regions with no overlap = 0
    full_range_zero_plus <- full_range[!id %in% overlap_plus$name]
    full_range_zero_plus[, V1 := 0]
    setnames(full_range_zero_plus, 'id', 'name')

    overlap_plus <- rbind(overlap_plus, full_range_zero_plus)
    setkey(overlap_plus, chr, start, end)

    overlap_plus[, name := factor(name, levels = name)]
    overlap_plus[, num_region := 1:nrow(overlap_plus)]

    # give name
    overlap_plus[, condition := name_rna_bdg]
    # region center
    overlap_plus[, region_center := round(mean(c(start, end)), digits = 0), by = 1:nrow(overlap_plus)]


    # MINUS #
    # add regions with no overlap = 0
    full_range_zero_minus <- full_range[!id %in% overlap_minus$name]
    full_range_zero_minus[, V1 := 0]
    setnames(full_range_zero_minus, 'id', 'name')

    overlap_minus <- rbind(overlap_minus, full_range_zero_minus)
    setkey(overlap_minus, chr, start, end)

    overlap_minus[, name := factor(name, levels = name)]
    overlap_minus[, num_region := 1:nrow(overlap_minus)]

    # give name
    overlap_minus[, condition := name_rna_bdg]
    # region center
    overlap_minus[, region_center := round(mean(c(start, end)), digits = 0), by = 1:nrow(overlap_minus)]

    if(make_minus_strand_negative == T) overlap_minus[, V1 := -V1]



    # create y title
    full_title = paste0(name_rna_bdg, '\n', y_title)

    # create breaks
    calculated_breaks <- seq(start_loc, end_loc, length.out = breaks_wanted)

    if(show_labels == T) {
      mylabels <- paste0(round(calculated_breaks/1000, digits = 0), 'kb')
    } else {
      mylabels <- rep('', length(calculated_breaks))
    }


    # set y limits if required
    my_y_min <- ifelse(is.null(forced_y_min), min(overlap_minus[['V1']]), forced_y_min)
    my_y_max <- ifelse(is.null(forced_y_max), max(overlap_plus[['V1']]), forced_y_max)

    # breaks for y-axis
    y_breaks <- seq(my_y_min, my_y_max, length.out = 3)
    my_y_labels <- round(y_breaks, digits = 1)


    ## Reverse strand ##
    if(reverse_strand == TRUE) {
      my_y_min = -my_y_min
      my_y_max = -my_y_max
      overlap_plus[, V1 := -V1]
      overlap_minus[, V1 := -V1]
      y_breaks = -y_breaks
      my_y_labels <- round(y_breaks, digits = 1)
    }



    plus_color = color_rna_bdg[[1]]
    minus_color = color_rna_bdg[[2]]


    pl <- ggplot()
    # barplot vs lineplot
    if(use_barplot == TRUE) {

      if(reverse_strand == TRUE) {
        pl <- pl + geom_bar(data = overlap_plus, aes(x = rev(region_center), y = V1), stat = 'identity', width = 1, color = plus_color)
        pl <- pl + geom_bar(data = overlap_minus, aes(x = rev(region_center), y = V1), stat = 'identity', width = 1, color = minus_color)

      } else {
        pl <- pl + geom_bar(data = overlap_plus, aes(x = region_center, y = V1), stat = 'identity', width = 1, color = plus_color)
        pl <- pl + geom_bar(data = overlap_minus, aes(x = region_center, y = V1), stat = 'identity', width = 1, color = minus_color)

      }

    } else {

      if(reverse_strand == TRUE) {
        pl <- pl + geom_line(data = overlap_plus, aes(x = rev(region_center), y = V1), stat = 'identity', color = plus_color)
        pl <- pl + geom_ribbon(data = overlap_plus, aes(x = rev(region_center), ymax = V1, ymin = 0), fill = plus_color, alpha = alpha_fill)
        pl <- pl + geom_line(data = overlap_minus, aes(x = rev(region_center), y = V1), stat = 'identity', color = minus_color)
        pl <- pl + geom_ribbon(data = overlap_minus, aes(x = rev(region_center), ymax = 0, ymin = V1), fill = minus_color, alpha = alpha_fill)

      } else {
        pl <- pl + geom_line(data = overlap_plus, aes(x = region_center, y = V1), stat = 'identity', color = plus_color)
        pl <- pl + geom_ribbon(data = overlap_plus, aes(x = region_center, ymax = V1, ymin = 0), fill = plus_color, alpha = alpha_fill)
        pl <- pl + geom_line(data = overlap_minus, aes(x = region_center, y = V1), stat = 'identity', color = minus_color)
        pl <- pl + geom_ribbon(data = overlap_minus, aes(x = region_center, ymax = 0, ymin = V1), fill = minus_color, alpha = alpha_fill)
      }

    }

    pl <- pl + geom_hline(yintercept = 0, color = 'black', size = 0.5)
    pl <- pl + theme_bw() + theme(panel.background = element_blank(),
                                  panel.border = element_blank(),
                                  axis.line = element_line(color = 'black', size = 0.2),
                                  panel.grid = element_blank(),
                                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                                  plot.margin=unit(marginvec_mm,"mm"))
    pl <- pl + scale_x_continuous(expand = c(0, 0), limits = c(start_loc, end_loc),
                                  breaks = calculated_breaks, labels = mylabels)
    pl <- pl + scale_y_continuous(breaks = y_breaks, labels = my_y_labels)
    pl <- pl + coord_cartesian(ylim=c(my_y_min, my_y_max))
    pl <- pl + labs(x = NULL, y = full_title)


    if(print_plot == T) print(pl)

    return(pl)




  }

  else if(is.na(strand_col)) {

    cat('\n only one strand data \n')

    # no strand difference
    # TO DO

    rna_bdg <- rna_bdg[, c(chrom_col, start_col, end_col, counts_col), with = F]
    colnames(rna_bdg) <- c('chr', 'start', 'end', 'counts')

    bdg_subset_plus <- rna_bdg[chr == mychr & start >= start_loc & end <= end_loc]
    bdg_subset_plus[, name := paste0('name_', 1:nrow(bdg_subset_plus))]
    bdg_subset_plus <- bdg_subset_plus[, .(chr, start, end, name, counts)]


    # calculate overlap
    overlap_plus <- bedOverlap(bdg_subset_plus, full_range, minOverlap = 2)
    overlap_plus <- overlap_plus[, .(chr, start, end, name, counts)]; setkey(overlap_plus, chr, start, end)
    overlap_plus <- overlap_plus[, sum(counts), by = c('chr', 'start', 'end', 'name')]



    # PLUS #
    # add regions with no overlap = 0
    full_range_zero_plus <- full_range[!id %in% overlap_plus$name]
    full_range_zero_plus[, V1 := 0]
    setnames(full_range_zero_plus, 'id', 'name')

    overlap_plus <- rbind(overlap_plus, full_range_zero_plus)
    setkey(overlap_plus, chr, start, end)

    overlap_plus[, name := factor(name, levels = name)]
    overlap_plus[, num_region := 1:nrow(overlap_plus)]

    # give name
    overlap_plus[, condition := name_rna_bdg]
    # region center
    overlap_plus[, region_center := round(mean(c(start, end)), digits = 0), by = 1:nrow(overlap_plus)]


    if(make_minus_strand_negative == T) overlap_plus[, V1 := -V1]



    # create y title
    full_title = paste0(name_rna_bdg, '\n', y_title)

    # create breaks
    calculated_breaks <- seq(start_loc, end_loc, length.out = breaks_wanted)

    if(show_labels == T) {
      mylabels <- paste0(round(calculated_breaks/1000, digits = 0), 'kb')
    } else {
      mylabels <- rep('', length(calculated_breaks))
    }


    # set y limits if required
    my_y_min <- ifelse(is.null(forced_y_min), min(overlap_plus[['V1']]), forced_y_min)
    my_y_max <- ifelse(is.null(forced_y_max), max(overlap_plus[['V1']]), forced_y_max)

    # breaks for y-axis
    y_breaks <- seq(my_y_min, my_y_max, length.out = 3)
    my_y_labels <- round(y_breaks, digits = 1)

    plus_color = color_rna_bdg[[1]]


    pl <- ggplot()
    # choose barplot vs lineplot
    if(use_barplot == TRUE) {

      if(reverse_strand == TRUE) {
        pl <- pl + geom_bar(data = overlap_plus, aes(x = rev(region_center), y = V1), stat = 'identity', width = 1, color = plus_color)
      } else {
        pl <- pl + geom_bar(data = overlap_plus, aes(x = region_center, y = V1), stat = 'identity', width = 1, color = plus_color)
      }

    } else {

      if(reverse_strand == TRUE) {
        pl <- pl + geom_line(data = overlap_plus, aes(x = rev(region_center), y = V1), stat = 'identity', color = plus_color)
        pl <- pl + geom_ribbon(data = overlap_plus, aes(x = rev(region_center), ymax = V1, ymin = 0), fill = plus_color, alpha = alpha_fill)
      } else {
        pl <- pl + geom_line(data = overlap_plus, aes(x = region_center, y = V1), stat = 'identity', color = plus_color)
        pl <- pl + geom_ribbon(data = overlap_plus, aes(x = region_center, ymax = V1, ymin = 0), fill = plus_color, alpha = alpha_fill)
      }

    }

    pl <- pl + geom_hline(yintercept = 0, color = 'black', size = 0.5)
    pl <- pl + theme_bw() + theme(panel.background = element_blank(),
                                  panel.border = element_blank(),
                                  axis.line = element_line(color = 'black', size = 0.2),
                                  panel.grid = element_blank(),
                                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                                  plot.margin=unit(marginvec_mm,"mm"))
    pl <- pl + scale_x_continuous(expand = c(0, 0), limits = c(start_loc, end_loc),
                                  breaks = calculated_breaks, labels = mylabels)
    pl <- pl + scale_y_continuous(breaks = y_breaks, labels = my_y_labels)
    pl <- pl + coord_cartesian(ylim=c(my_y_min, my_y_max))
    pl <- pl + labs(x = NULL, y = full_title)


    if(print_plot == T) print(pl)

    return(pl)




  }


}




#' @title DNA bedgraph plotter
#' @description This function plots DNA bedgraph data
#' @param bdg bedgraph file for DNA coverage data
#' @param name_bdg name of sample
#' @param color_bdg color of track
#' @param mychr chromosome
#' @param start_loc start location on chromomse
#' @param end_loc end location on chromosome
#' @param binsize size of bin
#' @param breaks_wanted number of breaks on x-axis
#' @param y_title title for y-axis
#' @param chrom_col column name for chromosomes
#' @param start_col column name for start
#' @param end_col column name for end
#' @param counts_col column name with counts
#' @param forced_y_min hard minimum y-axis value
#' @param forced_y_max hard maximum y-axis value
#' @param axis_line_size line size of x-axis
#' @param marginvec_mm margin around plot, default = c(0,0,0,0)
#' @param use_barplot plot track as barplot, default = F
#' @param alpha_fill alpha for lineplot, default = 1
#' @param reverse_strand boolean: flip strand x-axis direction
#' @param show_labels boolean: show x-axis lables
#' @param print_plot boolean: print individual plot
#' @import ggplot2 data.table
#' @return ggplot object
#' @keywords DNA, coverage, bedgraph
#' @examples
#'     plot_bedgraph(bdg)
#'
#'
#' @export
plot_bedgraph <- function(bdg, name_bdg, color_bdg,
                          mychr, start_loc, end_loc, binsize = 100,
                          breaks_wanted = 4, y_title = 'RPM/bp',
                          chrom_col = 'chr', start_col = 'start', end_col = 'end', counts_col = 'RPKM',
                          forced_y_min = NULL, forced_y_max = NULL,
                          axis_line_size = 0.2, marginvec_mm = c(0,0,0,0),
                          use_barplot = F, alpha_fill = 1,
                          reverse_strand = FALSE,
                          show_labels = F, print_plot = F) {

  # subset of bedgraph
  bdg_subset <- bdg[, c(chrom_col, start_col, end_col, counts_col), with = F]
  colnames(bdg_subset) <- c('chr', 'start', 'end', 'counts')
  bdg_subset <- bdg_subset[chr == mychr & start >= start_loc & end <= end_loc]
  bdg_subset[, name := paste0('name_', 1:nrow(bdg_subset))]
  bdg_subset <- bdg_subset[, .(chr, start, end, name, counts)]
  #print(bdg_subset)

  # get full range of bins
  full_range <- data.table(chr = mychr, start = seq(start_loc, end_loc, by = binsize))
  full_range[, end := start + (binsize-1)]
  full_range[, id := paste0('id_', 1:nrow(full_range))]
  #print(full_range)

  # calculate overlap
  overlap <- bedOverlap(bdg_subset, full_range, minOverlap = 2)
  overlap <- overlap[, .(chr, start, end, name, counts)]; setkey(overlap, chr, start, end)
  overlap <- overlap[, sum(counts), by = c('chr', 'start', 'end', 'name')]


  # add regions with no overlap = 0
  full_range_zero <- full_range[!id %in% overlap$name]
  full_range_zero[, V1 := 0]
  setnames(full_range_zero, 'id', 'name')

  overlap <- rbind(overlap, full_range_zero)
  setkey(overlap, chr, start, end)

  overlap[, name := factor(name, levels = name)]
  overlap[, num_region := 1:nrow(overlap)]

  # give name
  overlap[, condition := name_bdg]

  # region center
  overlap[, region_center := round(mean(c(start, end)), digits = 0), by = 1:nrow(overlap)]

  # create y title
  full_title = paste0(name_bdg, '\n', y_title)

  # create breaks
  calculated_breaks <- seq(start_loc, end_loc, length.out = breaks_wanted)

  if(show_labels == T) {
    mylabels <- paste0(round(calculated_breaks/1000, digits = 0), 'kb')
  } else {
    mylabels <- rep('', length(calculated_breaks))
  }


  # set y limits if required
  my_y_min <- ifelse(is.null(forced_y_min), 0, forced_y_min)
  my_y_max <- ifelse(is.null(forced_y_max), max(overlap[['V1']]), forced_y_max)

  # breaks for y-axis
  y_breaks <- seq(my_y_min, my_y_max, length.out = 3)
  my_y_labels <- round(y_breaks, digits = 0)

  # create plot
  pl <- ggplot()


  ## TEST lineplot for DNA bedgraph data ##
  # choose barplot vs lineplot
  if(use_barplot == TRUE) {

    if(reverse_strand == TRUE) {
      pl <- pl + geom_bar(data = overlap, aes(x = rev(region_center), y = V1), stat = 'identity', width = 1, color = color_bdg)
    } else {
      pl <- pl + geom_bar(data = overlap, aes(x = region_center, y = V1), stat = 'identity', width = 1, color = color_bdg)
    }

  } else {

    if(reverse_strand == TRUE) {
      pl <- pl + geom_line(data = overlap, aes(x = rev(region_center), y = V1), stat = 'identity', color = color_bdg)
      pl <- pl + geom_ribbon(data = overlap, aes(x = rev(region_center), ymax = V1, ymin = 0), fill = color_bdg, alpha = alpha_fill)
    } else {
      pl <- pl + geom_line(data = overlap, aes(x = region_center, y = V1), stat = 'identity', color = color_bdg)
      pl <- pl + geom_ribbon(data = overlap, aes(x = region_center, ymax = V1, ymin = 0), fill = color_bdg, alpha = alpha_fill)
    }

  }

  ## END TEST ##

  #if(reverse_strand == TRUE) {
  #  pl <- pl + geom_bar(data = overlap, aes(x = rev(region_center), y = V1), stat = 'identity', width = 1, color = color_bdg)
  #} else {
  #  pl <- pl + geom_bar(data = overlap, aes(x = region_center, y = V1), stat = 'identity', width = 1, color = color_bdg)
  #}

  pl <- pl + theme_bw() + theme(panel.background = element_blank(),
                                panel.border = element_blank(),
                                axis.line = element_line(color = 'black', size = axis_line_size),
                                panel.grid = element_blank(),
                                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                                plot.margin=unit(marginvec_mm,"mm"))
  pl <- pl + scale_x_continuous(expand = c(0, 0), limits = c(start_loc, end_loc),
                                breaks = calculated_breaks, labels = mylabels)
  pl <- pl + scale_y_continuous(breaks = y_breaks, labels = my_y_labels)
  pl <- pl + coord_cartesian(ylim=c(my_y_min, my_y_max))
  pl <- pl + labs(x = NULL, y = full_title)

  if(print_plot == T) print(pl)

  return(pl)

}




#' @title Deprecated: loops or interactions plotter
#' @description This function plots data in bedpe format
#' @param bedpe bedpe file for DNA loops or interaction data
#' @param name_bedpe name of sample
#' @param color_bedpe color of track
#' @param mychr chromosome
#' @param start_loc start location on chromomse
#' @param end_loc end location on chromosome
#' @param breaks_wanted number of breaks on x-axis
#' @param y_title title for y-axis
#' @param chrom_col column name for chromosomes
#' @param start_col column name for start
#' @param end_col column name for end
#' @param show_partial_overlap include loops that do not entirely overlap the coordinates
#' @param axis_line_size line size of x-axis
#' @param marginvec_mm margin around plot, default = c(0,0,0,0)
#' @param show_labels boolean: show x-axis lables
#' @param print_plot boolean: print individual plot
#' @import ggplot2 data.table
#' @return ggplot object
#' @keywords loops, bedpe
#' @examples
#'     plot_loops_old(bedpe)
#'
#'
#' @export
plot_loops_old <- function(bedpe, name_bedpe, color_bedpe,
                       mychr, start_loc, end_loc,
                       breaks_wanted = 4, y_title = 'prob',
                       chrom_col = 'chr', start_col = 'start', end_col = 'end',
                       show_partial_overlap = T,
                       axis_line_size = 0.2, marginvec_mm = c(0,0,0,0),
                       reverse_strand = FALSE,
                       show_labels = F, print_plot = F) {


  # bedpe configuration is currently strict
  colnames(bedpe) <- paste0('V', 1:ncol(bedpe))

  # get start and end locations of anchor regions
  bedpe[, start_left_anchor := V2]
  bedpe[, end_left_anchor := V2 + as.numeric(strsplit(V11, split = ',')[[1]][1]), by = 1:nrow(bedpe)]

  bedpe[, start_right_anchor := V3 - as.numeric(strsplit(V11, split = ',')[[1]][2]), by = 1:nrow(bedpe)]
  bedpe[, end_right_anchor := V3]


  # take subset
  bedpe_subset <- bedpe[V1 == mychr & V2 >= start_loc & V3 <= end_loc]

  if(show_partial_overlap == TRUE) {

    # only the start of loop is within selected genomic region
    bedpe_subset_startOnly <- bedpe[V1 == mychr & V2 >= start_loc & V2 <= end_loc]
    bedpe_subset_startOnly <- bedpe_subset_startOnly[!V4 %in% bedpe_subset$V4]
    bedpe_subset_startOnly[, V3 := end_loc]
    bedpe_subset_startOnly[, end_right_anchor := ifelse(start_right_anchor <= end_loc, end_loc, end_right_anchor)]

    # only the end of loop is within selected genomic region
    bedpe_subset_endOnly <- bedpe[V1 == mychr & V3 >= start_loc & V3 <= end_loc]
    bedpe_subset_endOnly <- bedpe_subset_endOnly[!V4 %in% bedpe_subset$V4]
    bedpe_subset_endOnly[, V2 := start_loc]
    bedpe_subset_endOnly[, start_left_anchor := ifelse(end_left_anchor >= start_loc, start_loc, start_left_anchor)]

    # merge together
    bedpe_all_overlapping_subset <- do.call('rbind', list(bedpe_subset, bedpe_subset_startOnly, bedpe_subset_endOnly))
    bedpe_subset <- bedpe_all_overlapping_subset
  }




  # give y-values
  bedpe_subset[, segm_y_values := 1:nrow(bedpe_subset)]


  # create breaks
  calculated_breaks <- seq(start_loc, end_loc, length.out = breaks_wanted)

  if(show_labels == T) {
    mylabels <- paste0(round(calculated_breaks/1000, digits = 0), 'kb')
  } else {
    mylabels <- rep('', length(calculated_breaks))
  }

  # maximum y
  max_y = max(bedpe_subset$segm_y_values)+1

  # create y title
  full_title = paste0(name_bedpe, '\n', y_title)


  int_pl <- ggplot()
  int_pl <- int_pl + geom_segment(data = bedpe_subset,
                                  aes(x = end_left_anchor, y = segm_y_values, xend = start_right_anchor, yend = segm_y_values),
                                  color = color_bedpe)
  int_pl <- int_pl + geom_segment(data = bedpe_subset,
                                  aes(x = start_left_anchor, y = segm_y_values, xend = end_left_anchor, yend = segm_y_values),
                                  size = 1.2, color = color_bedpe)
  int_pl <- int_pl + geom_segment(data = bedpe_subset,
                                  aes(x = start_right_anchor, y = segm_y_values, xend = end_right_anchor, yend = segm_y_values),
                                  size = 1.2, color = color_bedpe)

  int_pl <- int_pl + scale_x_continuous(expand = c(0, 0), limits = c(start_loc, end_loc),
                                        breaks = calculated_breaks, labels = mylabels)
  int_pl <- int_pl + ylim(c(0,max_y))
  int_pl <- int_pl + theme_bw() + theme(panel.background = element_blank(),
                                        panel.border = element_blank(),
                                        axis.line = element_line(color = 'black', size = axis_line_size),
                                        panel.grid = element_blank(),
                                        axis.text.x = element_text(angle = 45, hjust = 45, vjust = 45),
                                        axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        plot.margin=unit(marginvec_mm,"mm"))
  int_pl <- int_pl + labs(list(x = NULL, y = full_title))

  # start: not yet tested #
  if(reverse_strand == TRUE) {
    int_pl <- int_pl + scale_x_reverse()
  }
  # stop: not yet tested #

  if(print_plot == T) print(int_pl)

  return(int_pl)

}




#' @title loops or interactions plotter
#' @description This function plots data in bedpe format
#' @param bedpe bedpe file for DNA loops or interaction data
#' @param name_bedpe name of sample
#' @param color_bedpe color of track
#' @param mychr chromosome
#' @param start_loc start location on chromomse
#' @param end_loc end location on chromosome
#' @param breaks_wanted number of breaks on x-axis
#' @param y_title title for y-axis
#' @param show_partial_overlap include loops that do not entirely overlap the coordinates
#' @param axis_line_size line size of x-axis
#' @param marginvec_mm margin around plot, default = c(0,0,0,0)
#' @param reverse_strand reverse the x-axis, default = FALSE
#' @param show_labels boolean: show x-axis lables
#' @param print_plot boolean: print individual plot
#' @import ggplot2 data.table
#' @return ggplot object
#' @keywords loops, bedpe
#' @examples
#'     plot_loops(bedpe)
#'
#'
#' @export
plot_loops = function (bedpe, name_bedpe, color_bedpe,
                       mychr, start_loc, end_loc,
                       breaks_wanted = 4,
                       y_title = "prob",
                       show_partial_overlap = T,
                       axis_line_size = 0.2,
                       marginvec_mm = c(0, 0, 0, 0),
                       reverse_strand = FALSE,
                       show_labels = F,
                       print_plot = F) {
  colnames(bedpe)[1:8] <- c('left_chrom', 'start_left_anchor', 'end_left_anchor',
                            'right_chrom', 'start_right_anchor', 'end_right_anchor',
                            'name_loop', 'counts')

  # 1. select all loops that are overlapping with the selected region
  bedpe_subset <- bedpe[left_chrom == mychr & start_left_anchor <= end_loc & end_right_anchor >=  start_loc]

  # 2. show loops that only partially overlap OR that completely overlap but for which the anchors are outside the selected region
  if (show_partial_overlap == TRUE) {
    bedpe_subset[, start_left_anchor := ifelse(start_left_anchor < start_loc, start_loc, start_left_anchor)]
    bedpe_subset[, end_left_anchor := ifelse(end_left_anchor < start_loc, start_loc, end_left_anchor)]

    bedpe_subset[, start_right_anchor := ifelse(start_right_anchor > end_loc, end_loc, start_right_anchor)]
    bedpe_subset[, end_right_anchor := ifelse(end_right_anchor > end_loc, end_loc, end_right_anchor)]

  } else {
    bedpe_subset = bedpe_subset[start_left_anchor >= start_loc & end_right_anchor <= end_loc]
  }

  bedpe_subset[, `:=`(segm_y_values, 1:nrow(bedpe_subset))]
  calculated_breaks <- seq(start_loc, end_loc, length.out = breaks_wanted)
  if (show_labels == T) {
    mylabels <- paste0(round(calculated_breaks/1000, digits = 0),
                       "kb")
  }
  else {
    mylabels <- rep("", length(calculated_breaks))
  }
  max_y = max(bedpe_subset$segm_y_values) + 1

  # create y-title
  full_title = paste0(name_bedpe, "\n", y_title)


  int_pl <- ggplot()
  int_pl <- int_pl + geom_segment(data = bedpe_subset, aes(x = end_left_anchor,
                                                           y = segm_y_values, xend = start_right_anchor, yend = segm_y_values),
                                  color = color_bedpe)
  int_pl <- int_pl + geom_segment(data = bedpe_subset, aes(x = start_left_anchor,
                                                           y = segm_y_values, xend = end_left_anchor, yend = segm_y_values),
                                  size = 1.2, color = color_bedpe)
  int_pl <- int_pl + geom_segment(data = bedpe_subset, aes(x = start_right_anchor,
                                                           y = segm_y_values, xend = end_right_anchor, yend = segm_y_values),
                                  size = 1.2, color = color_bedpe)
  int_pl <- int_pl + scale_x_continuous(expand = c(0, 0), limits = c(start_loc,
                                                                     end_loc), breaks = calculated_breaks, labels = mylabels)
  int_pl <- int_pl + ylim(c(0, max_y))
  int_pl <- int_pl + theme_bw() + theme(panel.background = element_blank(),
                                        panel.border = element_blank(), axis.line = element_line(color = "black",
                                                                                                 size = axis_line_size), panel.grid = element_blank(),
                                        axis.text.x = element_text(angle = 45, hjust = 45, vjust = 45),
                                        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                        plot.margin = unit(marginvec_mm, "mm"))
  int_pl <- int_pl + labs(x = NULL, y = full_title)
  if (reverse_strand == TRUE) {
    int_pl <- int_pl + scale_x_reverse()
  }
  if (print_plot == T)
    print(int_pl)
  return(int_pl)
}





#' @title bed region plotter
#' @description This function plots data from bed file
#' @param bed bedpe file for DNA loops or interaction data
#' @param name_bed name of sample
#' @param color_bed color of track
#' @param size_bed size of bed (height)
#' @param mychr chromosome
#' @param start_loc start location on chromomse
#' @param end_loc end location on chromosome
#' @param breaks_wanted number of breaks on x-axis
#' @param y_title title for y-axis
#' @param show_partial_overlap include loops that do not entirely overlap the coordinates
#' @param same_y_level plot all bed regions on same height (y-axis)
#' @param axis_line_size line size of x-axis
#' @param marginvec_mm margin around plot, default = c(0,0,0,0)
#' @param show_labels boolean: show x-axis lables
#' @param print_plot boolean: print individual plot
#' @import ggplot2 data.table
#' @return ggplot object
#' @keywords loops, bedpe
#' @examples
#'     plot_bed(bed)
#'
#'
#' @export
plot_bed <- function(bed, name_bed, color_bed, size_bed = 1,
                     mychr, start_loc, end_loc,
                     breaks_wanted = 4, y_title = 'SE',
                     show_partial_overlap = T,
                     same_y_level = T,
                     axis_line_size = 1, marginvec_mm = c(0,0,0,0),
                     reverse_strand = FALSE,
                     show_labels = F, print_plot = F) {


  # rename columns of UCSC BED format
  colnames(bed)[1:6] <- c('chr', 'start', 'stop', 'name', 'score', 'strand')
  setorder(bed, chr, start, stop)

  # if feature needs to be completely within coordinates
  bed_subset <- bed[chr == mychr & start >= start_loc & stop <= end_loc]

  # any feature within coordines
  # if feature needs to be completely within coordinates
  if(show_partial_overlap == TRUE) {

    # cap features for which only the start is within the selected genomic region
    bed_subset_onlyStart <- bed[chr == mychr & start >= start_loc & start <= end_loc]
    bed_subset_onlyStart <- bed_subset_onlyStart[!name %in% bed_subset$name]
    bed_subset_onlyStart[, stop := end_loc]

    # cap features for which only the stop is within the selected genomic region
    bed_subset_onlyStop <- bed[chr == mychr & stop >= start_loc & stop <= end_loc]
    bed_subset_onlyStop <- bed_subset_onlyStop[!name %in% bed_subset$name]
    bed_subset_onlyStop[, start := start_loc]

    bed_subset <- do.call('rbind', list(bed_subset, bed_subset_onlyStart, bed_subset_onlyStop))
    setorder(bed_subset, chr, start, stop)
  }

  # give y-values
  if(same_y_level == TRUE) {
    bed_subset[, segm_y_values := 1]
  } else {
    bed_subset[, segm_y_values := 1:nrow(bed_subset)]
  }


  # create breaks
  calculated_breaks <- seq(start_loc, end_loc, length.out = breaks_wanted)

  if(show_labels == T) {
    mylabels <- paste0(round(calculated_breaks/1000, digits = 0), 'kb')
  } else {
    mylabels <- rep('', length(calculated_breaks))
  }

  # maximum y
  max_y = max(bed_subset$segm_y_values)+4

  # create y title
  full_title = paste0(name_bed, '\n', y_title)


  bed_pl <- ggplot()
  bed_pl <- bed_pl + geom_segment(data = bed_subset,
                                  aes(x = start, y = segm_y_values, xend = stop, yend = segm_y_values),
                                  color = color_bed, size  = size_bed)

  bed_pl <- bed_pl + scale_x_continuous(expand = c(0, 0), limits = c(start_loc, end_loc),
                                        breaks = calculated_breaks, labels = mylabels)
  bed_pl <- bed_pl + ylim(c(0,max_y))
  bed_pl <- bed_pl + theme_bw() + theme(panel.background = element_blank(),
                                        panel.border = element_blank(),
                                        axis.line = element_line(color = 'black', size = axis_line_size),
                                        panel.grid = element_blank(),
                                        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                                        axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        plot.margin=unit(marginvec_mm,"mm"))
  bed_pl <- bed_pl + labs(x = NULL, y = full_title)

  # start: not yet tested #
  if(reverse_strand == TRUE) {
    bed_pl <- bed_pl + scale_x_reverse()
  }
  # stop: not yet tested #

  if(print_plot == T) print(bed_pl)

  return(bed_pl)

}




#' @title genomic coordinate plotter
#' @description This function plots genomic coordinates
#' @param mychr chromosome
#' @param start_loc start location on chromomse
#' @param end_loc end location on chromosome
#' @param breaks_wanted number of breaks on x-axis
#' @param label_size label size for x-axis text, default is 6
#' @param axis_line_size line size of x-axis
#' @param marginvec_mm margin around plot, default = c(0,0,0,0)
#' @param show_labels boolean: show x-axis lables
#' @param print_plot boolean: print individual plot
#' @import ggplot2 data.table
#' @return ggplot object
#' @keywords loops, bedpe
#' @examples
#'     plot_coord(mychr, ...)
#'
#'
#' @export
plot_coord <- function(mychr, start_loc, end_loc,
                       breaks_wanted = 4, label_size = 6,
                       axis_line_size = 0.5, marginvec_mm = c(0,0,0,0),
                       reverse_strand = FALSE,
                       show_labels = T, print_plot = F) {

  # create empty data.frame
  df <- data.frame()

  # create breaks
  calculated_breaks <- seq(start_loc, end_loc, length.out = breaks_wanted)

  if(show_labels == T) {
    mylabels <- paste0(round(calculated_breaks/1000, digits = 0), 'kb')
    #if(reverse_strand == TRUE) {
    #  mylabels <- rev(mylabels)
    #}
    print(mylabels)
  } else {
    mylabels <- rep('', length(calculated_breaks))
    #if(reverse_strand == TRUE) {
    #  mylabels <- rev(mylabels)
    #}
    print(mylabels)
  }

  coord_pl <- ggplot()
  # transcript length
  coord_pl <- coord_pl + geom_segment(data = df)

  if(reverse_strand == TRUE) {
    coord_pl <- coord_pl + scale_x_continuous( expand = c(0, 0), limits = c(end_loc, start_loc), breaks = calculated_breaks, labels = mylabels, trans = 'reverse')
  } else {
    coord_pl <- coord_pl + scale_x_continuous( expand = c(0, 0), limits = c(start_loc, end_loc), breaks = calculated_breaks, labels = mylabels)
  }

  coord_pl <- coord_pl + theme_bw() + theme(panel.background = element_blank(),
                                            panel.border = element_blank(),
                                            axis.line = element_line(color = 'black', axis_line_size),
                                            panel.grid = element_blank(),
                                            axis.text.x = element_text(angle = 45, size = label_size, hjust = 1, vjust = 1),
                                            axis.text.y = element_blank(),
                                            axis.ticks.y = element_blank(),
                                            plot.margin = unit(marginvec_mm, 'mm'))
  coord_pl <- coord_pl + labs(x = NULL, y = '')



  if(print_plot == T) print(coord_pl)

  return(coord_pl)




}



#' @title genome model plotter
#' @description This function plots the genomic model
#' @param transcript transcript in bed format
#' @param exon exon coordinates in bed format
#' @param five_UTR five_UTR coordinates in bed format
#' @param three_UTR three_UTR coordinates in bed format
#' @param mychr chromosome
#' @param start_loc start location on chromomse
#' @param end_loc end location on chromosome
#' @param breaks_wanted number of breaks on x-axis
#' @param show_partial_overlap include loops that do not entirely overlap the coordinates
#' @param y_coord_gene_dist distance between genes
#' @param y_top_space additional space on top
#' @param y_genename_space space between gene name and gene model
#' @param genename_size size of gene name
#' @param marginvec_mm margin around plot, default = c(0,0,0,0)
#' @param show_labels boolean: show x-axis lables
#' @param print_plot boolean: print individual plot
#' @import ggplot2 data.table
#' @return ggplot object
#' @keywords loops, bedpe
#' @examples
#'     plot_genome(mychr, ...)
#'
#'
#' @export
plot_genome <- function(transcript, exon, five_UTR, three_UTR,
                        mychr = 'chr14', start_loc = 61100000 , end_loc = 61195000,
                        breaks_wanted = 4,
                        show_partial_overlap = T, marginvec_mm = c(0,0,0,0),
                        y_coord_gene_dist = 0.5, y_top_space = 2,
                        y_genename_space = 1, genename_size = 4,
                        exon_size = 0.5, UTR_size = 0.5,
                        reverse_strand = FALSE,
                        show_labels = T, print_plot = T) {


  ## TO DO ##
  # if there are only partially overlapping genes, then it does not work: FIX

  # UCSC BED format is strict! and first 6 are minimally required



  # 1. select transcripts that are within range limits
  colnames(transcript)[1:6] <- c('chr', 'start', 'end', 'name', 'score', 'strand')


  subset_transcript <- transcript[chr == mychr & start >= start_loc & end <= end_loc]

  if(nrow(subset_transcript) > 0) {
    selected_transcripts <- unique(subset_transcript[['name']])
    subset_transcript[, c('xstart', 'xend') := list(start, end)]
    subset_transcript[, gene_name_center := round(mean(c(xstart, xend)), digits = 0), by = 1:nrow(subset_transcript)]
  } else {
    subset_transcript <- NULL
  }


  #cat('1 \n')



  if(show_partial_overlap == TRUE) {

    subset_transcript_OnlyStart <- transcript[chr == mychr & start >= start_loc & start <= end_loc]
    if(!is.null(subset_transcript)) {
      subset_transcript_OnlyStart <- subset_transcript_OnlyStart[!name %in% subset_transcript$name]
    }
    subset_transcript_OnlyStart[, end := end_loc]
    subset_transcript_OnlyStart[, c('xstart', 'xend') := list(start, end)]

    if(nrow(subset_transcript_OnlyStart) > 0) {
      subset_transcript_OnlyStart[, gene_name_center := round(mean(c(xstart, xend)), digits = 0), by = 1:nrow(subset_transcript_OnlyStart)]
    } else {
      subset_transcript_OnlyStart[, gene_name_center := NA]
    }


    subset_transcript_OnlyEnd <- transcript[chr == mychr & end >= start_loc & end <= end_loc]
    if(!is.null(subset_transcript)) {
      subset_transcript_OnlyEnd <- subset_transcript_OnlyEnd[!name %in% subset_transcript$name]
    }
    subset_transcript_OnlyEnd[, start := start_loc]
    subset_transcript_OnlyEnd[, c('xstart', 'xend') := list(start, end)]

    if(nrow(subset_transcript_OnlyEnd) > 0) {
      subset_transcript_OnlyEnd[, gene_name_center := round(mean(c(xstart, xend)), digits = 0), by = 1:nrow(subset_transcript_OnlyEnd)]
    } else {
      subset_transcript_OnlyEnd[, gene_name_center := NA]
    }

    subset_transcript_all_overlap <- do.call('rbind', list(subset_transcript, subset_transcript_OnlyStart, subset_transcript_OnlyEnd))
    subset_transcript <- subset_transcript_all_overlap
    selected_transcripts <- unique(subset_transcript[['name']])

  }

  #cat('2 \n')

  setorder(subset_transcript, chr, start, end)
  subset_transcript_plus <- subset_transcript[strand == '+']
  subset_transcript_minus <- subset_transcript[strand == '-']


  # calculate y-coord for each gene
  # y-distance between genes and size of exon / UTR
  # number of start positions can be adjusted
  plus_rows <- nrow(subset_transcript_plus)
  minus_rows <- nrow(subset_transcript_minus)


  y_dist = y_coord_gene_dist

  if(plus_rows > 0) {
    plus_y = c(1*y_dist, 4*y_dist)
    times_plus <- ceiling(plus_rows/(length(plus_y)))
    plus_starts <- rep(plus_y, times_plus)
    plus_starts <- plus_starts[1:plus_rows]
    subset_transcript_plus[, ypos := plus_starts]
  } else {
    subset_transcript_plus <- NULL
  }

  if(minus_rows > 0) {
    minus_y = c(7*y_dist, 10*y_dist)
    times_minus <- ceiling(minus_rows/(length(minus_y)))
    minus_starts <- rep(minus_y, times_minus)
    minus_starts <- minus_starts[1:minus_rows]
    subset_transcript_minus[, ypos := minus_starts]
  } else {
    subset_transcript_minus <- NULL
  }


  subset_transcript <- rbind(subset_transcript_plus, subset_transcript_minus)


  #cat('3 \n')






  # 2. get exons for selected transcripts
  colnames(exon)[1:6] <- c('chr', 'start', 'end', 'name', 'score', 'strand')

  subset_exon <- exon[name %in% selected_transcripts]
  subset_exon[, c('xstart', 'xend') := list(start, end)]

  # remove exons that are outside the window
  subset_exon <- subset_exon[!xend < start_loc]
  subset_exon <- subset_exon[!start > end_loc]


  #cat('4 \n')

  if(show_partial_overlap == TRUE) {

    subset_exon[, xstart := ifelse(xstart >= start_loc & xstart <= end_loc, xstart,
                                   ifelse(xend <= start_loc, xstart,
                                          ifelse(xstart <= start_loc & xend <= end_loc, start_loc, xstart)))]

    subset_exon[, xend := ifelse(xend >= start_loc & xend <= end_loc, xend,
                                 ifelse(xstart >= end_loc, xend,
                                        ifelse(xstart <= end_loc & xend >= end_loc, end_loc, xend)))]

  }

  subset_exon <- merge(subset_exon, subset_transcript[, .(name, ypos)], by = 'name')



  # 3. get 5'UTR and 3'UTR regions
  colnames(five_UTR)[1:6] <- c('chr', 'start', 'end', 'name', 'score', 'strand')
  colnames(three_UTR)[1:6] <- c('chr', 'start', 'end', 'name', 'score', 'strand')

  # 5 UTR
  subset_five_UTR <- five_UTR[name %in% selected_transcripts]
  subset_five_UTR[, c('xstart', 'xend') := list(start, end)]

  # remove exons that are outside the window
  subset_five_UTR <- subset_five_UTR[!xend < start_loc]
  subset_five_UTR <- subset_five_UTR[!start > end_loc]




  if(show_partial_overlap == TRUE) {

    subset_five_UTR[, xstart := ifelse(xstart >= start_loc & xstart <= end_loc, xstart,
                                       ifelse(xend <= start_loc, xstart,
                                              ifelse(xstart <= start_loc & xend <= end_loc, start_loc, xstart)))]

    subset_five_UTR[, xend := ifelse(xend >= start_loc & xend <= end_loc, xend,
                                     ifelse(xstart >= end_loc, xend,
                                            ifelse(xstart <= end_loc & xend >= end_loc, end_loc, xend)))]
  }


  subset_five_UTR <- merge(subset_five_UTR, subset_transcript[, .(name, ypos)], by = 'name')


  # three UTR
  subset_three_UTR <- three_UTR[name %in% selected_transcripts]
  subset_three_UTR[, c('xstart', 'xend') := list(start, end)]

  # remove exons that are outside the window
  subset_three_UTR <- subset_three_UTR[!xend < start_loc]
  subset_three_UTR <- subset_three_UTR[!start > end_loc]


  if(show_partial_overlap == TRUE) {

    subset_three_UTR[, xstart := ifelse(xstart >= start_loc & xstart <= end_loc, xstart,
                                        ifelse(xend <= start_loc, xstart,
                                               ifelse(xstart <= start_loc & xend <= end_loc, start_loc, xstart)))]

    subset_three_UTR[, xend := ifelse(xend >= start_loc & xend <= end_loc, xend,
                                      ifelse(xstart >= end_loc, xend,
                                             ifelse(xstart <= end_loc & xend >= end_loc, end_loc, xend)))]
  }
  subset_three_UTR <- merge(subset_three_UTR, subset_transcript[, .(name, ypos)], by = 'name')


  # create breaks
  calculated_breaks <- seq(start_loc, end_loc, length.out = breaks_wanted)

  if(show_labels == T) {
    mylabels <- paste0(round(calculated_breaks/1000, digits = 0), 'kb')
  } else {
    mylabels <- rep('', length(calculated_breaks))
  }


  # maximum y
  max_y = max(subset_transcript$ypos)+y_top_space



  ### visualization ###
  gene_model_pl <- ggplot()

  ## transcript length ##
  gene_model_pl <- gene_model_pl + geom_segment(data = subset_transcript, aes(x = xstart, y = ypos, xend = xend, yend = ypos))
  gene_model_pl <- gene_model_pl + annotate('text', x = subset_transcript$gene_name_center, y = subset_transcript$ypos+y_genename_space, label = subset_transcript$name, size = genename_size)

  ## exon ##
  gene_model_pl <- gene_model_pl + geom_rect(data = subset_exon, aes(xmin = xstart, ymin = ypos-exon_size, xmax = xend, ymax = ypos+exon_size),fill = 'black' ,color = 'black', size  = 0.5)

  ## UTR coord ##
  if(!is.null(subset_five_UTR)) {
    gene_model_pl <- gene_model_pl + geom_rect(data = subset_five_UTR, aes(xmin = xstart, ymin = ypos-UTR_size, xmax = xend, ymax = ypos+UTR_size), size = 0.5, fill = 'grey')
  }

  if(!is.null(subset_three_UTR)) {
    gene_model_pl <- gene_model_pl + geom_rect(data = subset_three_UTR, aes(xmin = xstart, ymin = ypos-UTR_size, xmax = xend, ymax = ypos+UTR_size), size = 0.5, fill = 'grey')
  }


  if(reverse_strand == TRUE) {
    gene_model_pl <- gene_model_pl +scale_x_continuous( expand = c(0, 0), limits = c(end_loc, start_loc), breaks = calculated_breaks, labels = mylabels, trans = 'reverse')
  } else {
    gene_model_pl <- gene_model_pl + scale_x_continuous( expand = c(0, 0), limits = c(start_loc, end_loc), breaks = calculated_breaks, labels = mylabels)
  }

  gene_model_pl <- gene_model_pl + ylim(c(0, max_y))
  gene_model_pl <- gene_model_pl + theme_bw() + theme(panel.background = element_blank(),
                                                      panel.border = element_blank(),
                                                      axis.line = element_blank(),
                                                      panel.grid = element_blank(),
                                                      axis.text.x = element_blank(),
                                                      axis.text.y = element_blank(),
                                                      axis.ticks.y = element_blank(),
                                                      axis.ticks.x = element_blank(),
                                                      plot.margin=unit(marginvec_mm,"mm"))
  gene_model_pl <- gene_model_pl + labs(list(x = NULL, y = 'genes'))


  if(print_plot == T) print(gene_model_pl)

  return(gene_model_pl)


}





#' @title Genomic tracks plotter
#' @description Main function which is a wrapper for individual plotters and combinations of tracks
#' @param format_vec transcript in bed format
#' @param figure_list exon coordinates in bed format
#' @param uniq_name_vec five_UTR coordinates in bed format
#' @param color_vec three_UTR coordinates in bed format
#' @param mychr chromosome
#' @param start_loc start location on chromomse
#' @param end_loc end location on chromosome
#' @param coord_axis_line_size line size for coordinate axis
#' @param marginvec_mm margin around plot, default = c(0,0,0,0)
#' @param model_show_partial_overlap include loops that do not entirely overlap the coordinates
#' @param transcript transcript in bed format
#' @param exon exon in bed format
#' @param five_UTR five_UTR in bed format
#' @param three_UTR three_UTR in bed format
#' @param breaks_wanted number of breaks on x-axis
#' @param label_size size of labels
#' @param rna_bdg_binsize = 100,
#' @param rna_bdg_y_title = 'RPM/bp',
#' @param rna_bdg_forced_y_min_vec = NULL
#' @param rna_bdg_forced_y_max_vec = NULL,
#' @param rna_bdg_axis_line_size = 0.2,
#' @param rna_bdg_chrom_col = 'chr',
#' @param rna_bdg_start_col = 'start',
#' @param rna_bdg_end_col = 'end',
#' @param rna_bdg_counts_col = 'RPKM',
#' @param rna_bdg_strand_col = NA,
#' @param rna_min_strand_sign = '-',
#' @param rna_plus_strand_sign = '+',
#' @param rna_make_minus_strand_negative = F,
#' @param bdg_binsize = 100,
#' @param bdg_y_title = 'RPM/bp',
#' @param bdg_forced_y_min_vec = NULL,
#' @param bdg_forced_y_max_vec = NULL,
#' @param bdg_axis_line_size = 0.2,
#' @param bdg_chrom_col = 'chr',
#' @param bdg_start_col = 'start',
#' @param bdg_end_col = 'end',
#' @param bdg_counts_col = 'RPKM',
#' @param loop_y_title = 'prob',
#' @param loop_show_partial_overlap = T,
#' @param loops_axis_line_size = 0.2,
#' @param bed_y_title = 'SE',
#' @param bed_size_bed = 4,
#' @param bed_show_partial_overlap = T,
#' @param bed_same_y_level = T,
#' @param bed_axis_line_size = 1,
#' @param show_labels boolean: show x-axis lables
#' @param print_plot boolean: print individual plot
#' @import cowplot ggplot2 data.table
#' @return ggplot object
#' @keywords loops, bedpe, bed, bedgraph, dna, rna, main
#' @examples
#'     Genomic_tracks_plot(format_vec, ...)
#'
#'
#' @export
Genomic_tracks_plot <- function(format_vec, figure_list, uniq_name_vec, color_vec,
                                mychr = 'chr14', start_loc = 61100000 , end_loc = 61195000,
                                coord_axis_line_size = 0.5, marginvec_mm = c(0,0,0,0),
                                reverse_strand = FALSE,
                                # gene model
                                model_show_partial_overlap = T,
                                transcript, exon, five_UTR, three_UTR,
                                breaks_wanted = 4,  label_size = 10,
                                exon_size = 0.5, UTR_size = 0.5,
                                y_coord_gene_dist = 0.5, y_top_space = 2,
                                y_genename_space = 1, genename_size = 4,
                                # RNA bedgraph specific params
                                rna_bdg_binsize = 100,
                                rna_bdg_y_title = 'RPM/bp',
                                rna_bdg_forced_y_min_vec = NULL, rna_bdg_forced_y_max_vec = NULL,
                                rna_bdg_axis_line_size = 0.2,
                                rna_bdg_chrom_col = 'chr', rna_bdg_start_col = 'start', rna_bdg_end_col = 'end',
                                rna_bdg_counts_col = 'RPKM', rna_bdg_strand_col = NA,
                                rna_min_strand_sign = '-', rna_plus_strand_sign = '+',
                                rna_make_minus_strand_negative = F,
                                rna_use_barplot = F, rna_alpha_fill = 0.7,
                                # bedgraph specific params
                                bdg_binsize = 100,
                                bdg_y_title = 'RPM/bp',
                                bdg_forced_y_min_vec = NULL, bdg_forced_y_max_vec = NULL,
                                bdg_axis_line_size = 0.2,
                                bdg_chrom_col = 'chr', bdg_start_col = 'start', bdg_end_col = 'end', bdg_counts_col = 'RPKM',
                                # bedpe specific params
                                loop_y_title = 'prob',
                                loop_show_partial_overlap = T,
                                loops_axis_line_size = 0.2,
                                # bed specific params
                                bed_y_title = 'SE', bed_size_bed = 4,
                                bed_show_partial_overlap = T, bed_same_y_level = T,
                                bed_axis_line_size = 1,
                                # cowplot: plotgrid params
                                show_labels = F, print_plot = T, ...) {


  # libraries
  library(ggplot2)
  library(data.table)
  library(cowplot)

  # preparation
  names(format_vec) <- uniq_name_vec
  names(figure_list) <- uniq_name_vec
  names(color_vec) <- uniq_name_vec


  # to create different limits for each bedgraph
  if(!is.null(bdg_forced_y_min_vec)) {
    min_limit_vec_names <- names(format_vec[format_vec == 'bedgraph'])
    names(bdg_forced_y_min_vec) <- min_limit_vec_names
  }


  if(!is.null(bdg_forced_y_max_vec)) {
    max_limit_vec_names <- names(format_vec[format_vec == 'bedgraph'])
    names(bdg_forced_y_max_vec) <- max_limit_vec_names
  }




  ## START TEST ##
  # to create different limits for each rna bedgraph
  if(!is.null(rna_bdg_forced_y_min_vec)) {
    min_limit_vec_names <- names(format_vec[format_vec == 'rna'])
    names(rna_bdg_forced_y_min_vec) <- min_limit_vec_names
    print(rna_bdg_forced_y_min_vec)
  }


  if(!is.null(rna_bdg_forced_y_max_vec)) {
    max_limit_vec_names <- names(format_vec[format_vec == 'rna'])
    names(rna_bdg_forced_y_max_vec) <- max_limit_vec_names
    print(rna_bdg_forced_y_max_vec)
  }

  ## STOP TEST ##






  save_figure_list <- list()

  for(plot in 1:length(uniq_name_vec)) {

    plotname = uniq_name_vec[[plot]]
    plotformat = format_vec[[plotname]]
    plotcolor = color_vec[[plotname]]
    plotdata = figure_list[[plotname]]




    if(plotformat == 'rna') {


      cat('\n for ',plotname,' create RNA bedgraph \n \n')


      if(!is.null(rna_bdg_forced_y_min_vec)) {
        rna_bdg_forced_y_min = rna_bdg_forced_y_min_vec[[plotname]]
        print(rna_bdg_forced_y_min)
      } else {
        rna_bdg_forced_y_min = NULL
      }


      if(!is.null(rna_bdg_forced_y_max_vec)) {
        rna_bdg_forced_y_max = rna_bdg_forced_y_max_vec[[plotname]]
        print(rna_bdg_forced_y_max)
      } else {
        rna_bdg_forced_y_max = NULL
      }




      newplot <- plot_RNA_bedgraph(rna_bdg = plotdata, name_rna_bdg = plotname, color_rna_bdg = plotcolor,
                                   mychr = mychr, start_loc = start_loc , end_loc = end_loc, binsize = rna_bdg_binsize,
                                   breaks_wanted = breaks_wanted, y_title = rna_bdg_y_title,
                                   chrom_col = rna_bdg_chrom_col, start_col = rna_bdg_start_col, end_col = rna_bdg_end_col,
                                   counts_col = rna_bdg_counts_col, strand_col = rna_bdg_strand_col,
                                   min_strand_sign = rna_min_strand_sign, plus_strand_sign = rna_plus_strand_sign,
                                   make_minus_strand_negative = rna_make_minus_strand_negative,
                                   forced_y_min = rna_bdg_forced_y_min, forced_y_max = rna_bdg_forced_y_max,
                                   axis_line_size = rna_bdg_axis_line_size, marginvec_mm = marginvec_mm,
                                   show_labels = show_labels, print_plot = print_plot,
                                   reverse_strand = reverse_strand)


    }




    if(plotformat == 'bedgraph') {

      cat('\n for ',plotname,' create bedgraph \n \n')


      if(!is.null(bdg_forced_y_min_vec)) {
        bdg_forced_y_min = bdg_forced_y_min_vec[[plotname]]
      } else {
        bdg_forced_y_min = NULL
      }


      if(!is.null(bdg_forced_y_max_vec)) {
        bdg_forced_y_max = bdg_forced_y_max_vec[[plotname]]
      } else {
        bdg_forced_y_max = NULL
      }



      newplot <- plot_bedgraph(bdg = plotdata, name_bdg = plotname, color_bdg = plotcolor,
                               mychr = mychr, start_loc = start_loc, end_loc = end_loc, binsize = bdg_binsize,
                               breaks_wanted = breaks_wanted, y_title = bdg_y_title,
                               chrom_col = bdg_chrom_col, start_col = bdg_start_col, end_col = bdg_end_col, counts_col = bdg_counts_col,
                               forced_y_min = bdg_forced_y_min, forced_y_max = bdg_forced_y_max,
                               axis_line_size = bdg_axis_line_size, marginvec_mm = marginvec_mm,
                               reverse_strand = reverse_strand,
                               show_labels = show_labels, print_plot = print_plot)

    } else if(plotformat == 'bedpe') {

      cat('\n for ',plotname,' create bedpe/loops \n \n')

      newplot <- plot_loops(bedpe = plotdata, name_bedpe = plotname, color_bedpe = plotcolor,
                            mychr = mychr, start_loc = start_loc, end_loc = end_loc,
                            breaks_wanted = breaks_wanted, y_title = loop_y_title,
                            show_partial_overlap = loop_show_partial_overlap,
                            axis_line_size = loops_axis_line_size, marginvec_mm = marginvec_mm,
                            reverse_strand = reverse_strand,
                            show_labels = show_labels, print_plot = print_plot)

    } else if(plotformat == 'bed') {

      cat('\n for ',plotname,' create bed \n \n')

      newplot <- plot_bed(bed = plotdata, name_bed = plotname,
                          color_bed = plotcolor, size_bed = bed_size_bed,
                          mychr = mychr, start_loc = start_loc , end_loc = end_loc,
                          breaks_wanted = breaks_wanted, y_title = bed_y_title,
                          show_partial_overlap = bed_show_partial_overlap,
                          same_y_level = bed_same_y_level,
                          axis_line_size = bed_axis_line_size, marginvec_mm = marginvec_mm,
                          reverse_strand = reverse_strand,
                          show_labels = show_labels, print_plot = print_plot)


    } else if(plotformat == 'model') {

      cat('\n for ',plotname,' create gene model \n \n')

      newplot <- plot_genome(transcript = transcript, exon = exon, five_UTR = five_UTR, three_UTR = three_UTR,
                             mychr = mychr, start_loc = start_loc , end_loc = end_loc,
                             breaks_wanted = breaks_wanted,
                             show_partial_overlap = model_show_partial_overlap, marginvec_mm = marginvec_mm,
                             exon_size = exon_size, UTR_size = UTR_size,
                             y_coord_gene_dist = y_coord_gene_dist, y_top_space = y_top_space,
                             y_genename_space = y_genename_space, genename_size = genename_size,
                             reverse_strand = reverse_strand,
                             show_labels = show_labels, print_plot = print_plot)

    } else if(plotformat == 'coord') {

      cat('\n for ',plotname,' create coordinate system \n \n')

      newplot <- plot_coord(mychr = mychr, start_loc = start_loc, end_loc = end_loc,
                            breaks_wanted = breaks_wanted, label_size = label_size,
                            axis_line_size = coord_axis_line_size, marginvec_mm = marginvec_mm,
                            reverse_strand = reverse_strand,
                            show_labels = T, print_plot = print_plot)

    } else if(!plotformat %in% c('rna', 'bed', 'bedgraph', 'bedpe', 'model', 'coord')) {
      cat('\n \n format is not know, use bedgraph, bedpe, model or coord \n \n')
    }


    save_figure_list[[plotname]] = newplot


  }

  # return together

  if(print_plot == T) {
    final_plot <- cowplot::plot_grid(plotlist = save_figure_list, ...)
    print(final_plot)
  }
  return(save_figure_list)

}



