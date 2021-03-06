% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeTrackVIz2_main.R
\name{Genomic_tracks_plot}
\alias{Genomic_tracks_plot}
\title{Genomic tracks plotter}
\usage{
Genomic_tracks_plot(format_vec, figure_list, uniq_name_vec, color_vec,
  mychr = "chr14", start_loc = 61100000, end_loc = 61195000,
  coord_axis_line_size = 0.5, marginvec_mm = c(0, 0, 0, 0),
  reverse_strand = FALSE, model_show_partial_overlap = T, transcript,
  exon, five_UTR, three_UTR, breaks_wanted = 4, label_size = 10,
  exon_size = 0.5, UTR_size = 0.5, y_coord_gene_dist = 0.5,
  y_top_space = 2, y_genename_space = 1, genename_size = 4,
  rna_bdg_binsize = 100, rna_bdg_y_title = "RPM/bp",
  rna_bdg_forced_y_min_vec = NULL, rna_bdg_forced_y_max_vec = NULL,
  rna_bdg_axis_line_size = 0.2, rna_bdg_chrom_col = "chr",
  rna_bdg_start_col = "start", rna_bdg_end_col = "end",
  rna_bdg_counts_col = "RPKM", rna_bdg_strand_col = NA,
  rna_min_strand_sign = "-", rna_plus_strand_sign = "+",
  rna_make_minus_strand_negative = F, rna_use_barplot = F,
  rna_alpha_fill = 0.7, bdg_binsize = 100, bdg_y_title = "RPM/bp",
  bdg_forced_y_min_vec = NULL, bdg_forced_y_max_vec = NULL,
  bdg_axis_line_size = 0.2, bdg_chrom_col = "chr",
  bdg_start_col = "start", bdg_end_col = "end",
  bdg_counts_col = "RPKM", loop_y_title = "prob",
  loop_show_partial_overlap = T, loops_axis_line_size = 0.2,
  loops_color_column = NULL, loops_color_bedpe = "black",
  bed_y_title = "SE", bed_size_bed = 4, bed_show_partial_overlap = T,
  bed_same_y_level = T, bed_axis_line_size = 1, show_labels = F,
  print_plot = T, ...)
}
\arguments{
\item{format_vec}{transcript in bed format}

\item{figure_list}{exon coordinates in bed format}

\item{uniq_name_vec}{five_UTR coordinates in bed format}

\item{color_vec}{three_UTR coordinates in bed format}

\item{mychr}{chromosome}

\item{start_loc}{start location on chromomse}

\item{end_loc}{end location on chromosome}

\item{coord_axis_line_size}{line size for coordinate axis}

\item{marginvec_mm}{margin around plot, default = c(0,0,0,0)}

\item{model_show_partial_overlap}{include loops that do not entirely overlap the coordinates}

\item{transcript}{transcript in bed format}

\item{exon}{exon in bed format}

\item{five_UTR}{five_UTR in bed format}

\item{three_UTR}{three_UTR in bed format}

\item{breaks_wanted}{number of breaks on x-axis}

\item{label_size}{size of labels}

\item{rna_bdg_binsize}{= 100,}

\item{rna_bdg_y_title}{= 'RPM/bp',}

\item{rna_bdg_forced_y_min_vec}{= NULL}

\item{rna_bdg_forced_y_max_vec}{= NULL,}

\item{rna_bdg_axis_line_size}{= 0.2,}

\item{rna_bdg_chrom_col}{= 'chr',}

\item{rna_bdg_start_col}{= 'start',}

\item{rna_bdg_end_col}{= 'end',}

\item{rna_bdg_counts_col}{= 'RPKM',}

\item{rna_bdg_strand_col}{= NA,}

\item{rna_min_strand_sign}{= '-',}

\item{rna_plus_strand_sign}{= '+',}

\item{rna_make_minus_strand_negative}{= F,}

\item{bdg_binsize}{= 100,}

\item{bdg_y_title}{= 'RPM/bp',}

\item{bdg_forced_y_min_vec}{= NULL,}

\item{bdg_forced_y_max_vec}{= NULL,}

\item{bdg_axis_line_size}{= 0.2,}

\item{bdg_chrom_col}{= 'chr',}

\item{bdg_start_col}{= 'start',}

\item{bdg_end_col}{= 'end',}

\item{bdg_counts_col}{= 'RPKM',}

\item{loop_y_title}{= 'prob',}

\item{loop_show_partial_overlap}{= T,}

\item{loops_axis_line_size}{= 0.2,}

\item{loops_color_column}{= NULL,}

\item{loops_color_bedpe}{= 'black',}

\item{bed_y_title}{= 'SE',}

\item{bed_size_bed}{= 4,}

\item{bed_show_partial_overlap}{= T,}

\item{bed_same_y_level}{= T,}

\item{bed_axis_line_size}{= 1,}

\item{show_labels}{boolean: show x-axis lables}

\item{print_plot}{boolean: print individual plot}
}
\value{
ggplot object
}
\description{
Main function which is a wrapper for individual plotters and combinations of tracks
}
\examples{
    Genomic_tracks_plot(format_vec, ...)


}
\keyword{bed,}
\keyword{bedgraph,}
\keyword{bedpe,}
\keyword{dna,}
\keyword{loops,}
\keyword{main}
\keyword{rna,}
