#' Workflow: Post-alignment to altered enhancers and promoters
#'
#' This workflow takes alignment data and hotspot calls from
#' assays of open chromatin (ATAC-seq and Dnase-seq) and identifies
#' regions (enhancers and promoters) that differ based on cell and
#' tissue type.
#'
#' \pkg{ALTRE} requires sample information CSV file, peak files (bed format), and alignment bam files (BAM) as input
#'
#' The following is the order in which the functions should be used:
#' (Click on function to get more detailed information.)
#'
#' 1. \code{\link{loadPeaks}}
#'
#' Takes in a sample information file (CSV), loads peak files, and outputs a GRangesList object that holds all peaks for each sample type.  
#'
"_PACKAGE"
