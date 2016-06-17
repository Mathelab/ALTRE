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
#' 2.\code{\link{getConsensusPeaks}}
#'
#' Takes in a sample peaks list (output from loadPeaks in step 1), and outputs consensus peaks.  Consensus peaks are those present in at least N replicates, where is defined by the user.
#' The function outputs a list with the following slots: 1) GRangesList with consensus peaks for each sample type; 2) A statistics table summarizing how many peaks are in replicate samples and are called as consensus peaks.
#'
#' 3.\code{\link{combineAnnotatePeaks}}
#'
#' The GRanges for all sample types are combined and annotated with
#' type specificity (which cell types the hotspot is present in) and
#' whether each region represented in the GRanges is a promoter (default: <1500bp from
#' a transcription start site) or an enhancer (>1500bp from a transcription
#' start site). Function can also merge regulatory regions that are within a specified
#' distance from each other. 
