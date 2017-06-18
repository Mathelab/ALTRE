#' THIS FUNCTION FOR USE WHEN RUNNING ALTRE ON MACOS/LINUX ONLY!
#' Generates count data for regulatory regions.
#'
#' Counts the number of reads in each regulatory region for each sample type --
#' read count is derived from user-input BAM files, and regions of interest
#' are supplied in a GRanges object, ideally output of combineAnnotatePeaks.
#' The function getCounts generates count data using
#' the featureCounts function from the R package Rsubreads, which is the fastest
#' way available on R to count. The getCountsFast function CANNOT be used when
#' running ALTRE on a Windows computer. Windows computers must use the function
#' getCounts (also available in the ALTRE package), which is significantly
#' slower, but ultimately will give the exact same results. For high-thoughput
#' experiments (many samples need to be analyzed), it is highly suggested that
#' a non-Windows computer is used (MacOS/Linux).
#'

#'
#' @param annotpeaks list output from combineAnnotatePeaks() function
#' @param sampleinfo dataframe as returned from loadCSVFile() function
#' @param reference name of sample type to be
#' considered 'reference' in DESeq2 analysis
#' @param singleEnd whether input data is single-end (default is TRUE)
#' @param chrom optional, only chromosome chrom will be evaluated
#'
#' @return List containing three items:
#' (1) DESeqDataSet: contains count information for all replicates of all samples
#' (2) Matrix: contains number of TSS-distal and TSS-proximal
#'  before and after filtering (if applicable)
#' (3) Data frame for creating a density plot (use function plotGetCounts()
#'
#'
#' @examples
#' \dontrun{
#' csvfile <- loadCSVFile("DNaseEncodeExample.csv")
#' samplePeaks <- loadBedFiles(csvfile)
#' consensusPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps = 2)
#' TSSannot <- getTSS()
#' consensusPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consensusPeaks,
#'    TSS = TSSannot,
#'    merge = TRUE,
#'    regionspecific = TRUE,
#'    distancefromTSSdist = 1500,
#'    distancefromTSSprox = 1000)
#' consensusPeaksCounts <- getCountsFast(annotpeaks = consensusPeaksAnnotated,
#'    sampleinfo = csvfile,
#'    reference = 'SAEC',
#'    chrom = 'chr21')
#' }
#' @export

getCountsFast <- function(annotpeaks,
                      sampleinfo,
                      reference,
                      singleEnd = TRUE,
                      chrom = NULL) {

  if (!requireNamespace("Rsubread", quietly = TRUE)) {
    stop("The Rsubread package is required for this function to work.
Rsubread is only available for computers running MacOS or Linux.
If you are running MacOS or Linux please install Rsubread:
https://bioconductor.org/packages/release/bioc/html/Rsubread.html.
If you are using a Windows computer please use the function getCounts
(also available in the ALTRE R package) in place of the current function.
The results of the analysis will be exactly the same, but the processing
time may be much slower. If you plan to analyze a very large number of samples
it is highly suggested that a non-Windows computer is used (MacOS/Linux).",
         call. = FALSE)
  }

#  datapaths <- paste(sampleinfo$datapath, sampleinfo$bamfiles, sep = "/")
  datapaths <- file.path(sampleinfo$datapath, sampleinfo$bamfiles)
  #subset by chromosome if necessary
  if (is.null(chrom) == FALSE) {
    regions <- annotpeaks[[1]][seqnames(annotpeaks[[1]]) == chrom,]
  } else {
    regions <- annotpeaks[[1]]
  }

  # create correct format (SAF) in order to use featureCounts
  regionsdataframe = as.data.frame(regions)
  SAFformat = as.data.frame(paste(regionsdataframe[,1], regionsdataframe[,2], regionsdataframe[,3], sep = ":"))
  SAFformat = cbind(SAFformat, regionsdataframe[,c(1:3)])
  SAFformat$Strand = "-"
  colnames(SAFformat) = c("GeneID", "Chr", "Start", "End", "Strand")
  sampleinfo$bamfiles
  results = Rsubread::featureCounts(files = datapaths, annot.ext = SAFformat, isPairedEnd = !singleEnd)
  counts = results$counts
  colnames(counts) = sampleinfo$bamfiles

  # get region/peak size
  originaldata <- grangestodataframe(regions)
  regionsize <- originaldata$stop - originaldata$start

  # Calculate RPKM for plotting densities multiply by 10^6 and divide by
  # regions size to get rpkm
  myrpkm <- as.data.frame((counts[, 1] * 10 ^ 6)/regionsize)
  for (i in 2:ncol(counts)) {
    myrpkm[, i] <- counts[, i] * 10 ^ 6/regionsize
  }
  # take the log2 so that it is a normalized distribution
  myrpkmlog2 <- log2(as.matrix(myrpkm) + 1)
  colnames(myrpkmlog2) <- unlist(lapply(paste(sampleinfo$sample,
                                              sampleinfo$replicate,
                                              sep = "_"),
                                        as.character)
                                 )

  #########################################
  # Create stats matrix originaldata is created ~ 10 lines lines above
  colnames(originaldata) <- unlist(lapply(colnames(originaldata), gsub,
                                          pattern = "meta.", replacement = ""))
  tssdistnum <- length(which(originaldata$region == "TSS-distal"))
  tssproxnum <- length(which(originaldata$region == "TSS-proximal"))

  statdf <- data.frame(Num_TSSdistals = tssdistnum,
                       Num_TSSproximals = tssproxnum)

  #########################################
  # Create densityplot
  region <- originaldata$region
  forplotdf <- cbind(myrpkmlog2, as.data.frame(region))

  ##########################################
  # Create a summerizedExperiment Object
  sampleinfo$status <- stats::relevel(as.factor(sampleinfo$sample), reference)
  colnames(counts) <- NULL
  SEcounts <- DESeq2::DESeqDataSetFromMatrix(counts, as.data.frame(sampleinfo[c(1:4,6)]), design = ~status)
  SummarizedExperiment::rowRanges(SEcounts) = regions
  colnames(SummarizedExperiment::rowData(SEcounts)) <-
    gsub("meta.","", colnames(SummarizedExperiment::rowData(SEcounts)))


  return(list(regioncounts = SEcounts, regioncountstats = statdf,
              regioncountsforplot = forplotdf, reference))
}
