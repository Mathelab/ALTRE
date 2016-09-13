#' Generates count data for regulatory regions
#'
#' Counts the number of reads in each regulatory region for each sample type --
#' read count is derived from user-input BAM filex, and regions of interest
#' are supplied in a Granges object, ideally output of combineannotatepeaks.R.
#'
#' @param annotpeaks list output from combineannotatepeaks function
#' @param sampleinfo dataframe as returned from loadCSVFile() function
#' @param reference name of sample type to be
#' considered 'reference' in DESeq2 analysis
#' @param chrom optional, only chromosome chrom will be evaluated
#'
#' @return List containing three items:
#' (1) DESeqDataSet: contains count information for all replicates of all samples
#' (2) Matrix: contains number of TSS-distal and TSS-proximal
#'  before and after filtering (if applicable)
#' (3) Data frame for creating a density plot (use function plotgetcounts()
#'
#'
#' @examples
#' \dontrun{
#' csvfile <- loadCSVFile("DNAseEncodeExample.csv")
#' samplePeaks <- loadBedFiles(csvfile)
#' consensusPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps = 2)
#' TSSannot <- getTSS()
#' consensusPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consensusPeaks,
#'    TSS = TSSannot,
#'    merge = TRUE,
#'    regionspecific = TRUE,
#'    distancefromTSSdist = 1500,
#'    distancefromTSSprox = 1000)
#' consensusPeaksCounts <- getCounts(annotpeaks = consensusPeaksAnnotated,
#'    sampleinfo = csvfile,
#'    reference = 'SAEC',
#'    chrom = 'chr21')
#' }
#' @export

getCounts <- function(annotpeaks, sampleinfo, reference, chrom = NULL) {
  bamfileslist <- loadBamFiles(sampleinfo)

  if (is.null(chrom) == FALSE) {
    inputgranges <- annotpeaks[[1]][seqnames(annotpeaks[[1]]) == chrom,
                                    ]
  } else {
    inputgranges <- annotpeaks[[1]]
  }

  # Count number of reads overlapping each annotated peak
  countsse <- GenomicAlignments::summarizeOverlaps(features = inputgranges,
                                                   reads = bamfileslist,
                                                   mode = "Union",
                                                   singleEnd = TRUE,
                                                   ignore.strand = TRUE)
  # add column labels
  SummarizedExperiment::colData(countsse) <- DataFrame(sampleinfo[, c(1:4)])
  countsse$sample <- as.factor(countsse$sample)

  countsse$status <- stats::relevel(countsse$sample, reference)
  countssedds <- DESeq2::DESeqDataSet(countsse, design = ~status)

  # Optional filtering out of lowcount regions As part of the DESeq2
  # algorithm, more stringent filtering will be applied subsequently
  # countssedds[ rowSums(counts(countssedds)) > 1, ]

  # get counts referenceized by librarysize
  normcountssedds <- SummarizedExperiment::assay(countssedds, norm = T)

  # get region/peak size
  originaldata <- grangestodataframe(inputgranges)
  regionsize <- originaldata$stop - originaldata$start

  # Calculate RPKM for plotting densities multiply by 10^6 and divide by
  # regions size to get rpkm
  myrpkm <- as.data.frame(normcountssedds[, 1] * 10 ^ 6/regionsize)
  for (i in 2:ncol(normcountssedds)) {
    myrpkm[, i] <- normcountssedds[, i] * 10 ^ 6/regionsize
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

  colnames(SummarizedExperiment::rowData(countssedds)) <-
    gsub("meta.","", colnames(SummarizedExperiment::rowData(countssedds)))


  return(list(regioncounts = countssedds, regioncountstats = statdf,
              regioncountsforplot = forplotdf, reference))
}
