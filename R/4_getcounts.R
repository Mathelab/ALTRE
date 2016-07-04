#' Generates count data for regulatory regions
#'
#' Counts the number of reads in each regulatory region for each sample type --
#' read count is derived from user-input BAM filex, and regions of interest
#' are supplied in a Granges object, ideally output of combineannotatepeaks.R.
#'
#' @param annotpeaks list output from combineannotatepeaks function
#' @param sampleinfo name of CSV file, with complete path, that contains the
#'  following 4 columns for each sample
#'  1) complete filepath;
#'  2) name of bamfiles;
#'  3) name of peak files;
#'  4) name or type of sample;
#'  5) name of replicate
#' @param reference name of sample type to be
#' considered 'reference' in DESeq2 analysis
#' @param chrom optional, only chromosome chrom will be evaluated
#'
#' @return List containing three items:
#' (1) DESeqDataSet: contains count information for all replicates of all samples
#' (2) Matrix: contains number of enhancers and promoters
#'  before and after filtering (if applicable)
#' (3) Data frame for creating a density plot (use function plotgetcounts()
#'
#'
#' @examples
#' \dontrun{
#' TSSannot <- getTSS()
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps=2)
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           mergedistenh = 1500,
#'                                           mergedistprom = 1000 )
#' counts_consPeaks <- getcounts(annotpeaks = consPeaksAnnotated,
#'                               sampleinfo = sampleinfo,
#'                               reference = 'SAEC',
#'                               chrom = 'chr21')
#' }
#' @export

getcounts <- function(annotpeaks, sampleinfo, reference, chrom = NULL) {
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
  countssedds <- DESeq2::DESeqDataSet(countsse, design = ~sample)

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
  myrpkm <- as.data.frame(normcountssedds[, 1] * 10^6/regionsize)
  for (i in 2:ncol(normcountssedds)) {
    myrpkm[, i] <- normcountssedds[, i] * 10^6/regionsize
  }
  # take the log2 so that it is a normalized distribution
  myrpkmlog2 <- log2(as.matrix(myrpkm) + 1)
  colnames(myrpkmlog2) <- unlist(lapply(sampleinfo$sample, as.character))

  #########################################
  # Create stats matrix originaldata is created ~ 10 lines lines above
  colnames(originaldata) <- unlist(lapply(colnames(originaldata), gsub,
                                          pattern = "meta.", replacement = ""))
  enhancernum <- length(which(originaldata$region == "enhancer"))
  promoternum <- length(which(originaldata$region == "promoter"))

  statdf <- data.frame(Num_Enhancers = enhancernum,
                       Num_Promoters = promoternum)

  #########################################
  # Create densityplot
  region <- originaldata$region
  forplotdf <- cbind(myrpkmlog2, as.data.frame(region))

  colnames(SummarizedExperiment::rowData(countssedds)) <-
    gsub("meta.", "", colnames(SummarizedExperiment::rowData(countssedds)))

  return(list(regioncounts = countssedds, regioncountstats = statdf,
              regioncountsforplot = forplotdf))
}
