#' Determines which regulatory regions are signifigantly altered
#' between sample types
#'
#' Altered regions are those that show differences in chromatin accessibility
#' (using DESeq2 algorithm)
#'
#' @param counts counts for each region
#' @param pval optional, pvalue considered significant (0.05, 0.01, etc.)
#' @param lfcvalue optional, logfold change value considered significant
#' (value reported on a log scale base 2 so log2fold change of 1.5 means
#' difference in peaks increased by 2^1.5 or 2.8)
#'
#' @return list containing:
#' 1) DESeq2 results table
#' 2) some statistics
#' 3) data.frame used for plotting
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
#' consensusPeaksCounts <- getCounts(annotpeaks = consensusPeaksAnnotated,
#'    sampleinfo = csvfile,
#'    reference = 'SAEC',
#'    chrom = 'chr21')
#' alteredPeaks <- countanalysis(counts = consensusPeaksCounts,
#'    pval = 0.01,
#'    lfcvalue = 1)
#' }
#' @export

countanalysis <- function(counts,
                          pval = 0.01,
                          lfcvalue = 1) {
  countsonly <- counts[[1]]
  # run differential analysis
  countsdiffexp <- suppressMessages(DESeq2::DESeq(countsonly))
  # define the threshold
  countsresultsalphalfc <- DESeq2::results(countsdiffexp,
                                           alpha = pval,
                                           lfcThreshold = lfcvalue)

   resultsdataframe <- as.data.frame(countsresultsalphalfc)
   originalgranges <- SummarizedExperiment::rowRanges(countsonly)
   grangesdataframe <- grangestodataframe(originalgranges)

   fulldataframe <- cbind(resultsdataframe, grangesdataframe)
   #Rename some columns
   metacols <- grep("meta",colnames(fulldataframe))
   colnames(fulldataframe)[metacols] <-
     gsub("metab.","", colnames(fulldataframe)[metacols])

  return(list(results <- as.data.frame(fulldataframe), counts[[4]]))
}

