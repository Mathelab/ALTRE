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
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           mergedistenh = 1500,
#'                                           mergedistprom = 1000)
#' counts_consPeaks <- getcounts(annotpeaks = consPeaksAnnotated,
#'                               csvfile = csvfile,
#'                               reference = 'SAEC')
#' altre_peaks <- countanalysis(counts=counts_consPeaks,
#'                              pval=0.01,
#'                              lfcvalue=1)
#' }
#' @export

countanalysis <- function(counts, pval = 0.01,
                          lfcvalue = 1) {
  countsdds <- counts[[1]]
  errortest <- try(DESeq2::counts(countsdds),
                   silent = TRUE)
  if (inherits(errortest, "try-error") == TRUE) {
    stop("The input is not a summarized experiment object!")
  }

  # run differential analysis
  countsdiffexp <- suppressMessages(DESeq2::DESeq(countsdds))
  # define the threshold
  countsresultsalphalfc <- DESeq2::results(countsdiffexp,
                                           alpha = pval,
                                           lfcThreshold = lfcvalue)

   resultsdataframe <- as.data.frame(countsresultsalphalfc)
   originalgranges <- SummarizedExperiment::rowRanges(countsdds)
   grangesdataframe <- grangestodataframe(originalgranges)

   fulldataframe <- cbind(resultsdataframe,
	grangesdataframe)
   #Rename some columns
   metacols=grep("meta",colnames(fulldataframe))
   colnames(fulldataframe)[metacols]=gsub("metab.","",colnames(fulldataframe)[metacols])   

  return(list(results = as.data.frame(fulldataframe)))  #,
   #           stats = stats,
  #             dftoplot = list(toplot = toplot,
                             # pval = pval,
                             # lfcvalue = lfcvalue)))
}

