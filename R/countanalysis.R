#' Determines which regulatory regions are signifigantly altered between sample types
#'
#' Altered regions are those that show differences in chromatin accessibility
#' (using DESeq2 algorithm)
#'
#' @param counts counts for each region
#' @param pval optional, pvalue considered significant (0.05, 0.01, etc.)
#' @param lfcvalue optional, logfold change value considered significant (value reported on a
#' log scale base 2 so log2fold change of 1.5 means difference in peaks
#' increased by 2^1.5)
#'
#' @return list containing:
#' 	1) DESeq2 results table
#'	2) some statistics
#'	3) data.frame used for plotting
#'
#' @examples
#' TSSannot <- getTSS()
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot,
#'	regionspecific=TRUE,mergedistenh=1500,mergedistprom=1000 )
#' counts_consPeaks <-getcounts(annotpeaks=consPeaksAnnotated, csvfile=csvfile,
#' reference="SAEC", chrom="chr21")
#' altre_peaks=countanalysis(counts=counts_consPeaks, pval=0.01, lfcvalue=1)
#'
#' @export

countanalysis <- function(counts,
                          pval = 0.01,
                          lfcvalue = 1) {
  countsdds = counts[[1]]
  errortest = try(DESeq2::counts(countsdds), silent = TRUE)
  if (inherits(errortest, 'try-error') == TRUE) {
    stop("The input is not a summarized experiment object!")
  }

  #run differential analysis
  countsdiffexp <- suppressMessages(DESeq2::DESeq(countsdds))
  #define the threshold
  countsresultsalphalfc <-
    DESeq2::results(countsdiffexp, alpha = pval, lfcThreshold = lfcvalue)

  #create a dataframe with results and position and filter by pval and lfcvalue
  resultsdataframe = as.data.frame(countsresultsalphalfc)
  keepers = intersect(which(resultsdataframe$padj < pval),
                      which(abs(resultsdataframe$log2FoldChange) > lfcvalue))

  #change row labels to dataframe
  originalgranges = SummarizedExperiment::rowRanges(countsdds)
  grangesdataframe = grangestodataframe(originalgranges)

  #create a dataframe with results and position
  fulldataframe = cbind(resultsdataframe[keepers, ], grangesdataframe[keepers, ])

  #######################
  # Get some stats
  tempstats = table(fulldataframe$meta.region[which(fulldataframe$log2FoldChange <
                                                      0)])
  downproms = as.numeric(tempstats["promoter"])
  downenh = as.numeric(tempstats["enhancer"])
  tempstats = table(fulldataframe$meta.region[which(fulldataframe$log2FoldChange >
                                                      0)])
  upproms = as.numeric(tempstats["promoter"])
  upenh = as.numeric(tempstats["enhancer"])

  stats = matrix(c(upproms, downproms,
                   upenh, downenh), nrow = 2)
  colnames(stats) = c("Promoters", "Enhancers")
  rownames(stats) = c("UP", "DOWN")

  #######################
  # Produce a volcano plot and a barplot showing differential REs
  # Find NAs to remove them, thereby avoiding printing warnings
  nonas = which(!is.na(resultsdataframe$padj))
  toplot = resultsdataframe[nonas, ]
  keepers = intersect(which(toplot$padj < pval),
                      which(abs(toplot$log2FoldChange) > lfcvalue))
  toplot$col = rep("non-altered", nrow(toplot))
  toplot$col[keepers] = "altered"

  return(list(
    results = fulldataframe,
    stats = stats,
    dftoplot = list(toplot=toplot,pval=pval,lfcvalue=lfcvalue)
  ))
}
