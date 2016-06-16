#' Determines which regulatory regions are signifigantly altered between sample types
#'
#' Altered regions are those that show differences in chromatin accessibility 
#' (using DESeq2 algorithm)
#'
#' @param counts counts for each region
#' @param pval pvalue considered signifigant (0.05, 0.01, etc.).
#' @param lfcvalue logfold change value considered significant (value reported on a
#' log scale base 2 so log2fold change of 1.5 means difference in peaks
#' increased by 2^1.5)
#'
#' @return dataframe containing logfold2 changes and p-values of all regions
#' analyzed
#'
#' @examples
#' dir=system.file("extdata", package="ALTRE", mustWork=TRUE)
#' TSSpath=file.path(dir,"Homosapiens_GRCh37.75_TSS.bed")
#' TSSannot=read.table(TSSpath, header=TRUE)
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' consPeaksAnnotated=CombineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot)
#' counts_consPeaks=getcounts(annotpeaks=consPeaksAnnotated, csvfile=csvfile, reference="SAEC", chrom="chr21")
#' altre_peaks=countanalysis(counts=counts_consPeaks, pval=0.01, lfcvalue=1)
#'
#' @export

countanalysis<-function(counts, pval=0.01, lfcvalue=1){
  countsdds=counts[[1]]
  errortest=try(DESeq2::counts(countsdds), silent=TRUE)
  if (inherits(errortest, 'try-error')==TRUE){
    stop("The input is not a summarized experiment object!")
  }

  #run differential analysis
  countsdiffexp <- suppressMessages(DESeq2::DESeq(countsdds))
  #define the threshold
  countsresultsalphalfc <- DESeq2::results(countsdiffexp, alpha=pval, lfcThreshold=lfcvalue)

  #change row labels to dataframe
  originalgranges=SummarizedExperiment::rowRanges(countsdds)
  grangesdataframe=grangestodataframe(originalgranges)

  #create a dataframe with results and position
  resultsdataframe=as.data.frame(countsresultsalphalfc)
  fulldataframe=cbind(resultsdataframe, grangesdataframe)
  return(list(altre_peaks=fulldataframe))
}

