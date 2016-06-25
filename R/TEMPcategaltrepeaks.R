#' #' Categorize altered peaks as experiment-specific, reference-specific, or shared
#' #'
#' #' @param analysisresults results from analysis of countanalysis()
#' #' @param lfctypespecific log2fold change (of chromatin accessibility) for type specific enhancers/promoters
#' #' @param lfcshared log2fold change (of chromatin accessibility) for shared enhancers/promoters
#' #' @param pvaltypespecific p-value (of chromatin accessibility) for type specific enhancers/promoters
#' #' @param pvalshared p-value (of chromatin accessibility) for shared enhancers/promoters
#' #'
#' #' @examples
#' #' \dontrun{
#' #' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' #' csvfile <- file.path(dir, 'lung.csv')
#' #' samplePeaks <- loadPeaks(csvfile)
#' #' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' #' TSSannot <- getTSS()
#' #' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot,merge=TRUE,
#' #'\tregionspecific=TRUE,mergedistenh=1500,mergedistprom=1000 )
#' #' # Need to run getcounts on all chromosomes
#' #' counts_consPeaks=getcounts(annotpeaks=consPeaksAnnotated, csvfile=csvfile, reference='SAEC')
#' #' altre_peaks=countanalysis(counts=counts_consPeaks, pval=0.01, lfcvalue=1)
#' #' categaltre_peaks=categAltrePeaks(altrepeaks, lfctypespecific=1.5,lfcshared=1.2,
#' #'\tpvaltypespecific=0.01,pvalshared=0.05)
#' #'}
#' #' @return list 1) stats (number of experiment-specific, reference-specific, and shared REs);
#' #'\t\t2) a volcano plot
#' #'
#' #' @export
#'
#' pathenrich<-function(analysisresults,
#' \tlfctypespecific=1.5,
#' \tlfcshared=1.2,
#' \tpvaltypespecific=0.01,
#' \tpvalshared=0.05){
#'
#'   analysisresults=analysisresults[[1]]
#'
#'   if (is.data.frame(analysisresults)==FALSE) {
#' \tstop('analysisresults parameter is not in the correct format, make sure you are using the output from countanalysis()')
#'   }
#'
#'   # Define regions that are more open, less open, or shared
#'   exptind = which(!(is.na(analysisresults$padj)) &
#'         analysisresults$log2FoldChange > lfctypespecific &
#'         analysisresults$padj < pvaltypespecific)
#'   refind=which(!(is.na(analysisresults$padj)) &
#'         analysisresults$log2FoldChange < -lfctypespecific &
#'         analysisresults$padj < pvaltypespecific)
#'   sharedind=which((analysisresults$log2FoldChange <= lfcshared &
#'         analysisresults$log2FoldChange >= -lfcshared) &
#'         (analysisresults$padj >= pvalshared | is.na(analysisresults$padj)))
#'   ambigind=setdiff(1:nrow(analysisresults),c(exptind,refind,sharedind))
#'
#'   mycateg=rep(NA,nrow(analysisresults))
#'   mycateg[exptind]='Experiment Specific'
#'   mycateg[refind]='Reference Specific'
#'   mycateg[sharedind]='Shared'
#'   mycateg[ambigind]='Ambiguous'
#'
#'   if(length(which(is.na(mcateg)))>0) {
#' \tstop('Categorization failed, some REs are not categorized')
#'   }
#'   else {
#' \tanalysisresults$altrecateg=mycateg
#' \treturn(analysisresults)
#'   }
#'
#' }
