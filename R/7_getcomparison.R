#' Comparison of methods for identifying altered regulatory regions
#'
#' Creates a table to compare two methods of identifying altered regulatory regions
#' – one based on peak intensity, the other on peak presence as determined by
#' hotspot calling algorithms.
#'
#' @param analysisresults analysisresults of countanalysis.
#' @param samplenames vector of sample types
#' @param reference cell type to be considered "reference" or "reference" to which
#' other cell types will be compared
#' @param lfctypespecific log2fold change for type specific enhancers/promoters
#' @param lfcshared log2fold chance for shared enhancers/promoters
#' @param pvaltypespecific p-value for type specific enhancers/promoters
#' @param pvalshared p-value for shared enhancers/promoters
#'
#' @return matrix comparing the two methods of identifying altered regulatory
#' regions – one based on peak intensity, the other on peak presence as
#' determined by hotspot calling algorithms.
#' @examples
#' \dontrun{
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps = 2)
#' plotConsensusPeaks(samplepeaks = consPeaks)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                            TSS = TSSannot,
#'                                            merge = TRUE,
#'                                            regionspecific = TRUE,
#'                                            mergedistenh = 1500,
#'                                            mergedistprom = 1000 )
#' counts_consPeaks <- getcounts(annotpeaks = consPeaksAnnotated,
#'                               csvfile = csvfile,
#'                               reference = 'SAEC')
#' altre_peaks <- countanalysis(counts = counts_consPeaks,
#'                              pval = 0.01,
#'                              lfcvalue = 1)
#' categaltre_peaks <- categAltrePeaks(altre_peaks,
#'                                     lfctypespecific = 1.5,
#'                                     lfcshared = 1.2,
#'                                     pvaltypespecific = 0.01,
#'                                     pvalshared = 0.05)
#' analysisresults <- comparePeaksAltre(categaltre_peaks, reference= 'SAEC')
#' }
#' @export
#'
comparePeaksAltre <- function(analysisresults,
	samplenames,
	reference,
	lfctypespecific=1.5,
	lfcshared=1.2,
	pvaltypespecific=0.01,
	pvalshared=0.05){

  analysisresults = analysisresults[[1]]

  if (!is.data.frame(analysisresults) ||
      !all.equal(
        colnames(analysisresults),
        c(
          "baseMean",
          "log2FoldChange",
          "lfcSE",
          "stat",
          "pvalue",
          "padj",
          "chr",
          "start",
          "stop",
          "region",
          "A549",
          "SAEC",
          "REaltrecateg"
        )
      )) {
    stop("Make sure the output of the analysis is from categAltrePeaks() function")
  }

  # removed this bc prone to errors: samplenames <-
  #colnames(analysisresults[,11:ncol(analysisresults)])
  analysisresultsmatrix <-  matrix(nrow = 9, ncol = 2)
  nonreference <-  samplenames[!(samplenames %in% reference)]

  if (length(nonreference) == 1) {
    string <- nonreference
    nonreferencestr <- nonreference
  } else {
    string <- c()
    for (i in 1:length(nonreference)) {
      or <- "and/or"
      if (i == length(nonreference)) {
        string <- paste(string, or, nonreference[i], sep = " ")
      } else {
        string <- paste(string, nonreference[i], sep = ", ")
      }
    }
    nonreferencestr <- substr(string, 3, nchar(string))
  }

  one <- c(rep(nonreferencestr, 3), rep(reference, 3), rep("Shared", 3))
  two <- rep(c("enhancers", "promoters", "all"), 3)

  matrixrownames <- paste(one, two)
  rownames(analysisresultsmatrix) <- matrixrownames
  colnames(analysisresultsmatrix) <- c("intensity", "peak")

  # create a matrix for my analysisresults which says how mant of one-specific
  # total, two-specific promoters, etc.

  analysisresultsmatrix[1, 1] <-
    length(which(analysisresults$log2FoldChange > lfctypespecific &
                   analysisresults$padj < pvaltypespecific &
                   analysisresults$region == "enhancer"))
  analysisresultsmatrix[2, 1] <-
    length(which(analysisresults$log2FoldChange > lfctypespecific &
                   analysisresults$padj <  pvaltypespecific &
                   analysisresults$region == "promoter"))
  analysisresultsmatrix[3, 1] <-
    length(which(analysisresults$log2FoldChange > lfctypespecific &
                   analysisresults$padj <  pvaltypespecific))
  analysisresultsmatrix[4, 1] <-
    length(which(analysisresults$log2FoldChange < -lfctypespecific &
                   analysisresults$padj <  pvaltypespecific &
                   analysisresults$region ==  "enhancer"))
  analysisresultsmatrix[5, 1] <-
    length(which(analysisresults$log2FoldChange <  -lfctypespecific &
                   analysisresults$padj <   pvaltypespecific &
                   analysisresults$region ==  "promoter"))
  analysisresultsmatrix[6, 1] <-
    length(which(analysisresults$log2FoldChange <   -lfctypespecific &
                   analysisresults$padj <   pvaltypespecific))
  analysisresultsmatrix[7, 1] <-
    length(which((analysisresults$log2FoldChange >= -lfcshared &
                    analysisresults$log2FoldChange <= lfcshared) &
                   (analysisresults$padj >=   pvalshared |
                      is.na(analysisresults$padj)) &
                   analysisresults$region == "enhancer"))
  analysisresultsmatrix[8, 1] <-
    length(which((analysisresults$log2FoldChange >= -lfcshared &
                    analysisresults$log2FoldChange <= lfcshared) &
                   (analysisresults$padj >= pvalshared |
                      is.na(analysisresults$padj)) &
                   analysisresults$region == "promoter"))
  analysisresultsmatrix[9, 1] <-
    length(which((analysisresults$log2FoldChange >= -lfcshared &
                    analysisresults$log2FoldChange <= lfcshared) &
                   analysisresults$padj >= pvalshared))

  typespecsub <- analysisresults[, 11:ncol(analysisresults)]

  referencesub <- as.matrix(analysisresults[,  reference])
  nonreferencesub <- as.matrix(analysisresults[,  nonreference])

  analysisresultsmatrix[1, 2] <-
    enhancernonreference <-
    length(which(apply(nonreferencesub,
                       1,
                       function(x) {any(is.na(x) == FALSE)}) &
                   apply(referencesub, 1, function(x) {any(is.na(x))}) &
                   analysisresults$region == "enhancer"))
  analysisresultsmatrix[2, 2] <-
    promoternonreference <-
    length(which(apply(nonreferencesub,
                       1,
                       function(x) {any(is.na(x) == FALSE)}) &
                   apply(referencesub, 1, function(x) {any(is.na(x))}) &
                   analysisresults$region == "promoter"))
  analysisresultsmatrix[3, 2] <-
    allnonreference <-
    length(which(apply(nonreferencesub,
                       1,
                       function(x) {any(is.na(x) == FALSE)}) &
                   apply(referencesub,
                         1,
                         function(x) {any(is.na(x))})))
  analysisresultsmatrix[4, 2] <-
    enhancerreference <-
    length(which(apply(referencesub,
                       1,
                       function(x) { any(is.na(x) == FALSE)}) &
                   apply(nonreferencesub,
                         1,
                         function(x) { all(is.na(x))}) &
                   analysisresults$region == "enhancer"))
  analysisresultsmatrix[5, 2] <-
    promoterreference <-
    length(which(apply(referencesub,
                       1,
                       function(x) { any(is.na(x) == FALSE)}) &
                   apply(nonreferencesub,
                         1,
                         function(x) { all(is.na(x))}) &
                   analysisresults$region == "promoter"))
  analysisresultsmatrix[6, 2] <-
    allreference <-
    length(which(apply(referencesub,
                       1,
                       function(x) { any(is.na(x) == FALSE)}) &
                   apply(nonreferencesub,
                         1,
                         function(x) { all(is.na(x))})))
  analysisresultsmatrix[7, 2] <-
    enhancershared <-
    length(which(apply(referencesub,
                       1, function(x) {
                         any(is.na(x) == FALSE)
                       }) & apply(nonreferencesub, 1, function(x) {
                         any(is.na(x) == FALSE)
                       }) & analysisresults$region == "enhancer"))
  analysisresultsmatrix[8, 2] <-
    promotershared <-
    length(which(apply(referencesub,
                       1, function(x) {
                         any(is.na(x) == FALSE)
                       }) & apply(nonreferencesub, 1, function(x) {
                         any(is.na(x) == FALSE)
                       }) & analysisresults$region == "promoter"))
  analysisresultsmatrix[9, 2] <-
    allshared <-
    length(which(apply(referencesub,
                       1, function(x) {
                         any(is.na(x) == FALSE)
                       }) & apply(nonreferencesub, 1, function(x) {
                         any(is.na(x) == FALSE)
                       })))

  return(list(analysisresultsmatrix))
}

