#' Comparison of methods for identifying altered regulatory regions
#'
#' Creates a table to compare two methods of identifying altered regulatory
#' regions, one based on peak intensity, the other on peak presence as
#'  determined by hotspot calling algorithms.
#'
#' @param analysisresults analysisresults of countanalysis.
#' @param reference cell type to be considered "reference" or "reference" to which
#' other cell types will be compared
#' @param lfctypespecific log2fold change for type specific TSS-dists/TSS-proxs
#' @param lfcshared log2fold chance for shared TSS-dists/TSS-proxs
#' @param pvaltypespecific p-value for type specific TSS-dists/TSS-proxs
#' @param pvalshared p-value for shared TSS-dists/TSS-proxs
#'
#' @return matrix comparing the two methods of identifying altered regulatory
#' regions, one based on peak intensity, the other on peak presence as
#' determined by hotspot calling algorithms.
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
#' alteredPeaks <- countanalysis(counts = consensusPeaksCounts,
#'    pval = 0.01,
#'    lfcvalue = 1)
#' alteredPeaksCategorized <- categAltrePeaks(alteredPeaks,
#'    lfctypespecific = 1.5,
#'    lfcshared = 1.2,
#'    pvaltypespecific = 0.01,
#'    pvalshared = 0.05)
#' comparePeaksAnalysisResults <- comparePeaksAltre(alteredPeaksCategorized, reference = "SAEC")
#' }
#' @export
#'
comparePeaksAltre <- function(analysisresults,
	reference,
	lfctypespecific=1.5,
	lfcshared=1.2,
	pvaltypespecific=0.01,
	pvalshared=0.05){

  analysisresults <- analysisresults[[1]]

  if (!is.data.frame(analysisresults)) {
    stop("Make sure the output of the analysis is from categAltrePeaks() function")
  }

  if (!reference %in% colnames(analysisresults)) {
	stop("Make sure the reference sample exists!")
  }

  samplenames <- colnames(analysisresults)[11:(ncol(analysisresults) - 1)]
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
  two <- rep(c("TSS-dists", "TSS-proxs", "all"), 3)

  matrixrownames <- paste(one, two)
  rownames(analysisresultsmatrix) <- matrixrownames
  colnames(analysisresultsmatrix) <- c("intensity", "peak")

  # create a matrix for my analysisresults which says how mant of one-specific
  # total, two-specific TSS-proxs, etc.

  analysisresultsmatrix[1, 1] <-
    length(which(analysisresults$log2FoldChange > lfctypespecific &
                   analysisresults$padj < pvaltypespecific &
                   analysisresults$region == "TSS-distal"))
  analysisresultsmatrix[2, 1] <-
    length(which(analysisresults$log2FoldChange > lfctypespecific &
                   analysisresults$padj <  pvaltypespecific &
                   analysisresults$region == "TSS-proximal"))
  analysisresultsmatrix[3, 1] <-
    length(which(analysisresults$log2FoldChange > lfctypespecific &
                   analysisresults$padj <  pvaltypespecific))
  analysisresultsmatrix[4, 1] <-
    length(which(analysisresults$log2FoldChange < -lfctypespecific &
                   analysisresults$padj <  pvaltypespecific &
                   analysisresults$region ==  "TSS-distal"))
  analysisresultsmatrix[5, 1] <-
    length(which(analysisresults$log2FoldChange <  -lfctypespecific &
                   analysisresults$padj <   pvaltypespecific &
                   analysisresults$region ==  "TSS-proximal"))
  analysisresultsmatrix[6, 1] <-
    length(which(analysisresults$log2FoldChange <   -lfctypespecific &
                   analysisresults$padj <   pvaltypespecific))
  analysisresultsmatrix[7, 1] <-
    length(which((analysisresults$log2FoldChange >= -lfcshared &
                    analysisresults$log2FoldChange <= lfcshared) &
                   (analysisresults$padj >=   pvalshared |
                      is.na(analysisresults$padj)) &
                   analysisresults$region == "TSS-distal"))
  analysisresultsmatrix[8, 1] <-
    length(which((analysisresults$log2FoldChange >= -lfcshared &
                    analysisresults$log2FoldChange <= lfcshared) &
                   (analysisresults$padj >= pvalshared |
                      is.na(analysisresults$padj)) &
                   analysisresults$region == "TSS-proximal"))
  analysisresultsmatrix[9, 1] <-
    length(which((analysisresults$log2FoldChange >= -lfcshared &
                    analysisresults$log2FoldChange <= lfcshared) &
                   analysisresults$padj >= pvalshared))

  typespecsub <- analysisresults[, 11:ncol(analysisresults)]

  referencesub <- as.matrix(analysisresults[,  reference])
  nonreferencesub <- as.matrix(analysisresults[,  nonreference])

  analysisresultsmatrix[1, 2] <-
    tssdistnonreference <-
    length(which(apply(nonreferencesub,
                       1,
                       function(x) {any(is.na(x) == FALSE)}) &
                   apply(referencesub, 1, function(x) {any(is.na(x))}) &
                   analysisresults$region == "TSS-distal"))
  analysisresultsmatrix[2, 2] <-
    tssproxnonreference <-
    length(which(apply(nonreferencesub,
                       1,
                       function(x) {any(is.na(x) == FALSE)}) &
                   apply(referencesub, 1, function(x) {any(is.na(x))}) &
                   analysisresults$region == "TSS-proximal"))
  analysisresultsmatrix[3, 2] <-
    allnonreference <-
    length(which(apply(nonreferencesub,
                       1,
                       function(x) {any(is.na(x) == FALSE)}) &
                   apply(referencesub,
                         1,
                         function(x) {any(is.na(x))})))
  analysisresultsmatrix[4, 2] <-
    tssdistreference <-
    length(which(apply(referencesub,
                       1,
                       function(x) { any(is.na(x) == FALSE)}) &
                   apply(nonreferencesub,
                         1,
                         function(x) { all(is.na(x))}) &
                   analysisresults$region == "TSS-distal"))
  analysisresultsmatrix[5, 2] <-
    tssproxreference <-
    length(which(apply(referencesub,
                       1,
                       function(x) { any(is.na(x) == FALSE)}) &
                   apply(nonreferencesub,
                         1,
                         function(x) { all(is.na(x))}) &
                   analysisresults$region == "TSS-proximal"))
  analysisresultsmatrix[6, 2] <-
    allreference <-
    length(which(apply(referencesub,
                       1,
                       function(x) { any(is.na(x) == FALSE)}) &
                   apply(nonreferencesub,
                         1,
                         function(x) { all(is.na(x))})))
  analysisresultsmatrix[7, 2] <-
    tssdistshared <-
    length(which(apply(referencesub,
                       1, function(x) {
                         any(is.na(x) == FALSE)
                       }) & apply(nonreferencesub, 1, function(x) {
                         any(is.na(x) == FALSE)
                       }) & analysisresults$region == "TSS-distal"))
  analysisresultsmatrix[8, 2] <-
    tssproxshared <-
    length(which(apply(referencesub,
                       1, function(x) {
                         any(is.na(x) == FALSE)
                       }) & apply(nonreferencesub, 1, function(x) {
                         any(is.na(x) == FALSE)
                       }) & analysisresults$region == "TSS-proximal"))
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

