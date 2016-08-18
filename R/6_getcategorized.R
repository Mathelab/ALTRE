#' Categorize altered peaks as experiment-specific, reference-specific,
#'  or shared
#'
#' @param analysisresults output of countanalysis()
#' @param lfctypespecific log2fold change cutoff (of chromatin accessibility)
#' for type specific enhancers/promoters
#' @param lfcshared log2fold change cutoff (of chromatin accessibility)
#' for shared enhancers/promoters
#' @param pvaltypespecific p-value cutoff (of chromatin accessibility)
#' for type specific enhancers/promoters
#' @param pvalshared p-value cutoff (of chromatin accessibility)
#' for shared enhancers/promoters
#'
#' @examples
#' \dontrun{
#' csvfile <- file.path(dir="yourfilepath", 'sampleinfo.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo )
#' consPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps = 2)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           mergedistenh = 1500,
#'                                           mergedistprom = 1000)
#' counts_consPeaks <- getCounts(annotpeaks = consPeaksAnnotated,
#'                               sampleinfo = sampleinfo,
#'                               reference = 'SAEC')
#' altre_peaks <- countanalysis(counts = counts_consPeaks,
#'                              pval = 0.01,
#'                              lfcvalue = 1)
#' categaltre_peaks <- categAltrePeaks(
#'   altre_peaks,
#'   lfctypespecific = 1.5,
#'   lfcshared = 1.2,
#'   pvaltypespecific = 0.01,
#'   pvalshared = 0.05
#' )
#'}
#' @return list 1) results of countanalysis,
#' 2) Frequences of sample type specific, shared,
#'  and ambiguous regulatory elements
#'
#' @export
 categAltrePeaks <- function(analysisresults,
 	lfctypespecific = 1.5,
 	lfcshared = 1.2,
 	pvaltypespecific = 0.01,
 	pvalshared = 0.05){

   #quick fix
   #names(analysisresults$results)[10:12] <- c("region","A549","SAEC")

   analysisresults <- analysisresults[[1]]

   if (is.data.frame(analysisresults) == FALSE ||
       is.null(analysisresults$padj) ||
       is.null(analysisresults$log2FoldChange)) {
     stop("analysisresults parameter is not in the correct format,
          make sure you are using the output from countanalysis()")
   }

   # Define regions that are more open, less open, or
   # shared
   exptind <- which(!(is.na(analysisresults$padj)) &
                      analysisresults$log2FoldChange > lfctypespecific &
                      analysisresults$padj < pvaltypespecific)
   refind <- which(!(is.na(analysisresults$padj)) &
                     analysisresults$log2FoldChange < -lfctypespecific &
                     analysisresults$padj < pvaltypespecific)
   sharedind <- which((analysisresults$log2FoldChange <= lfcshared &
                         analysisresults$log2FoldChange >= -lfcshared) &
                        (analysisresults$padj >= pvalshared |
                           is.na(analysisresults$padj)))
   ambigind <- setdiff(1:nrow(analysisresults), c(exptind, refind, sharedind))
   REaltrecateg <- rep(NA, nrow(analysisresults))
   REaltrecateg[exptind] <- "Experiment Specific"
   REaltrecateg[refind] <- "Reference Specific"
   REaltrecateg[sharedind] <- "Shared"
   REaltrecateg[ambigind] <- "Ambiguous"


   if (!all.equal(sum(table(REaltrecateg)), nrow(analysisresults))) {
     stop("Categorization failed, some REs are not categorized")
   }

   analysisresults$REaltrecateg <- REaltrecateg

   stats <- data.frame(
     REtype = c("Promoter", "Enhancer"),
     Expt_Specific = c(length(intersect(
       which(analysisresults$region == "promoter"),
       which(analysisresults$REaltrecateg ==
               "Experiment Specific")
     )),
     length(intersect(
       which(analysisresults$region == "enhancer"),
       which(analysisresults$REaltrecateg ==
               "Experiment Specific")
     ))),
     Ref_Specific = c(length(intersect(
       which(analysisresults$region ==  "promoter"),
       which(analysisresults$REaltrecateg == "Reference Specific")
     )),
     length(intersect(
       which(analysisresults$region == "enhancer"),
       which(analysisresults$REaltrecateg == "Reference Specific")
     ))),
     Shared = c(length(intersect(
       which(analysisresults$region == "promoter"),
       which(analysisresults$REaltrecateg == "Shared")
     )),
     length(intersect(
       which(analysisresults$region == "enhancer"),
       which(analysisresults$REaltrecateg == "Shared")
     ))),
     Ambiguous = c(length(intersect(
       which(analysisresults$region == "promoter"),
       which(analysisresults$REaltrecateg == "Ambiguous")
     )), length(intersect(
       which(analysisresults$region == "enhancer"),
       which(analysisresults$REaltrecateg == "Ambiguous")
     )))
   )


   return(list(analysisresults = analysisresults, stats = stats))
} # end function

