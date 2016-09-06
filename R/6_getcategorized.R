#' Categorize altered peaks as experiment-specific, reference-specific,
#'  or shared
#'
#' @param analysisresults output of countanalysis()
#' @param lfctypespecific log2fold change cutoff (of chromatin accessibility)
#' for type specific TSS-distals/TSS-proximals
#' @param lfcshared log2fold change cutoff (of chromatin accessibility)
#' for shared TSS-distals/TSS-proximals
#' @param pvaltypespecific p-value cutoff (of chromatin accessibility)
#' for type specific TSS-distals/TSS-proximals
#' @param pvalshared p-value cutoff (of chromatin accessibility)
#' for shared TSS-distals/TSS-proximals
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
#' alteredPeaks <- countanalysis(counts = consensusPeaksCounts,
#'    pval = 0.01,
#'    lfcvalue = 1)
#' alteredPeaksCategorized <- categAltrePeaks(alteredPeaks,
#'    lfctypespecific = 1.5,
#'    lfcshared = 1.2,
#'    pvaltypespecific = 0.01,
#'    pvalshared = 0.05)
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
    
    #Make sure to names things are from the user-entered sample names 
    reference <- analysisresults[[2]] 
    allSamples <- colnames(analysisresults[[1]])[11:length(analysisresults[[1]])]
    experimentSpecific <- allSamples[which(!(allSamples %in% reference))]
    referenceSpecific <- paste0(reference, "SpecificByIntensity")
    experimentSpecific <- paste0(experimentSpecific, "SpecificByIntensity")
    
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
   REaltrecateg <- REaltrecategplot <- rep(NA, nrow(analysisresults))
   REaltrecategplot[exptind] <- experimentSpecific
   REaltrecategplot[refind] <- referenceSpecific
   REaltrecateg[exptind] <- "ExperimentSpecificByIntensity"
   REaltrecateg[refind] <- "ReferenceSpecificByIntensity"
   REaltrecategplot[sharedind] <- REaltrecateg[sharedind] <- "Shared"
   REaltrecategplot[ambigind] <- REaltrecateg[ambigind] <- "Ambiguous"


   if (!all.equal(sum(table(REaltrecateg)), nrow(analysisresults))) {
     stop("Categorization failed, some REs are not categorized")
   }

   analysisresults$REaltrecateg <- REaltrecateg
   analysisresults$REaltrecategplot <- REaltrecategplot

   stats <- data.frame(
     REtype = c("TSS-proximal", "TSS-distal"),
     Expt_Specific = c(length(intersect(
       which(analysisresults$region == "TSS-proximal"),
       which(analysisresults$REaltrecategplot ==
               experimentSpecific)
     )),
     length(intersect(
       which(analysisresults$region == "TSS-distal"),
       which(analysisresults$REaltrecategplot ==
               experimentSpecific)
     ))),
     Ref_Specific = c(length(intersect(
       which(analysisresults$region ==  "TSS-proximal"),
       which(analysisresults$REaltrecategplot == referenceSpecific)
     )),
     length(intersect(
       which(analysisresults$region == "TSS-distal"),
       which(analysisresults$REaltrecategplot == referenceSpecific)
     ))),
     Shared = c(length(intersect(
       which(analysisresults$region == "TSS-proximal"),
       which(analysisresults$REaltrecategplot == "Shared")
     )),
     length(intersect(
       which(analysisresults$region == "TSS-distal"),
       which(analysisresults$REaltrecategplot == "Shared")
     ))),
     Ambiguous = c(length(intersect(
       which(analysisresults$region == "TSS-proximal"),
       which(analysisresults$REaltrecategplot == "Ambiguous")
     )), length(intersect(
       which(analysisresults$region == "TSS-distal"),
       which(analysisresults$REaltrecategplot == "Ambiguous")
     )))
   )

    listToReturn <-list(analysisresults = analysisresults, stats = stats, reference)
    names(listToReturn)[3] <- c("reference")
    
   return(listToReturn)
} # end function

