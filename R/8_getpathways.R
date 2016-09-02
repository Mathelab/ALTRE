#' Enrichment analysis to identify putative pathways of interest for further
#' investigation
#'
#' Determine which pathways are overrepresented in
#' altered TSS-proximal and TSS-distal regions. Pathways are determined by linking
#' the TSS-proximal/distal regions to the nearest gene, and then linking genes to pathways
#' using the gene ontology database. The 'gene' argument limits how few genes
#' a pathway can contain, while the 'offspring' argument limits how many
#' offspring a pathway can contain. Pathways with low gene counts are less
#' reliable (often false positives), while pathways with many offspring are
#' vague and unlikely to be of much use -- the enrichment of their more precise
#'  offspring is the more interesting question.
#'
#' ************CHANGE ANALYSISRESULTS TO ALTREPEAKS!!!!******************
#' ****************CHANGE REGIONSUBSET TO REGIONTYPE!!!!*****************
#' @param analysisresults Results from analysis of counts, categaltre_peaks.
#' @param ontoltype One of three categories: 'MF' (molecular function),
#'  'CC' (cellular component), 'BP' (biological process).
#' @param enrichpvalfilt Adjusted pval for enrichment to filter on
#'  (adjusted for multiple testing).
#' @param lfctypespecific Log2fold change (of chromatin accessibility)
#' for type specific TSS-proximal and TSS-distal regions.
#' @param lfcshared Log2fold change (of chromatin accessibility) for
#'  shared TSS-proximal and TSS-distal regions.
#' @param pvaltypespecific P-value (of chromatin accessibility) for
#'  type specific TSS-proximal and TSS-distal regions.
#' @param pvalshared P-value (of chromatin accessibility) for
#' shared TSS-proximal and TSS-distal regions.
#' @param genes Minimum number of genes allowable in a pathway.
#' @param offspring Maximum number of offspring allowable in a pathway.
#' @param regionsubset  'TSS-proximal' or 'TSS-distal'. Default is "TSS-distal".
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
#' comparePeaksAnalysisResults <- comparePeaksAltre(alteredPeaksCategorized)
#' MFenrich <- pathenrich(analysisresults = alteredPeaksCategorized,
#'    ontoltype = 'MF',
#'    enrichpvalfilt = 0.99)}
#' @return dataframe identifying p-values for enriched pathways --
#' pathways also annotated with additional information
#'
#'
#' @export

pathenrich <- function(analysisresults,
                       ontoltype = "MF",
                       enrichpvalfilt = 0.01,
                       lfctypespecific = 1.5,
                       lfcshared = 1.2,
                       pvaltypespecific = 0.01,
                       pvalshared = 0.05,
                       genes = 20,
                       offspring = 300,
                       regionsubset = "TSS-distal") {

  analysisResults <- analysisresults$analysisresults

  if (is.data.frame(analysisResults) == FALSE) {
    stop("analysisresults parameter is not in the correct format,
         make sure you are using the output from countanalysis()")
  }

  #make sure the sample names are the user-entered names 
  reference <- analysisresults[[3]] 
  allSamples <- colnames(analysisResults)[12:length(analysisResults)-1]
  nonreference<-allSamples[which(!(allSamples %in% reference))]
  referenceSpecific <- paste0(reference, "SpecificByIntensity")
  experimentSpecific <- paste0(nonreference, "SpecificByIntensity")
  

  if (regionsubset == "all") {
    newanalysisresults <- analysisResults
  } else {
    newanalysisresults <- analysisResults[analysisResults$region ==
                                                regionsubset, ]
  }

  up <- newanalysisresults[which(newanalysisresults$REaltrecateg ==
                                     experimentSpecific), ]
  down <- newanalysisresults[which(newanalysisresults$REaltrecateg ==
                                       referenceSpecific), ]
  shared <- newanalysisresults[which(newanalysisresults$REaltrecateg ==
                                       "Shared"), ]
  all <- rbind(up, down, shared)
  subsets <- list(up, down, shared, newanalysisresults)
  names(subsets) <- c(experimentSpecific, referenceSpecific, "shared", "all")

  message("finding expt-specific...")
  if (nrow(subsets[[experimentSpecific]]) == 0) {
    expt <- as.data.frame("No REs higher in experiment group")
  } else {
    expt <- rundose(set = subsets[[experimentSpecific]],
                    background = subsets[["all"]],
                    log2FoldChange = lfctypespecific,
                    ontoltype = ontoltype,
                    pvalfilt = enrichpvalfilt,
                    genes = genes,
                    offspring = offspring)
    if (nrow(shared) == 0) {
      expt <- as.data.frame("No enrichment found for experiment REs")
    }
  }

  message("finding reference-specific...")
  if (nrow(subsets[[referenceSpecific]]) == 0) {
    reference <- as.data.frame("No REs higher in reference group")
  } else {
    reference <- rundose(set = subsets[[referenceSpecific]],
                         background = subsets[["all"]],
                         log2FoldChange = lfctypespecific,
                         ontoltype = ontoltype,
                         pvalfilt = enrichpvalfilt,
                         genes = genes,
                         offspring = offspring)
    if (nrow(reference) == 0) {
      reference <- as.data.frame("No enrichment found for reference REs")
    }
  }

  message("finding shared...")
  if (nrow(subsets[["shared"]]) == 0) {
    shared <- as.data.frame("No shared REs")
  } else {
    shared <- rundose(set = subsets[["shared"]],
                      background = subsets[["all"]],
                      log2FoldChange = lfcshared,
                      ontoltype = ontoltype,
                      pvalfilt = enrichpvalfilt,
                      genes = genes,
                      offspring = offspring)
    print(paste("Number of rows:", nrow(shared)))
    if (nrow(shared) == 0) {
      shared <- as.data.frame("No enrichment found for shared REs")
    }
  }
  

 
  allthree <- list(expt = as.data.frame(expt),
                   reference = as.data.frame(reference),
                   shared = as.data.frame(shared))
  names(allthree) <- c(experimentSpecific, referenceSpecific, "shared")

  return(allthree)
}
