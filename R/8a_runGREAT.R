#' Enrichment analysis using GREAT package
#' to identify putative pathways of interest for further
#' investigation
#' @param peaks list, output of categAltrePeaks() function
##' @param peaktype character, "Experiment Specific", "Reference Specific",
##' 	"Ambiguous", "Shared", or "All" (All is default)
#' @param species default hg19
#' @param rule character, "basalPlusExt", "twoClosest", "oneClosest" rule that associates
#' 	genomic regions to genes (default is "basalPlusExt").
#' 	See https://bioconductor.org/packages/release/bioc/html/chipenrich.html for more detail.
#' @param adv_upstream kb, extension to upstream (if rule is basalPlusExt), default 5
#' @param adv_downstream kb, extension to downstream (if rule is basalPlusExt), default 1.0
#' @param adv_span kb, max extension (if rule is basalPlusExt), default 1000.0
#' @param adv_twoDistance kb, max extension (if rule is twoClosest), default 1000.0
#' @param adv_oneDistance kb, max extension (if rule is oneClosest), default 1000.0
#' @param pathway_category character, "GO", "Pathway Data", "Regulatory Motifs",
#'	"Phenotype Data and Human Disease", "Gene Expression", "Gene Families"
#'	(default is "GO")
#' @examples
#' \dontrun{
#' runGREAT(peaks=categaltre_peaks)
#' }
#' @return ways --
#' pathways also annotated with additional information
# run with categaltre_peaks
#' @export
#'
runGREAT <- function(peaks,
                     #	peaktype="All",
                     species = "hg19",
                     rule = "basalPlusExt",
                     adv_upstream = 5.0,
                     adv_downstream = 1.0,
                     adv_span = 1000.0,
                     adv_twoDistance = 1000.0,
                     adv_oneDistance = 1000.0,
                     pathway_category = "GO") {
  # Check that peaktype entry is allowable and grab peaks for analysis
  #  if (is.na(match(peaktype,c("All","Experiment Specific", "Shared", "Ambiguous",
  #	"Reference Specific")))) {
  #	stop("peaktype should be either 'All', 'Experiment Specific', 'Shared',
  #		'Ambiguous', or 'Reference Specific'")
  #  }
  #  if (peaktype == "All") {
  #	mypeaks = as.data.frame(peaks$analysisresults)[,c("chr","start","stop")] }
  #  else {
  #	mypeaks = as.data.frame(peaks$analysisresults)[which(peaks$analysisresults==peaktype),
  #		c("chr","start","stop")]}
  if (is.na(match(rule, c(
    "basalPlusExt", "twoClosest", "oneClosest"
  )))) {
    stop("rule must be either 'basalPlusExt', 'twoClosest', 'oneClosest'")
  }

  mygreat = list()
  for (i in c("Experiment Specific", "Reference Specific", "Shared")) {
    print(paste("Running", i))
    mypeaks = as.data.frame(peaks$analysisresults)[which(peaks$analysisresults$REaltrecateg ==
                                                           i),
                                                   c("chr", "start", "stop")]
    ilabel = gsub(" ", "_", i)
    # Run GREAT
    if (rule == "basalPlusExt") {
      mygreat[[ilabel]] = rGREAT::submitGreatJob(
        mypeaks,
        species = species,
        adv_span = adv_span,
        rule = "basalPlusExt",
        adv_upstream = adv_upstream,
        adv_downstream = adv_downstream,
        request_interval = 20
      )
    }
    if (rule == "twoClosest") {
      mygreat[[ilabel]] = rGREAT::submitGreatJob(
        mypeaks,
        species = species,
        rule = "twoClosest",
        adv_twoDistance = adv_twoDistance,
        request_interval = 20
      )
    }
    if (rule == "oneClosest") {
      mygreat[[ilabel]] = rGREAT::submitGreatJob(
        mypeaks,
        species = species,
        rule = "oneClosest",
        adv_twoDistance = adv_oneDistance,
        request_interval = 20
      )
    }
  } # end looping through peak types
  return(mygreat)
} # end function


#' Enrichment analysis using GREAT package
#' to identify putative pathways of interest for further
#' investigation
#'
#' @param GREATpath output of runGREAT()
#' @param pathway_category character, "GO", "Pathway Data", "Regulatory Motifs",
#'      "Phenotype Data and Human Disease", "Gene Expression", "Gene Families"
#'      (default is "GO")
#' @param enrichcutoff numeric, fold change enrichment cutoff to determine enriched pathways,
#' default is 2
#' @param adjpvalcutoff numeric, Bonferroni adjusted p-value cutoff to determine enriched pathways,
#'	default is 0.05
#' @param adjustby character, "fdr" or "bonferroni", default is "bonferroni"
#' @param test character, "Both" denotes hypergeometric and binomical tests are used to
#' 	determine enriched pathways, "Binom" denotes binomial tests used, "Hyper" denotes
#'	hypergeometric tests are used.  Default is "Binom"
#'
#' @return list of dataframes for enriched pathways - each dataframe in the list
#' represents one pathway type (e.g. "GO Molecular Function")
#' @export
#'
processPathways <- function(GREATpath,
                            pathway_category = "GO",
                            adjustby = "bonferroni",
                            test = "Binom",
                            enrichcutoff = 2,
                            adjpvalcutoff = 0.05) {
  finaloutput = list()
  for (job in names(GREATpath)) {
    if (!is(GREATpath[[job]], "GreatJob")) {
      stop(
        "GREATpath is not a list of 'GreatJob' objects.
        Input should be the output of runGREAT()"
      )
    }

    output = rGREAT::getEnrichmentTables(GREATpath[[job]], category = pathway_category)
    names(output) = gsub(" ", "_", names(output))
    stats = data.frame(Pathway = names(output),
                       NumSig = rep(NA, length(names(output))))

    for (i in names(output)) {
      output[[i]]$Binom_adj_PValue = stats::p.adjust(output[[i]]$"Binom_Raw_PValue",
                                                     method = adjustby)
      output[[i]]$Hyper_adj_PValue = stats::p.adjust(output[[i]]$"Hyper_Raw_PValue",
                                                     method = adjustby)
      if (test == "Both") {
        keepers = base::Reduce(intersect, list(
          which(output[[i]]$Binom_Fold_Enrichment > enrichcutoff),
          which(output[[i]]$Hyper_Fold_Enrichment >
                  enrichcutoff),
          which(output[[i]]$Binom_adj_PValue <=
                  adjpvalcutoff),
          which(output[[i]]$Hyper_adj_PValue <=
                  adjpvalcutoff)
        ))
      }
      else if (test == "Binom") {
        keepers = base::Reduce(intersect, list(
          which(output[[i]]$Binom_Fold_Enrichment > enrichcutoff),
          which(output[[i]]$Binom_adj_PValue <=
                  adjpvalcutoff)
        ))
      }
      else if (test == "Hyper") {
        keepers = base::Reduce(intersect, list(
          which(output[[i]]$Hyper_Fold_Enrichment > enrichcutoff),
          which(output[[i]]$Hyper_adj_PValue <=
                  adjpvalcutoff)
        ))
      }
      else {
        stop("test should be 'Both', 'Binom', or 'Hyper'")
      }
      print(length(keepers))
      output[[i]] = output[[i]][keepers, ]
      output[[i]] = output[[i]][base::order(output[[i]]$Binom_adj_PValue), ]
      stats$NumSig[which(stats$Pathway == i)] =
        length(which(output[[i]]$Binom_adj_PValue <= adjpvalcutoff))
    }
    finaloutput[[job]] = list(Sig_Pathways = output, stats = stats)
    } # end looping through each GREAT job
  return(finaloutput)
} # end function
