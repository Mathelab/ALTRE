#' Given a GRangeslist object with peaks for each samples, determine the consensus peaks (found in at least N replicates, where N is input by the user) for each sample type
#'
#' @param samplepeaks A GRangesList object comprising one GRanges object (peaks) for each sample (output of loadPeaks() function)

#' @param minreps minimum number of replicate samples that a peak should be contained in to be called as a consensus peak. This cutoff will be applied to both samples.
#'
#' @return a list comprising:
#' 	1) a GRangeslist with one GRange for each sample type which contains consensus peaks
#'	2) a summary statistic table
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#'}
#' @export

getConsensusPeaks <- function(samplepeaks, minreps) {
  if (class(samplepeaks) != "GRangesList")
    stop("Peaks must be a GRangesList Object")

  sampnames = names(samplepeaks)
  sampletypes = sort(unique(gsub("_.*", "", sampnames)))

  conspeaks = GRangesList()
  conspeaks_stats = list()

  for (mytype in order(sampletypes)) {
    mytypepeaks <- samplepeaks[grep(sampletypes[mytype], sampnames)]

    # Concatenate all peaks pertaining to the same sample type and merge peaks:
    allregregions = c(mytypepeaks[[1]])
    for (i in 2:length(mytypepeaks)) {
      allregregions = c(allregregions, mytypepeaks[[i]])
    }
    reducedallregregions = reduce(allregregions)

    # For each reduced peak, determine whether it was present in each sample type
    for (i in 1:length(mytypepeaks)) {
      typespecific = findOverlaps(reducedallregregions, mytypepeaks[[i]])
      newdataframe = data.frame(i = matrix(nrow = length(reducedallregregions)))
      newdataframe[queryHits(typespecific), 1] = "present"
      values(reducedallregregions) = cbind(values(reducedallregregions), newdataframe)
    }
    colnames(values(reducedallregregions)) = names(mytypepeaks)

    # Find regions that are present in at least N replicates (from user in put minreps)
    reducedallregionsdata = grangestodataframe(reducedallregregions)
    applymatrix = as.matrix(reducedallregionsdata[4:ncol(reducedallregionsdata)])
    keepers = which(apply(applymatrix, 1, function(x)
      length(which(x == "present")) >= minreps))

    reducedallregionsdatakeepers = reducedallregionsdata[keepers,]

    # Convert back to GRanges object
    finalgranges = GRanges(
      reducedallregionsdatakeepers$chr,
      IRanges(
        reducedallregionsdatakeepers$start,
        reducedallregionsdatakeepers$stop
      )
    )
    mcols(finalgranges)[1] = sampletypes[mytype]
    colnames(mcols(finalgranges)) = "sampletype"

    # Construct output
    conspeaks$mytype = finalgranges
    names(conspeaks)[mytype] = sampletypes[mytype]

    # Get some stats for the peaks (before/after merging)
    totconspeaks = NROW(finalgranges)
    names(totconspeaks) = sampletypes[mytype]
    totreppeaks = c()
    for (numreps in 1:length(mytypepeaks)) {
      totreppeaks = c(totreppeaks, NROW(mytypepeaks[[numreps]]))
      names(totreppeaks)[numreps] = names(mytypepeaks)[numreps]
    }
    conspeaks_stats[[mytype]] = c(totconspeaks, totreppeaks)
    names(conspeaks_stats)[[mytype]] = sampletypes[mytype]
  } # end looping through sampletypes

  maxreps = max(c(length(conspeaks_stats[[1]]) - 1, length(conspeaks_stats[[2]]) -
                    1))
  samp1 = conspeaks_stats[[1]]
  samp2 = conspeaks_stats[[2]]
  if (length(samp1) - 1 != maxreps) {
    samp1 = c(samp1, rep(0, maxreps - length(samp1) + 1))
  }
  if (length(samp2) - 1 != maxreps) {
    samp2 = c(samp2, rep(0, maxreps - length(samp2) + 1))
  }

  dfstats = as.data.frame(cbind(c(
    "ConsensusPeaks", paste0("rep", 1:maxreps)
  ), samp1, samp2))
  colnames(dfstats) = c("PeakType", names(conspeaks_stats))

  return(list(consPeaks = conspeaks, consPeaksStats = data.frame(dfstats)))
}
