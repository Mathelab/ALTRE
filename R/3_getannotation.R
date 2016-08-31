#' Combine and annotate peaks from different sample types. Optionally merge nearby regions.
#'
#' This function accomplishes three tasks:
#' (1) Combines peaks from different sample types into one master list,
#'  and annotates each peak with it's sample type-specificity (which cell or
#'   tissue types the peak can be found in)
#' (2) Categorizes peaks as either TSS-proximal or TSS-distal based on the distance
#' from a known transcription start site (TSS).  By defaul, distance is set to
#' 1500 bp but can be changed with the distancefromTSS argument).
#' (3) Optionally, regulatory regions that are within a certain distance of
#'  each other can be merged to form a larger regulatory region.
#'
#'
#' @param conspeaks list of GRanges objects for each sample type
#'  (output by getConsensusPeaks() function)
#' @param TSS file of transcription start sites
#' @param distancefromTSS in bp; peaks within distFromTSS of an annotated
#' Transcription Start Site (TSS) will be annotated as TSS-proximal (default=1500bp)
#' @param merge whether or not regions should be merged if they are within a
#' user set distance to each other (default is FALSE)
#' @param regionspecific logical to if TRUE, merging occurs within same type
#'  peaks (e.g. merge TSS-proximal, then merge TSS-distal)
#' @param mergedist merge TSS-proximal and TSS-distal if they are < mergedist apart
#'  (set when regionspecific is FALSE)
#' @param distancefromTSSprox merge TSS-proximal if they are < distancefromTSSprox apart
#'  (set when regionspecific is TRUE)
#' @param distancefromTSSdist merge TSS-distal peaks if they are < distancefromTSSdist apart
#'  (set when regionspecific is TRUE)
#'
#' @return List containing two items:
#' (1) GRanges object:  all sample types combined, regions annotated as
#'  type-specific and TSS-distal/TSS-proximal specific.
#' (2) Matrix: number and size of TSS-distal and TSS-proximal before and
#'  after merging nearby regulatory regions
#'
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
#'}
#' @export

combineAnnotatePeaks <- function(conspeaks,
                                 TSS,
                                 merge = FALSE,
                                 distancefromTSSdist = 0,
                                 distancefromTSSprox = 0,
                                 mergedist = 0,
                                 regionspecific = NA,
                                 distancefromTSS = 1500) {
  if (class(conspeaks[[1]])[1] != "GRangesList") {
    stop("The input conspeaks is not in the correct format!
         (try to rerun getConsensusPeaks())")
  }

  peaklist <- conspeaks[[1]]

  ################################# Aggregate and combine consensus peaks from
  ################################# all sample types
  allregregions <- c(peaklist[[1]])
  for (i in 2:length(peaklist)) {
    allregregions <- c(allregregions, peaklist[[i]])
  }
  reducedallregregions <- reduce(allregregions)

  TSSgranges <- tssannotgrange(reducedallregregions,
                               TSS, distancefromTSS)

  ################################# Annotate combined peaks (TSSgranges) by the
  ################################# cell type they came from

  # If merge is set to false, don't merge
  if (merge == FALSE)
  {
    if (distancefromTSSdist > 0 || distancefromTSSprox >
        0) {
      mergestats <- as.data.frame("No merging because merge is set to FALSE")
    } else {
      mergestats <- as.data.frame("No merging because distancefromTSSdist
                                  and distancefromTSSprox were set to zero")
    }
    namesvector <- c()
    for (i in 1:length(peaklist)) {
      newgranges <- peaklist[[i]]
      name <- as.character(unique(mcols(newgranges)[1])[1, 1])
      namesvector <- c(namesvector, name)
    }

    for (i in 1:length(peaklist)) {
      typespecific <- findOverlaps(TSSgranges,
                                   peaklist[[i]])
      newdataframe <- data.frame(matrix(nrow = length(TSSgranges)))
      newdataframe[queryHits(typespecific), 1] <- namesvector[i]
      values(TSSgranges) <- cbind(values(TSSgranges), newdataframe)
      colnames(mcols(TSSgranges))[i + 1] <- c(namesvector[i])
    }

    # this will annotate the regions with
    # type-specificity
    notmerging <- grangestodataframe(TSSgranges)
    listtoreturn <- list(consPeaksAnnotated =
                           GRanges(notmerging,
                                   meta = notmerging[, 4:ncol(notmerging)]),
                         mergestats = mergestats)
  }  # end if no merging
  else {
    # Do the merging
    if (is.na(regionspecific)) {
      stop("If merging, then the regionspecific paramater must
           be set to TRUE or FALSE")
    }
    # if merging TSS-distal and TSS-proximal regions
    # seperately, then run the merging function on
    # them seperately and then combine them
    # (WITHOUT reducing or you will lose the
    # annotation)
    if (regionspecific == TRUE) {
      if (is.na(distancefromTSSprox) || is.na(distancefromTSSdist)) {
        stop("If regionspecific is true,
             then distancefromTSSprox and distancefromTSSdist must be set")
      }
      dataframeformerge <- grangestodataframe(TSSgranges)
      TSSdistalbeforemergedata <- dataframeformerge[dataframeformerge$region ==
                                                     "TSS-distal", ]
      TSSproxbeforemergedata <- dataframeformerge[dataframeformerge$region ==
                                                     "TSS-proximal", ]

      # Merge TSS-distal and TSS-proximal independently
      # if they're within user defined distances
      TSSdistafter <- mergeclosepeaks(peaklist,
                                       TSSdistalbeforemergedata,
                                       mergedist = distancefromTSSdist,
                                       TSS, distancefromTSS)
      TSSproxafter <- mergeclosepeaks(peaklist,
                                       TSSproxbeforemergedata,
                                       mergedist = distancefromTSSprox,
                                       TSS, distancefromTSS)

      bothafter <- sort(GenomeInfoDb::sortSeqlevels(c(TSSdistafter,
                                        TSSproxafter)))
    }

    # if merging TSS-distal and TSS-proximal regions at
    # the same time, then you just need to run the
    # function once
    if (regionspecific == FALSE) {
      if (is.na(mergedist)) {
        stop("If regionspecific is FALSE, then mergedist must be set")
      }
      dataframeformerge <- grangestodataframe(TSSgranges)
      # create grange from dataframe

      bothafter <- mergeclosepeaks(peaklist,
                                   dataframeformerge,
                                   mergedist,
                                   TSS, distancefromTSS)
    }

    # Create matrix for comparison:
    resultuserinput <- grangestodataframe(bothafter)

    tableofinfo <- matrix(nrow = 4, ncol = 2)
    rownames(tableofinfo) <- c("TSS-distal_before_merging",
                               "TSS-distal_after_merging",
                               "TSS-proximal_before_merging",
                               "TSS-proximal_after_merging")
    colnames(tableofinfo) <- c("TotalNumber", "MeanLength")

    # Callstatscombineannotate internal function
    # to create table:
    tableofinfo <- statscombineannotate(tableofinfo, dataframeformerge, 1, 3)
    tableofinfo <- statscombineannotate(tableofinfo, resultuserinput, 2, 4)

    ### quick fix
    tableofinfo <- as.data.frame(tableofinfo)
    tableofinfo$Condition <- rownames(tableofinfo)
    tableofinfo <- tableofinfo[ , c(3,1,2)]
    rownames(tableofinfo) <- NULL
    ##
    #tableofinfo <- plyr::name_rows(as.data.frame(tableofinfo))[ , c(3,1,2)]
    #colnames(tableofinfo)[1] <- "Condition"
    ###

    listtoreturn <- list(consPeaksAnnotated =
                           GRanges(resultuserinput,
                                   meta = resultuserinput[ , 4:ncol(
                                     resultuserinput)]),
                         mergestats = tableofinfo)
  }
  return(listtoreturn)
}


