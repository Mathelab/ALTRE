#' Internal function: merges peaks within user input distance of each other
#'
#' Takes a data frame of GRanges and distance as input
#' Not for use by package user
#'
#' @param grange data frame of GRanges
#' @param mergedist numeric
#' @param peaklist peaklist
#' @param distancefromTSS distancefromTSS
#' @param TSS TSS
#'
#' @return GRanges object
#'

mergeclosepeaks <-
  function(peaklist,
           grange,
           mergedist,
           TSS,
           distancefromTSS) {
    chr = region = c()
    halfdistance = mergedist / 2

    #how to deal with odd distances entered by users
    if (halfdistance %% 1 != 0) {
      halfdistup = halfdistance + 0.5
      halfdistdown = halfdistance - 0.5
    }
    if (halfdistance %% 1 == 0) {
      halfdistup = halfdistance
      halfdistdown = halfdistance
    }

    #determine the size of each of the regions
    grange$sizeofregionbeforemerge = grange$stop - grange$start

    #add (approximately +/- 0.5) the same amount to both the beginning and end of each region
    #these will be subtracted later so it doesn't really matter what they are INDIVIDUALLY
    #It only matters what number they are TOGETHER ("merged distance")
    grange$start = grange$start - halfdistup
    grange$stop = grange$stop + halfdistdown

    #reduce the merged regions
    beforemerge <-
      GRanges(grange$chr, IRanges(grange$start, grange$stop), meta = region)
    aftermerge = reduce(beforemerge)

    #annotate with transcription start site (again)
    aftermerge = tssannotgrange(aftermerge, TSS, distancefromTSS)

    #get the names of the input files and place them in a vector
    namesvector = c()
    for (i in 1:length(peaklist)) {
      newgranges = peaklist[[i]]
      name = as.character(unique(mcols(newgranges)[1])[1, 1])
      namesvector = c(namesvector, name)
    }

    #this will annotate the regions with type-specificity
    for (i in 1:length(peaklist)) {
      typespecific = findOverlaps(aftermerge, peaklist[[i]])
      newdataframe = data.frame(matrix(nrow = length(aftermerge)))
      newdataframe[queryHits(typespecific), 1] = namesvector[i]
      values(aftermerge) = cbind(values(aftermerge), newdataframe)
      colnames(mcols(aftermerge))[i + 1] = c(namesvector[i])
    }

    aftermergedata = grangestodataframe(aftermerge)
    aftermergedata$start = aftermergedata$start + halfdistup
    aftermergedata$stop = aftermergedata$stop - halfdistdown
    #remove the extra length added to each region for merging purposes
    lengthdata = ncol(aftermergedata)

    #convert back to grange
    aftermerge <-
      GRanges(aftermergedata$chr,
              IRanges(aftermergedata$start, aftermergedata$stop),
              meta = aftermergedata[, c(4:lengthdata)])
    colnames(mcols(aftermerge)) = c("region", namesvector)
    return(aftermerge)
  } # end merging function
