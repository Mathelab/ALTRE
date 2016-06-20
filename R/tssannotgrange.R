#' Internal function: annotates peaks as promoters or enhancers
#'
#' Takes a Grange input
#' Not for use by package user
#'
#' @param grange GRanges object
#' @param TSS GRanges of Transcription Start Site Locations
#' @param distancefromTSS distancefromTSS
#' @return GRanges object
#'

tssannotgrange <- function(grange, TSS, distancefromTSS) {
  #read in TSS, make column header, change to Granges
  #Find distance to transcription start site
  distancetoTSS <- distanceToNearest(grange, TSS)
  #make dataframe from grange
  newdataframe = grangestodataframe(grange)
  newdataframe$distance = mcols(distancetoTSS)$distance
  newdataframe = within(newdataframe, {
    region = ifelse(distance <= distancefromTSS, "promoter", "enhancer")
  })
  #annotate anything <=1500 bp away as promoter, otherwise enhancer
  chr = c()
  annotatedgrange <-
    with(newdataframe, GRanges(chr, IRanges(start, stop), meta = newdataframe[, 5]))
  #create a grange
  colnames(mcols(annotatedgrange)) = c("region")
  return(annotatedgrange)
}
