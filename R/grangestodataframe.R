#' Internal function: converts granges to dataframe
#'
#' Takes a Grange input and converts it into a dataframe
#' Not for use by package user
#'
#' @param grange Grange
#'
#' @return dataframe
#'
#' Dataframe1=grangestodataframe(Grange1)

grangestodataframe <- function(grange) {
  chr = seqnames(grange)
  dataframe = data.frame(chr)
  dataframe = data.frame(lapply(dataframe, as.character), stringsAsFactors =
                           FALSE)
  dataframe$start = start(grange)
  dataframe$stop = end(grange)
  #create the beginning of the dataframe (chr, start, stop, etc...)
  if (length(mcols(grange)) != 0) {
    for (i in 1:length(mcols(grange))) {
      dataframe[, 3 + i] = lapply(mcols(grange)[i], as.character, stringsAsFactors =
                                    FALSE)
    }
  }
  #if the additional columns of the grange are not 0, do this....
  colnames(dataframe) = c("chr", "start", "stop", colnames(mcols(grange)))
  #name columns off grange
  return(dataframe)
}
