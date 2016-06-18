
#' plot barplot
#' @param sizePeaks consPeaksStats
#' @export
plotCountSummary <- function(sizePeaks) {
  dat <- tidyr::gather(sizePeaks,'sample','count',-1)
  p   <-
    ggplot(data = dat,
           aes_string(x = 'PeaksFile', y = 'count', fill='sample')) +
    geom_bar(colour = "black",
             stat = "identity",
             position = position_dodge())
  return(p)
}





