
#' plot barplot
#' @export
plotCountSummary <- function(sizePeaks) {
  dat <- tidyr::gather(sizePeaks,sample,count,-1)

  p   <-
    ggplot(data = dat, aes(x = PeaksFile, y = count, fill = sample)) +
    geom_bar(colour = "black",
             stat = "identity",
             position = position_dodge())
  return(p)
}












