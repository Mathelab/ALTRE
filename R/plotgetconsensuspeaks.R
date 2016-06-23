#' Given the output from getConsensusPeaks, generate a barplot of count statistics
#'
#' @param samplepeaks output generated from getConsensusPeaks
#'
#' @return a ggplot
#'
#' @examples
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#'
#' @export

plotConsensusPeaks <- function(samplepeaks) {
  dfstats = samplepeaks$consPeaksStats
  #  dfstats$Replicate=dfstats$PeakType
  #  dfstats=dfstats[,-which(colnames(dfstats)=="PeakType")]
  #  mydf=na.omit(reshape2::melt(dfstats))
  mydf = tidyr::gather(dfstats, 'CellType', 'count', 2:3)

  p <-
    ggplot(data = mydf,
           aes_string(
             x = 'CellType',
             y = 'count',
             fill = 'PeakType') ) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(
      aes_string(
        label = 'count',
        x = 'CellType',
        y = 'count',
        ymax = 'count'
      ),
      position = position_dodge(width = 1),
      size = 3,
      hjust = 0.5,
      vjust = -1.5
    ) +
    scale_colour_manual(values = c("red", "dark grey")) +
    labs(fill = "") +
    theme_bw(base_size = 15) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "Sample Type", y = "Number of Peaks (Regulatory Regions)")

  return(p)
}
