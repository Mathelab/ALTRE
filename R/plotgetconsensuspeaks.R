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
  dfstats=samplepeaks[[2]]
  dfstats$Replicate=rownames(dfstats)
  mydf=na.omit(reshape2::melt(dfstats))

  p <- ggplot(data=mydf,aes(x=mydf$variable, y=mydf$value,fill=mydf$Replicate)) +
	geom_bar(stat="identity",position=position_dodge()) +
	geom_text(aes(label = mydf$value,x=mydf$variable, y=mydf$value,ymax=mydf$value),
		position=position_dodge(width=1),size=3,hjust=0.5,vjust=-1.5) +
	scale_colour_manual(values = c("red", "dark grey")) +
	labs(fill="") +
	theme_bw(base_size = 15) +
	theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
	labs(x = "Sample Type", y = "Number of Peaks (Regulatory Regions)")

  return(p)
}


