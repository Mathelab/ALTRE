#' Given the output from getcounts, plot a density plot of log2 RPKM values of regulation regions
#'
#' @param countsconspeaks output generated from getcounts
#'
#' @return a ggplot
#'
#' @examples
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot=getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot,merge=TRUE,
#'	regionspecific=TRUE,mergedistenh=1500,mergedistprom=1000 )
#' counts_consPeaks <- getcounts(annotpeaks=consPeaksAnnotated, csvfile=csvfile,
#'	reference="SAEC", chrom="chr21")
#' plotgetcounts(counts_consPeaks)
#'
#' @export

plotgetcounts <- function(countsconspeaks) {
    mydf=countsconspeaks$regioncountsforplot
    varstack = suppressMessages(melt(mydf))
    varstack$concat = paste(varstack$region, varstack$variable, sep = ": ")
    varstack$concat = sub("librarysize.*", "", varstack$concat)
    densityplot = ggplot(varstack, aes(x = varstack$value)) +
      geom_density(aes(
        group = varstack$concat,
        color = varstack$concat,
        fill = varstack$concat
      ),
      alpha = 0.3) +
      theme_bw(base_size = 15, base_family = "") +
      theme(legend.title = element_blank()) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(
        #panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(x = "log2 read counts \n(normalized by library and region sizes)")

      return(densityplot)
}
