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
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot)
#' counts_consPeaks <-getcounts(annotpeaks=consPeaksAnnotated, csvfile=csvfile,
#'	reference="SAEC", chrom="chr21")
#' altre_peaks=countanalysis(counts=counts_consPeaks, pval=0.01, lfcvalue=1)
#' plotCountAnalysis(altre_peaks)
#'
#' @export

plotCountAnalysis <- function(altrepeaks) {
  toplot= altrepeaks$dftoplot$toplot
  pval=altrepeaks$dftoplot$pval
  lfcvalue=altrepeaks$dftoplot$lfcvalue
  plot = ggplot(toplot, aes(toplot$log2FoldChange, -log2(toplot$padj))) +
    geom_point(aes(col = factor(toplot$col))) +
    scale_colour_manual(values = c("red", "dark grey")) +
    theme_bw(base_size = 15) +
    theme(legend.title = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "log2FC", y = "-log2(pvalue)") +
    geom_hline(aes(yintercept = -log2(pval)), linetype = "dashed") +
    geom_vline(aes(xintercept = (-lfcvalue)), linetype = "dashed") +
    geom_vline(aes(xintercept = (lfcvalue)), linetype = "dashed")

    return(plot)
}
