

#' Given the output from getConsensusPeaks, generate a barplot of count statistics
#'
#' @param samplepeaks output generated from getConsensusPeaks
#'
#' @return a ggplot
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' }
#' @export

plotConsensusPeaks <- function(samplepeaks) {
  dfstats = samplepeaks$consPeaksStats
  #  dfstats$Replicate=dfstats$PeakType
  #  dfstats=dfstats[,-which(colnames(dfstats)=="PeakType")]
  #  mydf=na.omit(reshape2::melt(dfstats))
  mydf = tidyr::gather(dfstats, 'CellType', 'count', 2:3)

  p <-
    ggplot(data = mydf,
           aes_string(x = 'CellType',
                      y = 'count',
                      fill = 'PeakType')) +
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

####################################################################

#' Given the output from combineAnnotatePeaks,
#' plot a barplot showing number of peaks before/after merging
#' (only works if peaks were merged)
#'
#' @param conspeaks output generated from combineAnnotatePeaks
#'
#' @return a ggplot
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot=getTSS()
#'  consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks,
#'  TSS=TSSannot,merge=TRUE, regionspecific=TRUE,mergedistenh=1500,
#'   mergedistprom=1000)
#' plotCombineAnnotatePeaks(consPeaksAnnotated)
#' }
#' @export

plotCombineAnnotatePeaks <- function(conspeaks) {
  mydf = conspeaks$mergestats
  if (nrow(mydf) == 1) {
    stop(
      "No plot to show since merging was not performed when calling combineAnnotatePeaks function"
    )
  }
  else {
    graphics::par(
      mfrow = c(1, 2),
      oma = c(4, 1, 1, 1),
      mar = c(1, 3, 3, 3)
    )
    toplot = data.frame(Enhancers = mydf$total_number[1:2],
                        Promoters = mydf$total_number[3:4])
    rownames(toplot) = c("NoMerge", "Merge")

    b = graphics::barplot(
      as.matrix(toplot),
      ylab = "Number of Regulatory Regions",
      beside = TRUE,
      main = "Number of Regulatory Regions\nBefore/After merging",
      ylim = c(0, max(toplot) + max(toplot) * 0.05),
      col = c("blue4", "coral2")
    )
    graphics::text(
      x = b,
      y = as.matrix(toplot),
      label = as.matrix(toplot),
      pos = 3,
      cex = 0.8
    )

    toplot = data.frame(Enhancers = mydf$mean_length[1:2],
                        Promoters = mydf$mean_length[3:4])
    rownames(toplot) = c("NoMerge", "Merge")

    b = graphics::barplot(
      as.matrix(toplot),
      ylab = "Mean Length of Regulatory Regions",
      beside = TRUE,
      main = "Mean length of Regulatory Regions\nBefore/After merging",
      ylim = c(0, max(toplot) + max(toplot) * 0.05),
      col = c("blue4", "coral2")
    )
    graphics::text(
      x = b,
      y = as.matrix(toplot),
      label = as.matrix(toplot),
      pos = 3,
      cex = 0.8
    )

    graphics::par(
      fig = c(0, 1, 0, 1),
      oma = c(0, 0, 0, 0),
      mar = c(0, 0, 0, 0),
      new = TRUE
    )
    graphics::plot(
      0,
      0,
      type = "n",
      bty = "n",
      xaxt = "n",
      yaxt = "n"
    )
    graphics::legend(
      "bottom",
      c("No Merging", "After Merging"),
      fill = c("blue4", "coral2"),
      horiz = TRUE,
      inset = c(0, 0),
      xpd = TRUE,
      bty = "n"
    )
  }
}



#' Given the output from getcounts, plot a density plot of log2 RPKM values of regulation regions
#'
#' @param altrepeaks output generated from getcounts
#'
#' @return a ggplot
#'
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot=getTSS()
#'  consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot,
#'	regionspecific=TRUE,mergedistenh=1500,mergedistprom=1000 )
#' counts_consPeaks <-getcounts(annotpeaks=consPeaksAnnotated, csvfile=csvfile,
#'	reference="SAEC", chrom="chr21")
#' altre_peaks=countanalysis(counts=counts_consPeaks, pval=0.01, lfcvalue=1)
#' plotCountAnalysis(altre_peaks)
#' }
#' @export

plotCountAnalysis <- function(altrepeaks) {
  toplot = altrepeaks$dftoplot$toplot
  pval = altrepeaks$dftoplot$pval
  lfcvalue = altrepeaks$dftoplot$lfcvalue
  plot = ggplot(toplot, aes(toplot$log2FoldChange,-log2(toplot$padj))) +
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



###############################################################################


#' Given the output from getcounts, plot a density plot of log2 RPKM values of regulation regions
#'
#' @param countsconspeaks output generated from getcounts
#'
#' @return a ggplot
#'
#' @examples
#' \dontrun{
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
#' }
#' @export

plotgetcounts <- function(countsconspeaks) {
  mydf = countsconspeaks$regioncountsforplot
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
    theme(#panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
    labs(x = "log2 read counts \n(normalized by library and region sizes)")

  return(densityplot)
}


#' Multiple plot function
#'
#' Plots multiple ggplot objects in one window.
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' This function was not written by the authours of this package; it was found here:
#' http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#'
#' @param ... list of plots
#' @param cols number of columns in layout
#' @param layout matrix specifying layout, if present cols is ignored
#' @param plotlist plotlist
#' @param file file

multiplot <-
  function(...,
           plotlist = NULL,
           file,
           cols = 1,
           layout = NULL) {
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                       ncol = cols,
                       nrow = ceiling(numPlots / cols))
    }

    if (numPlots == 1) {
      print(plots[[1]])

    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

        print(plots[[i]],
              vp = viewport(
                layout.pos.row = matchidx$row,
                layout.pos.col = matchidx$col
              ))
      }
    }
  }
