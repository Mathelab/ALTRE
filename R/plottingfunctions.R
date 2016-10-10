#' Given the output from getConsensusPeaks, generate a barplot
#' of countstatistics
#'
#' @param samplepeaks output generated from getConsensusPeaks
#' @param palette RColorBrewer palette to change graph colors
#' @param xlabel label for x-axis (default, types of peaks - e.g. ConsensusPeaks, rep1, rep2, etc.)
#' @param ylabel label for y-axis (default, "Peak Counts")
#' @param xlabelsize size of xlabel (default, 15px)
#' @param ylabelsize size of ylabel (default, 15px)
#' @param maintitle main title (default, "Peak Counts by Cell Type")
#' @param maintitlesize main title size (default, 20px)
#' @param subtitle subtitle (default, "For bioreplicates and their merged consensus track")
#' @param subtitlesize subitle size (default 15px)
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' csvfile <- loadCSVFile("DNAseEncodeExample.csv")
#' samplePeaks <- loadBedFiles(csvfile)
#' consensusPeaks <- getConsensusPeaks(samplepeaks = samplePeaks,
#' minreps = 2)
#' plotConsensusPeaks(samplepeaks = consensusPeaks)}
#' @export
#'
plotConsensusPeaks <- function(samplepeaks,
                               palette = NULL,
                               xlabel = NULL,
                               ylabel = "Peak Counts",
                               xlabelsize = '15px',
                               ylabelsize ='15px',
                               maintitle = "Peak Counts by Cell Type",
                               subtitle = "For bioreplicates and their merged consensus track",
                               maintitlesize = "20px",
                               subtitlesize = "15px") {

  if ( !is.null(palette) ) {
    cols <- RColorBrewer::brewer.pal(3, palette)
  } else {
    cols <- RColorBrewer::brewer.pal(3, "Set1")
  }
  CellType <- NULL
  #without this R CMD check throws no visible binding for global variable error
  consPeaksStats <- samplepeaks$consPeaksStats
  row.names(consPeaksStats) <- NULL
  consPeaksStats[, 2] <-
    as.numeric(as.character(consPeaksStats[[2]]))
  consPeaksStats[, 3] <-
    as.numeric(as.character(consPeaksStats[[3]]))
  statsFormated <-
    tidyr::gather(consPeaksStats, "CellType", "Count", 2:3)

  plottingData <- statsFormated %>%
    split(levels(as.factor(statsFormated$CellType)))

  col1 <- cols[1]
  col2 <- cols[2]

  if ( is.null(xlabel) ) {
  xLabel <- as.character(samplepeaks$consPeaksStats$PeakType)
  } else {
    xLabel <- xlabel
  }

  p <- highchart(width = 520, height = 650) %>%
    hc_title(text = maintitle,
             style = list(color = '#2E1717',
                          fontSize = maintitlesize,
                          fontWeight = 'bold')) %>%
    hc_subtitle(text = subtitle,
                style = list(fontSize = subtitlesize)) %>%
    hc_add_series(
      color = col1,
      data = plottingData[[1]]$Count,
      name = names(plottingData[1]),
      type = "column",
      dataLabels = list(
        enabled = TRUE,
        rotation = 270,
        color = '#FFFFFF',
        y = 40
      )
    ) %>%
    hc_add_series(
      color = col2,
      data = plottingData[[2]]$Count,
      name = names(plottingData[2]),
      dataLabels = list(
        enabled = TRUE,
        rotation = 270,
        color = '#FFFFFF',
        y = 40
      ),
      type = "column"
    ) %>%
    hc_yAxis(title = list(text = ylabel,
                          style = list(fontSize = ylabelsize)),
             labels = list(format = "{value}")) %>%
    hc_xAxis(categories = xLabel,
	            labels = list(style = list(fontSize = xlabelsize))) %>%
    hc_legend(
      enabled = TRUE,
      layout = "horizontal",
      align = "center",
      horizontalAlign = "middle",
      verticalAlign = "top",
      floating = FALSE,
      maxHeight = 100,
      y = 50
    ) %>%
    hc_tooltip(
      headerFormat = "<b>{series.name}_{point.key}</b><br>",
      pointFormat = "{point.y}",
      valueSuffix = ' peaks'
    ) %>%
    hc_exporting(enabled = TRUE)
  return(p)
}

####################################################################

#' Given the output from combineAnnotatePeaks,
#' plot a barplot showing number of peaks before/after merging or length of
#' peaks before/after merging
#' (only works if peaks were merged)
#'
#' @param conspeaks output generated from combineAnnotatePeaks
#' @param palette RColorBrewer palette to change graph colors
#' @param viewer whether the plot should be displayed in the RStudio viewer or
#'        in Shiny/Knittr
#' @param xlabels labels for x-axis (default, Before merging/After merging)
#' @param leftylabel label for y-axis for number of REs plot (default, "Number of REs")
#' @param rightylabel label for y-axis for mean length plot (default, "Mean length of REs")
#' @param xlabelsize size of xlabel (default, 15px)
#' @param leftylabelsize size of leftylabel (default, 15px)
#' @param rightylabelsize size of rightylabel (default, 15px)
#' @param leftmaintitle main title for number of REs plot (default, "Number of REs")
#' @param rightmaintitle main title for mean length of REs plot (default, "Mean length of REs")
#' @param leftmaintitlesize main title size for number of REs plot (default, 20px)
#' @param rightmaintitlesize main title size for mean length of REs plot (default 15px)
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' csvfile <- loadCSVFile("DNAseEncodeExample.csv")
#' samplePeaks <- loadBedFiles(csvfile)
#' consensusPeaks <- getConsensusPeaks(samplepeaks = samplePeaks,
#' minreps = 2)
#' TSSannot <- getTSS()
#' consensusPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consensusPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           distancefromTSSdist = 1500,
#'                                           distancefromTSSprox = 1000)
#' plotCombineAnnotatePeaks(consensusPeaksAnnotated)
#' }
#' @export
#'
plotCombineAnnotatePeaks <- function(conspeaks,
                                     viewer = TRUE,
                                     palette = NULL,
                                     rightmaintitle = "Mean Length of REs",
                                     leftmaintitle = "Number of REs",
                                     xlabels = c("Before merging","After merging"),
                                     leftylabel = "Number of REs",
                                     rightylabel = "Mean Length of REs",
                                     xlabelsize = "15px",
                                     leftylabelsize = "15px",
                                     rightylabelsize = "15px",
                                     leftmaintitlesize = "20px",
                                     rightmaintitlesize = "20px") {
  if ( !is.null(palette) ) {
    cols <- RColorBrewer::brewer.pal(4, palette)
  } else {
    cols <- RColorBrewer::brewer.pal(4, "Set1")
  }

  CellType <- NULL
  #R CMD check throws no visible binding for global variable error

  mergeStats <- conspeaks$mergestats
  row.names(mergeStats) <- NULL
  mergeStats[, 2] <-  as.numeric(as.character(mergeStats[[2]]))
  mergeStats[, 3] <-  as.numeric(as.character(mergeStats[[3]]))
  mergeStatsFormatted <-
    tidyr::gather(mergeStats, "CellType", "Count", 2:3)


  if (nrow(mergeStatsFormatted) == 1) {
    stop(
      "No plot to show since merging was not performed
      when calling combineAnnotatePeaks function"
    )
  }


  feature <- "TotalNumber"
  if (feature == "TotalNumber") {
    mergeStatsTotal <- dplyr::filter(mergeStatsFormatted,
                                     CellType == "TotalNumber")
    thecondition <-
      matrix(unlist(strsplit(mergeStatsTotal$Condition, "_")),
             nrow = 3, ncol = 4)[2, ]
    mergeStatsBefore <- dplyr::filter(mergeStatsTotal,
                                      thecondition == "before")
    mergeStatsAfter <- dplyr::filter(mergeStatsTotal,
                                     thecondition == "after")

    p1 <- highchart(height = 400) %>%
      hc_title(text = leftmaintitle,
               style = list(color = '#2E1717',
                            fontSize = leftmaintitlesize,
                            fontWeight = 'bold')) %>%
      hc_add_series(
        data = mergeStatsBefore$Count,
        name = c("TSS-distal"),
        type = "column",
        dataLabels = list(
          enabled = TRUE,
          rotation = 270,
          color = '#FFFFFF',
          y = 40
        )
      ) %>%
      hc_add_series(
        data = mergeStatsAfter$Count,
        name = c("TSS-proximal"),
        type = "column",
        dataLabels = list(
          enabled = TRUE,
          rotation = 270,
          color = '#FFFFFF',
          y = 40
        )
      ) %>%
      hc_yAxis(title = list(text = leftylabel,
                            style = list(fontSize = leftylabelsize)),
               labels = list(format = "{value}")) %>%
      hc_xAxis(categories = xlabels,
               labels = list(style = list(fontSize = xlabelsize))) %>%
      hc_legend(
        enabled = TRUE,
        layout = "horizontal",
        align = "center",
        verticalAlign = "bottom",
        floating = FALSE,
        maxHeight = 100,
        x = 0,
        y = 17
      ) %>%
      hc_tooltip(
        headerFormat = "<b>{series.name}_{point.key}</b><br>",
        pointFormat = "{point.y}",
        valueSuffix = ' peaks'
      ) %>%
      hc_colors(cols) %>%
      hc_exporting(enabled = TRUE)

  }

  feature <- "MeanLength"
  if (feature == "MeanLength") {
    mergeStatsMean <- dplyr::filter(mergeStatsFormatted,
                                    CellType == "MeanLength")
    thecondition <-
      matrix(unlist(strsplit(mergeStatsMean$Condition, "_")),
             nrow = 3,
             ncol = 4)[2,]
    mergeStatsBefore <- dplyr::filter(mergeStatsMean,
                                      thecondition == "before")
    mergeStatsAfter <- dplyr::filter(mergeStatsMean,
                                     thecondition == "after")


    p2 <- highchart(height = 400) %>%
      hc_title(text = rightmaintitle,
               style = list(color = '#2E1717',
                            fontSize = rightmaintitlesize,
                            fontWeight = 'bold')) %>%
      hc_add_series(
        data = mergeStatsBefore$Count,
        name = c("TSS-distal"),
        type = "column",
        dataLabels = list(
          enabled = TRUE,
          rotation = 270,
          color = '#FFFFFF',
          y = 40
        )
      ) %>%
      hc_add_series(
        data = mergeStatsAfter$Count,
        name = c("TSS-proximal"),
        type = "column",
        dataLabels = list(
          enabled = TRUE,
          rotation = 270,
          color = '#FFFFFF',
          y = 40
        )
      ) %>%
      hc_yAxis(
        title = list(text = rightylabel,
                     style = list(fontSize = rightylabelsize)),
        labels = list(format = "{value}")
      ) %>%
      hc_xAxis(categories = xlabels,
               labels = list(style = list(fontSize = xlabelsize))) %>%
      hc_legend(
        enabled = TRUE,
        layout = "horizontal",
        align = "center",
        verticalAlign = "bottom",
        floating = FALSE,
        maxHeight = 100,
        x = 0,
        y = 17
      ) %>%
      hc_tooltip(
        headerFormat = "<b>{series.name}_{point.key}</b><br>",
        pointFormat = "{point.y}",
        valueSuffix = ' peaks'
      ) %>%
      hc_colors(cols) %>%
      hc_exporting(enabled = TRUE)
  }


  if (viewer == TRUE) {
    p <-
      htmltools::browsable(hw_grid(p1, p2, ncol = 2, rowheight = 550))
  }
  else {
    p <- hw_grid(p1, p2)
  }
  return(p)
  }

###############################################################################


#' Given the output from getCounts, plot a density plot of
#' log2 RPKM values of regulation regions
#'
#' @param countsConsPeaks output generated from getCounts
#' @param palette RColorBrewer palette to change graph colors
#' @param xlabel label for x-axis (default, "log2 normalized read counts")
#' @param ylabel label for y-axis (default, "Density")
#' @param xlabelsize size of xlabel (default, 15px)
#' @param ylabelsize size of ylabel (default, 15px)
#' @param maintitle main title (default, "Density of log2 read counts
#'      (normalized by library and region sizes)")
#' @param maintitlesize main title size (default, 20px)
#'
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' csvfile <- loadCSVFile("DNAseEncodeExample.csv")
#' samplePeaks <- loadBedFiles(csvfile)
#' consensusPeaks <- getConsensusPeaks(samplepeaks = samplePeaks,
#' minreps = 2)
#' TSSannot <- getTSS()
#' consensusPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consensusPeaks,
#' TSS = TSSannot,
#' merge = TRUE,
#' regionspecific = TRUE,
#' distancefromTSSdist = 1500,
#' distancefromTSSprox = 1000)
#' consensusPeaksCounts <- getCounts(annotpeaks = consensusPeaksAnnotated,
#'                                  sampleinfo = csvfile,
#'                                  reference = 'SAEC',
#'                                  chrom = 'chr21')
#' plotGetCounts(consensusPeaksCounts)}
#' @export

plotGetCounts <- function(countsConsPeaks,
            palette = NULL,
            xlabel = "Log2 Normalized Read Counts",
            ylabel = "Density",
            xlabelsize = "15px",
            ylabelsize = "15px",
            maintitle = "Density of log2 read counts\n(normalized by library and region sizes)",
            maintitlesize = "20px"
            ) {
  region <- NULL
  variable <- NULL
  value = Reference_specific = Shared = Experiment_specific = c()
  #set to null forR CMD Check error: Undefined global functions/variables

  if ( !is.null(palette) ) {
    cols <- RColorBrewer::brewer.pal(4, palette)
  } else {
    cols <- RColorBrewer::brewer.pal(4, "Set1")
  }


  mydf <- countsConsPeaks$regioncountsforplot
  #  varstack <- suppressMessages(reshape2::melt(mydf))
  varstack <-
    suppressMessages(tidyr::gather(mydf, variable, value, -region))
  varstack$variable = gsub("_.*", "", varstack$variable)

  mysamps <- unique(varstack$variable)

  samp1dist <- dplyr::filter(varstack,
                             region == "TSS-distal" &
                               variable == mysamps[1])
  samp2dist <- dplyr::filter(varstack,
                             region == "TSS-distal" &
                               variable == mysamps[2])
  samp1prox <- dplyr::filter(varstack,
                             region == "TSS-proximal" &
                               variable == mysamps[1])
  samp2prox <- dplyr::filter(varstack,
                             region == "TSS-proximal" &
                               variable == mysamps[2])


  round <- JS(
      "function() { return '<b>'+'log2 read count' +'</b>:'+
    Highcharts.numberFormat(this.x, 2) + ', <b>'+'density' +'</b>:' +
      Highcharts.numberFormat(this.y, 2); }"
  )

  p <- hchart(
    stats::density(samp1dist$value),
    area = TRUE,
    name = paste(mysamps[1], "TSS-distal")
  ) %>%
    hc_title(
      text = maintitle,
      style = list(color = '#2E1717',
                   fontSize = maintitlesize,
                   fontWeight = 'bold')
    ) %>%
    hc_yAxis(title = list(text = ylabel,
                          style = list(fontSize = ylabelsize))) %>%
    hc_xAxis(title = list(text = xlabel,
                          style = list(fontSize = xlabelsize))) %>%
    hc_add_series_density(
      stats::density(samp1prox$value),
      area = TRUE,
      name = paste(mysamps[1], "TSS-proximal")
    ) %>%
    hc_add_series_density(
      stats::density(samp2dist$value),
      area = TRUE,
      name = paste(mysamps[2], "TSS-distal")
    ) %>%
    hc_add_series_density(
      stats::density(samp2prox$value),
      area = TRUE,
      name = paste(mysamps[2], "TSS-proximal")
    ) %>%
    hc_colors(cols) %>%
    hc_tooltip(formatter = round) %>%
    hc_exporting(enabled = TRUE) %>%
    hc_legend(
      enabled = TRUE,
      layout = "horizontal",
      align = "center",
      horizontalAlign = "middle",
      verticalAlign = "top",
      floating = TRUE,
      maxHeight = 100,
      y = 50)
  return(p)
}


#' Create a volcano plot from the output of categAltrePeaks
#'
#' @param altrepeakscateg output generated from countanalysis() then
#' categAltrePeaks()
#' @param viewer whether the plot should be displayed in the RStudio viewer or
#'        in Shiny/Knittr
#' @param palette RColorBrewer palette to change graph colors
#' @param xlabel label for x-axis (default, "-log10 pvalue")
#' @param ylabel label for y-axis (default, "Density")
#' @param xlabelsize size of xlabel (default, 15px)
#' @param ylabelsize size of ylabel (default, 15px)
#' @param maintitlelefty main title (default, "TSS-distal")
#' @param maintitlerighty main title (default, "TSS-proximal")
#' @param maintitlesize main title size (default, 20px)
#'
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' csvfile <- loadCSVFile("DNAseEncodeExample.csv")
#' samplePeaks <- loadBedFiles(csvfile)
#' consensusPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps = 2)
#' TSSannot <- getTSS()
#' consensusPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consensusPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           distancefromTSSdist = 1500,
#'                                           distancefromTSSprox = 1000)
#' consensusPeaksCounts <- getCounts(annotpeaks = consensusPeaksAnnotated,
#'                               reference = 'SAEC',
#'                               sampleinfo = csvfile,
#'                               chrom = 'chr21')
#' alteredPeaks <- countanalysis(counts=consensusPeaksCounts,
#'                              pval=0.01,
#'                              lfcvalue=1)
#' alteredPeaksCategorized <- categAltrePeaks(alteredPeaks,
#'                              lfctypespecific = 1.5,
#'                              lfcshared = 1.2,
#'                              pvaltypespecific = 0.01,
#'                              pvalshared = 0.05)
#' plotCountAnalysis(altrepeakscateg = alteredPeaksCategorized)
#' }
#' @export

plotCountAnalysis <- function(altrepeakscateg, viewer = TRUE, palette = NULL,
                              maintitlelefty = "TSS-distal",
                              maintitlerighty = "TSS-proximal",
                              ylabel = "-log10 pvalue",
                              xlabel = "log2fold change",
                              xlabelsize = "15px",
                              ylabelsize = "15px",
                              maintitlesize = "20px") {

    if ( !is.null(palette) ) {
        cols <- RColorBrewer::brewer.pal(4, palette)
    } else {cols <- c("#C71585", "#d3d3d3", "#00E5EE", "#000080")}
                        #magenta (experiment-specific)
                        #grey (ambiguous)

                        #dark blue (shared)
                            #blue (ref)

  log2FoldChange <- NULL
  padj <- NULL
  REaltrecateg <- REaltrecategplot <- NULL
  #To prevent R CMD check error

  Referencespecificsamples <- altrepeakscateg[[3]]
  allsamples <- colnames(altrepeakscateg$analysisresults)[12:13]
  Experimentspecificsamples <- allsamples[which(!(allsamples %in% Referencespecificsamples))]

  Referencespecific <- paste0(Referencespecificsamples, "SpecificByIntensity")
  Experimentspecific <- paste0(Experimentspecificsamples, "SpecificByIntensity")

  Referencespecificlabels <- paste0(Referencespecificsamples, "-Specific (by intensity)")
  Experimentspecificlabels <- paste0(Experimentspecificsamples, "-Specific (by intensity)")


  toplot <- altrepeakscateg$analysisresults[ , c("region",
                                                "log2FoldChange",
                                                "padj",
                                                "REaltrecategplot")]
  replacement <- sub(Referencespecific, Referencespecificlabels, toplot$REaltrecategplot)
  replacement <- sub(Experimentspecific, Experimentspecificlabels, replacement)
  toplot$REaltrecategplot <- replacement

  tssdist <- toplot[which(toplot$region == "TSS-distal"), ]
  tssdist$padj <- -log10(tssdist$padj)
  tssprox <- toplot[which(toplot$region == "TSS-proximal"), ]
  tssprox$padj <- -log10(tssprox$padj)

  # remove the NAs
  tssdist <- tssdist[!is.na(tssdist$padj), -1]
  tssprox <- tssprox[!is.na(tssprox$padj), -1]

  # order values
  tssdist <- tssdist[order(tssdist$padj,
                           tssdist$log2FoldChange,
                           decreasing = TRUE), ]
  tssprox <- tssprox[order(tssprox$padj,
                           tssprox$log2FoldChange,
                           decreasing = TRUE), ]

  if ((dim(tssdist)[1] + dim(tssprox)[1]) > 6000) {

    message("Plot doesn't show all data points. Data around the origin has been trimmed")

    # remove shared
    tssdist <- tssdist[!(tssdist$REaltrecategplot == "Shared"), ]
    tssprox <- tssprox[!(tssprox$REaltrecategplot == "Shared"), ]

    #####################
    # trim the data

    cutoff_p <- 5
    cutoff_fold <- 4

    tssdistEdges <- tssdist[tssdist$padj > cutoff_p &
                              (abs(tssdist$log2FoldChange) > cutoff_fold), ]
    tssdistCenter <- tssdist[!(tssdist$padj > cutoff_p &
                                 (abs(tssdist$log2FoldChange) > cutoff_fold)), ]

    tssproxEdges <- tssprox[tssprox$padj > cutoff_p &
                              (abs(tssprox$log2FoldChange) > cutoff_fold), ]
    tssproxCenter <- tssprox[!(tssprox$padj > cutoff_p &
                                 (abs(tssprox$log2FoldChange) > cutoff_fold)), ]

    n1 <- dim(tssdistCenter)[1]
    n2 <- dim(tssproxCenter)[1]

    nsamp <- 1000

    idx1 <- sort(stats::rexp(min(nsamp ,  3*(n1 %/% 4)), 2))
    idx1 <- unique(floor(n1 * (idx1 / max(idx1))))

    idx2 <- sort(stats::rexp(min(nsamp ,  3*(n2 %/% 4)), 2))
    idx2 <- unique(floor(n2 * (idx2 / max(idx2))))

    tssdistCenter <- tssdistCenter[idx1, ]
    tssproxCenter <- tssproxCenter[idx2, ]

    tssdist <- rbind(tssdistEdges, tssdistCenter)
    tssprox <- rbind(tssproxEdges, tssproxCenter)

  } else{

    message("Plot shows all data points. Data around the origin has not been trimmed")
  }



  # upperThresh1 <- max((n1 - 1000),  3 * (n1 %/% 4))
  # upperThresh2 <- max((n2 - 1000),  3 *(n2 %/% 4))
  # lowerThresh1 <-   (n1 %/% 4)
  # lowerThresh2 <-   (n2 %/% 4)
  #
  # idx1 <- c(sample((lowerThresh1 + 1):upperThresh1,
  #                       min(3000, (upperThresh1 - lowerThresh1))),
  #               (upperThresh1 + 1):n1
  #                )
  # idx2 <- c(sample((lowerThresh2 + 1):upperThresh2,
  #                       min(3000, (upperThresh2 - lowerThresh2))),
  #               (upperThresh2 + 1):n2
  #               )


  ###########################

    p1 <- highchart() %>%
      hc_chart(type = "scatter") %>%
      hc_plotOptions(
        scatter = list(marker = list(radius = 2),
                       turboThreshold = 0)
      ) %>%
      hc_title(
        text = maintitlelefty,
        style = list(color = '#2E1717',
                     fontSize = maintitlesize,
                     fontWeight = 'bold')
      ) %>%
      hc_add_series_df(
        data = tssdist,
        x = log2FoldChange,
        y = padj,
        type = "scatter",
        group = REaltrecategplot
      )  %>%
      hc_yAxis(title = list(text = ylabel,
                            style = list(fontSize = ylabelsize))) %>%
      hc_xAxis(title = list(text = xlabel,
                            style = list(fontSize = xlabelsize))) %>%
      hc_tooltip(headerFormat = "",
                 pointFormat  = "<b>log2FC</b> = {point.x}<br> <b>-log10pvalue</b>
                = {point.y}<br>") %>%
      hc_colors(cols) %>%
      hc_exporting(enabled = TRUE)

    p2 <- highchart() %>%
      hc_plotOptions(
        scatter = list(marker = list(radius = 2))
      ) %>%
      hc_chart(type = "scatter") %>%
      hc_title(
        text = maintitlerighty,
        style = list(color = '#2E1717',
                     fontSize = maintitlesize,
                     fontWeight = 'bold')
      ) %>%
      hc_add_series_df(
        data = tssprox,
        x = log2FoldChange,
        y = padj,
        type = "scatter",
        group = REaltrecategplot
      )  %>%
      hc_yAxis(title = list(text = ylabel,
                            style = list(fontSize = ylabelsize))) %>%
      hc_xAxis(title = list(text = xlabel,
                            style = list(fontSize = xlabelsize))) %>%
      hc_tooltip(headerFormat = "",
                 pointFormat  = "<b>log2FC</b> = {point.x}
                 <br> <b>-log10pvalue</b> = {point.y}<br>",
                 valueDecimals = 2) %>%
      hc_colors(cols) %>%
      hc_exporting(enabled = TRUE)

    if (viewer == TRUE) {
      p <-
        htmltools::browsable(hw_grid(p1, p2, ncol = 2))
    }
    else {
      p <- hw_grid(p1, p2, ncol = 2)
    }
    return(p)
  }

###############################################################################
#' Creates a boxplot to see the distribution of read counts in type-specific and
#' shared TSS-proximal and TSS-distal regions.
#'
#' Takes the rlog transformation of the RRKM (Reads Per Kilobase of transcript
#' per Million) of the read counts of type-specific and shared regulatory
#' regions and plots the distribution of those read counts in all sample types
#' analyzed in the workflow.
#'
#' @param analysisresults output generated from countanalysis() then
#' categAltrePeaks()
#' @param counts output generated from getCounts()
#' @param palette RColorBrewer palette to change graph colors
#' @param ylabel label for y-axis (default, "Observations")
#' @param ylabelsize size of ylabel (default, 15px)
#' @param xlabelsize size of xlabel (default, 15px)
#' @param xlabel label for x-axis (default, sample names)

#' @param maintitle main title (default, "Distribution of Normalized Counts")
#' @param maintitlesize main title size (default, 20px)
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' csvfile <- loadCSVFile("DNAseEncodeExample.csv")
#' samplePeaks <- loadBedFiles(csvfile)
#' consensusPeaks <- getConsensusPeaks(samplepeaks = samplePeaks,
#' minreps = 2)
#' TSSannot <- getTSS()
#' consensusPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consensusPeaks,
#' TSS = TSSannot,
#' merge = TRUE,
#' regionspecific = TRUE,
#' distancefromTSSdist = 1500,
#' distancefromTSSprox = 1000)
#' consensusPeaksCounts <- getCounts(annotpeaks = consensusPeaksAnnotated,
#'                                  sampleinfo = csvfile,
#'                                  reference = 'SAEC',
#'                                  chrom = 'chr21')
#' alteredPeaks <- countanalysis(counts = consensusPeaksCounts,
#' pval = 0.01,
#' lfcvalue = 1)
#' alteredPeaksCategorized <- categAltrePeaks(alteredPeaks,
#'                                           lfctypespecific = 1.5,
#'                                           lfcshared = 1.2,
#'                                           pvaltypespecific = 0.01,
#'                                           pvalshared = 0.05)
#' plotDistCountAnalysis(alteredPeaksCategorized, consensusPeaksCounts)
#' }
#' @export
#'
plotDistCountAnalysis <-
  function(analysisresults,
           counts,
           palette = NULL,
           xlabelsize = "13px",
           ylabel = "log2(FPKM)",
           ylabelsize = "13px",
           maintitle = "Distribution of Normalized Counts (peaks types determine by intensity)",
           maintitlesize = "20px",
           xlabel = NULL) {
    altrecateg <- altrecategplot <- REaltrecategplot <- c()
    #Make sure to names things are from the user-entered sample names
    reference <- analysisresults$reference
    allSamples <- colnames(analysisresults$analysisresults)[12:13]
    nonreference <- allSamples[which(!(allSamples %in% reference))]
    Referencespecific <- paste0(reference, "SpecificByIntensity")
    Experimentspecific <- paste0(nonreference, "SpecificByIntensity")

    if (!is.null(palette)) {
      cols <- RColorBrewer::brewer.pal(4, palette)
    } else{
      cols <- c("#C71585", "#d3d3d3", "#000080", "#00E5EE")
    }
    #magenta (experiment-specific) #grey (ambiguous) #blue (shared))
    #blue (reference specific)

    readcounts <- counts$regioncounts
    analysisresults <- analysisresults$analysisresults
    errortest = try(SummarizedExperiment::assay(readcounts), silent = TRUE)
    if (inherits(errortest, 'try-error') == TRUE) {
      stop("The input for the readcounts arguement is
           not a summerized experiment object!")
    }

    if (is.data.frame(analysisresults) == FALSE)
    {
      stop("The input for the analysisresults arguement is not a dataframe!")

    }

    # Check that counts and analysisresults are in the same order
    countsinfo <-
      as.data.frame(SummarizedExperiment::rowRanges(readcounts))
    countcoord <-
      paste0(countsinfo$seqnames, countsinfo$start, countsinfo$end)
    analcoord <- paste0(analysisresults$chr,
                        analysisresults$start,
                        analysisresults$stop)

    if (!all.equal(analcoord, countcoord)) {
      stop("The peaks in the analysisresults and counts are not the same")
    }

    PEcateg <- analysisresults$region
    altrecategplot <- analysisresults$REaltrecategplot

    # Get log2FPM values:
    log2FPM <- log2(DESeq2::fpkm(readcounts, robust = TRUE) + 0.001)

    # Average log2FPM values over replicats:
    sampletypes <- SummarizedExperiment::colData(readcounts)$sample
    meanlog2FPM <- c()

  for (i in unique(sampletypes)) {
    samp <- which(sampletypes == i)
    meanlog2FPM <- cbind(meanlog2FPM,
                        as.numeric(apply(log2FPM[, samp], 1, mean)))
  }
  colnames(meanlog2FPM) <- unique(sampletypes)

  mydf <- data.frame(meanlog2FPM = meanlog2FPM,
                    PEcateg = PEcateg,
                    altrecateg = altrecategplot)
  #TSSdistal <- dplyr::filter(mydf, PEcateg == "TSS-distal")
  distal1 <- dplyr::filter(mydf, altrecateg == Experimentspecific)
  distal2 <- dplyr::filter(mydf, altrecateg == "Ambiguous")
  distal3 <- dplyr::filter(mydf, altrecateg == "Shared")
  distal4 <- dplyr::filter(mydf, altrecateg == Referencespecific)

  #TSSproximal <- dplyr::filter(mydf, PEcateg == "TSS-proximal")
  proximal1 <- dplyr::filter(mydf, altrecateg == Experimentspecific)
  proximal2 <- dplyr::filter(mydf, altrecateg == "Ambiguous")
  proximal3 <- dplyr::filter(mydf, altrecateg == "Shared")
  proximal4 <- dplyr::filter(mydf, altrecateg == Referencespecific)

  mysamps = as.character(unique(sampletypes))
  distal1_5num_samp1 <-
    stats::fivenum(distal1[[paste("meanlog2FPM", mysamps[1], sep = ".")]])
  proximal1_5num_samp1 <-
    stats::fivenum(proximal1[[paste("meanlog2FPM", mysamps[1], sep = ".")]])
  distal1_5num_samp2 <-
    stats::fivenum(distal1[[paste("meanlog2FPM", mysamps[2], sep = ".")]])
  proximal1_5num_samp2 <-
    stats::fivenum(proximal1[[paste("meanlog2FPM", mysamps[2], sep = ".")]])

  distal2_5num_samp1 <-
    stats::fivenum(distal2[[paste("meanlog2FPM", mysamps[1], sep = ".")]])
  proximal2_5num_samp1 <-
    stats::fivenum(proximal2[[paste("meanlog2FPM", mysamps[1], sep = ".")]])
  distal2_5num_samp2 <-
    stats::fivenum(distal2[[paste("meanlog2FPM", mysamps[2], sep = ".")]])
  proximal2_5num_samp2 <-
    stats::fivenum(proximal2[[paste("meanlog2FPM", mysamps[2], sep = ".")]])

  distal3_5num_samp1 <-
    stats::fivenum(distal3[[paste("meanlog2FPM", mysamps[1], sep = ".")]])
  proximal3_5num_samp1 <-
    stats::fivenum(proximal3[[paste("meanlog2FPM", mysamps[1], sep = ".")]])
  distal3_5num_samp2 <-
    stats::fivenum(distal3[[paste("meanlog2FPM", mysamps[2], sep = ".")]])
  proximal3_5num_samp2 <-
    stats::fivenum(proximal3[[paste("meanlog2FPM", mysamps[2], sep = ".")]])

  distal4_5num_samp1 <-
    stats::fivenum(distal4[[paste("meanlog2FPM", mysamps[1], sep = ".")]])
  proximal4_5num_samp1 <-
    stats::fivenum(proximal4[[paste("meanlog2FPM", mysamps[1], sep = ".")]])
  distal4_5num_samp2 <-
    stats::fivenum(distal4[[paste("meanlog2FPM", mysamps[2], sep = ".")]])
  proximal4_5num_samp2 <-
    stats::fivenum(proximal4[[paste("meanlog2FPM", mysamps[2], sep = ".")]])

  Experimentspecific_list <- list(distal1_5num_samp1,
                                  proximal1_5num_samp1,
                                  distal1_5num_samp2,
                                  proximal1_5num_samp2)
  Ambiguous_list <- list(distal2_5num_samp1,
                         proximal2_5num_samp1,
                         distal2_5num_samp2,
                         proximal2_5num_samp2)
  Shared_list <- list(distal3_5num_samp1,
                      proximal3_5num_samp1,
                      distal3_5num_samp2,
                      proximal3_5num_samp2)
  Referencespecific_list <- list(distal4_5num_samp1,
                                 proximal4_5num_samp1,
                                 distal4_5num_samp2,
                                 proximal4_5num_samp2)

  if (is.null(xlabel)) {
  categ <- c(paste0(mysamps[1],' TSS-distal'),
	paste0(mysamps[1],' TSS-proximal'),
        paste0(mysamps[2],' TSS-distal'),
	paste0(mysamps[2],' TSS-proximal'))
  }
  else(categ <- xlabel)

  explabel <- paste0(nonreference, "-specific (by intensity)")
  reflabel <- paste0(reference, "-specific (by intensity)")

    p <- highchart(width = 750, height = 750 ) %>%
      hc_title(text = maintitle,
               style = list(color = '#2E1717',
                            fontWeight = 'bold',
                            fontSize = maintitlesize)) %>%
      hc_plotOptions(
        boxplot = list(
          fillColor = '#ffffff',
          lineWidth = 2,
          medianColor = '#000000',
          medianWidth = 2,
          stemColor = '#000000',
          stemDashStyle = 'dot',
          stemWidth = 1,
          whiskerColor = '#000000',
          whiskerLength = '20%',
          whiskerWidth = 3
        )
      ) %>%
      hc_add_series(data = Experimentspecific_list,
                    fillColor = cols[1],
                    name = explabel,
                    type = "boxplot") %>%
      hc_add_series(data = Ambiguous_list,
                    fillColor = cols[2],
                    name = 'Ambiguous',
                    type = "boxplot") %>%
      hc_add_series(data = Shared_list,
                    fillColor = cols[3],
                    name = 'Shared',
                    type = "boxplot") %>%
      hc_add_series(data = Referencespecific_list,
                    fillColor = cols[4],
                    name = reflabel,
                    type = "boxplot") %>%
      hc_yAxis(title = list(text = ylabel,
                            style = list(fontSize = ylabelsize)),
               labels = list(format = "{value}")) %>%
      hc_xAxis(categories = categ,
               labels = list(style = list(fontSize = xlabelsize))) %>%
      hc_tooltip(valueDecimals = 2) %>%
      hc_colors(cols) %>%
      hc_exporting(enabled = TRUE)
    return(p)
    }


#' Plots a pie that compares altered regions as determined by peak
#' presence or by #' differential counts.  The type of regulatory region
#' (TSS-proximal, TSS-distal, or both) and type of peak comparison
#' (intensity or peak) must be specified.
#'
#' @param analysisresultsmatrix analysisresults of Intensity analysis place into
#' analysisresults matrix by the analyzeanalysisresults function
#' @param region pick a region, regions can be 'TSS-distal', 'TSS-proximal',
#' or 'both' -- INCLUDE quotes
#' @param method pick a method, methods can be 'Intensity' or 'Peak'
#' include quotes
#' @param palette RColorBrewer palette to change graph colors
#' @param maintitle main title (default generated from sample names)
#' @param maintitlesize main title size (default, 20px)
#' @return pie chart
#' @examples
#' \dontrun{
#' csvfile <- loadCSVFile("DNAseEncodeExample.csv")
#' samplePeaks <- loadBedFiles(csvfile)
#' consensusPeaks <- getConsensusPeaks(samplepeaks = samplePeaks,
#' minreps = 2)
#' TSSannot <- getTSS()
#' consensusPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consensusPeaks,
#' TSS = TSSannot,
#' merge = TRUE,
#' regionspecific = TRUE,
#' distancefromTSSdist = 1500,
#' distancefromTSSprox = 1000)
#' consensusPeaksCounts <- getCounts(annotpeaks = consensusPeaksAnnotated,
#'                                  sampleinfo = csvfile,
#'                                  reference = 'SAEC',
#'                                  chrom = 'chr21')
#' alteredPeaks <- countanalysis(counts = consensusPeaksCounts,
#' pval = 0.01,
#' lfcvalue = 1)
#' alteredPeaksCategorized <- categAltrePeaks(alteredPeaks,
#'                                           lfctypespecific = 1.5,
#'                                           lfcshared = 1.2,
#'                                           pvaltypespecific = 0.01,
#'                                           pvalshared = 0.05)
#' comparePeaksAnalysisResults <- comparePeaksAltre(alteredPeaksCategorized)
#' plotCompareMethods(comparePeaksAnalysisResults)
#'}

plotCompareMethods <- function(analysisresultsmatrix,
                               region = "both",
                               method = "Intensity",
                               palette = NULL,
                               maintitle = NULL,
                               maintitlesize = "16px") {



    analysisresultsmatrix <- as.matrix(analysisresultsmatrix$compareresults[,c("peak","intensity")])

    if (!is.null(palette)) {
    cols <- RColorBrewer::brewer.pal(3, palette)
  }
  else{
    cols <- c("#C71585", "#00E5EE", "#000080", "#d3d3d3")
  }

  if (region == "TSS-proximal") {
    feature <- c("TSS-proxs")
    coordinates <- c(2, 5, 8)
  }
  if (region == "TSS-distal") {
    feature <- c("TSS-dists")
    coordinates <- c(1, 4, 7)
  }
  if (region == "both") {
    region <- c("All")
    feature <- c("TSS-dists", "TSS-proxs")
    coordinates <- c(3, 6, 9)
  }
  # identifies the correct numbers from the
  # analysisresults matrix based on the
  # regulatory region of interest
  if (method == "Intensity") {
    case <- analysisresultsmatrix[coordinates[1], 1]
    reference <- analysisresultsmatrix[coordinates[2], 1]
    shared <- analysisresultsmatrix[coordinates[3], 1]
  }

  if (method == "Peak") {
    case <- analysisresultsmatrix[coordinates[1], 2]
    reference <- analysisresultsmatrix[coordinates[2], 2]
    shared <- analysisresultsmatrix[coordinates[3], 2]
  }
  # identifies the correct numbers from the
  # analysisresults matrix based on the
  # method of region
  # string <- paste(
  #   rownames(analysisresultsmatrix)[1],
  #   rownames(analysisresultsmatrix)[2],
  #   rownames(analysisresultsmatrix)[3],
  #   rownames(analysisresultsmatrix)[4],
  #   rownames(analysisresultsmatrix)[5],
  #   rownames(analysisresultsmatrix)[6],
  #   rownames(analysisresultsmatrix)[7],
  #   rownames(analysisresultsmatrix)[8],
  #   rownames(analysisresultsmatrix)[9]
  # )

  #stringsplit <- strsplit(string, " ")
  #uniquestringsplit <- unique(stringsplit[[1]])
  split <-
    unlist(strsplit(rownames(analysisresultsmatrix)[1], split = " "))
  names <- split[!(split %in% c("TSS-dists"))]
  names <- paste(names, collapse = " ")
  casename <- names

  split <-
    unlist(strsplit(rownames(analysisresultsmatrix)[4], split = " "))
  names <- split[!(split %in% c("TSS-dists"))]
  names <- paste(names, collapse = " ")
  referencename <- names

  # this is a way to the name of the 'case'
  # from the analysisresults matrix

  if (is.null(maintitle)) {
    mtitle <- paste(region, method)
  } else {
    mtitle <- maintitle
    }

  p <- highchart() %>%
    hc_chart(type = "pie") %>%
    hc_title(
      text = mtitle,
      style = list(color = '#2E1717',
                   fontWeight = 'bold',
                   fontSize = maintitlesize)) %>%
    hc_plotOptions(series = list(showInLegend = TRUE)) %>%
    hc_legend(
      enabled = TRUE,
      layout = "horizontal",
      align = "center",
      verticalAlign = "bottom",
      floating = FALSE,
      maxHeight = 100,
      x = 0,
      y = 16
    ) %>%
    hc_add_series(
      data = list(
        list(
        y = case,
        name = casename,
        dataLabels = FALSE
      ),
      list(
        y = reference,
        name = referencename,
        dataLabels = FALSE
      ),
      list(
        y = shared,
        name = "Shared",
        dataLabels = FALSE
      )
    ),
    name = paste(region, method)
    ) %>%
    hc_colors(cols)  %>%
    hc_exporting(enabled = TRUE)
  return(p)
}

#' Plots pie charts for comparison of two methods of identifying altered
#' regulatory regions.  Makes pie charts for TSS-proximal, TSS-distal, and
#' combined for both intensity-based peaks and for peaks identified by hotspot
#' calling algorithms.  There is no return value. Six pie charts swill be
#' plotted.
#' @param analysisresultsmatrix analysisresults of countanalysis function
#' place into a a analysisresults matrix by the analyzeanalysisresults function
#' @param viewer whether the plot should be displayed in the RStudio viewer or
#' in Shiny/Knittr
#' @param palette RColorBrewer palette to change graph colors
#' @param title11 title of the first graph in the first row
#' @param title12 title of the second graph in the first row
#' @param title13 title of the third graph in the first row
#' @param title21 title of the first graph in the second row
#' @param title22 title of the second graph in the second row
#' @param title23 title of the third graph in the second row
#' @param maintitlesize main title size (default, 20px)
#'
#' @examples
#' \dontrun{
#' csvfile <- loadCSVFile("DNAseEncodeExample.csv")
#' samplePeaks <- loadBedFiles(csvfile)
#' consensusPeaks <- getConsensusPeaks(samplepeaks = samplePeaks,
#' minreps = 2)
#' TSSannot <- getTSS()
#' consensusPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consensusPeaks,
#' TSS = TSSannot,
#' merge = TRUE,
#' regionspecific = TRUE,
#' distancefromTSSdist = 1500,
#' distancefromTSSprox = 1000)
#' consensusPeaksCounts <- getCounts(annotpeaks = consensusPeaksAnnotated,
#'                                  sampleinfo = csvfile,
#'                                  reference = 'SAEC',
#'                                  chrom = 'chr21')
#' alteredPeaks <- countanalysis(counts = consensusPeaksCounts,
#' pval = 0.01,
#' lfcvalue = 1)
#' alteredPeaksCategorized <- categAltrePeaks(alteredPeaks,
#'                                           lfctypespecific = 1.5,
#'                                           lfcshared = 1.2,
#'                                           pvaltypespecific = 0.01,
#'                                           pvalshared = 0.05)
#' comparePeaksAnalysisResults <- comparePeaksAltre(alteredPeaksCategorized)
#' plotCompareMethodsAll(comparePeaksAnalysisResults)
#' }
#' @export
#'
plotCompareMethodsAll <- function(analysisresultsmatrix,
                                  viewer = TRUE,
                                  palette = NULL,
                                  title11 = NULL,
                                  title12 = NULL,
                                  title13 = NULL,
                                  title21 = NULL,
                                  title22 = NULL,
                                  title23 = NULL,
                                  maintitlesize = "20px"
                                 ) {


#    if (is.matrix(analysisresultsmatrix[[1]]) == FALSE) {
#      stop("The input is not a matrix!")
#    }

    p1 <- plotCompareMethods(analysisresultsmatrix,
                             "TSS-proximal",
                             "Intensity",
                             palette = palette,
                             maintitle = title11,
                             maintitlesize = maintitlesize)
    p2 <- plotCompareMethods(analysisresultsmatrix,
                             "TSS-distal",
                             "Intensity",
                             palette = palette,
                             maintitle = title12,
                             maintitlesize = maintitlesize)
    p3 <- plotCompareMethods(analysisresultsmatrix,
                             "both",
                             "Intensity",
                             palette = palette,
                             maintitle = title13,
                             maintitlesize = maintitlesize)
    p4 <- plotCompareMethods(analysisresultsmatrix,
                             "TSS-proximal",
                             "Peak",
                             palette = palette,
                             maintitle = title21,
                             maintitlesize = maintitlesize)
    p5 <- plotCompareMethods(analysisresultsmatrix,
                             "TSS-distal",
                             "Peak",
                             palette = palette,
                             maintitle = title22,
                             maintitlesize = maintitlesize)
    p6 <- plotCompareMethods(analysisresultsmatrix,
                             "both",
                             "Peak",
                             palette = palette,
                             maintitle = title23,
                             maintitlesize = maintitlesize)

    if (viewer == TRUE) {
      p <- htmltools::browsable(hw_grid(p1, p2, p3, p4, p5, p6, ncol = 3,
                                        rowheight = 300))
    }
    else {
      p <- hw_grid(p1, p2, p3, p4, p5, p6, ncol = 3)
    }
    return(p)
  }



##############################################################################

#' Given the output from processPathways(), creates a heatmap from
#' the ouput of the GREAT enrichment analysis. Presence or absence of
#' the pathway in enrichment of both type-specific (increased or decreased
#' log2fold change, low p-value) and shared (no change, higher p-value)
#' regulatory regions is plotted.
#'
#' @param input results from GREAT enrichment analysis
#' @param pathwaycateg ontology, to see available ontologies in your input results (e.g. named
#'	GREATpathways, type getOntologies(GREATpathways)
#' @param test character, "Binom" uses binomial test restuls, "Hyper" uses
#'      hypergeometric test results.  Default is "Binom"
#' @param numshow number of top pathways (ranked according to p-value) of each type
#' 	(expt, reference, shared) to show in the plot (default=10)
#' @param maintitle main title (default, "GREAT Enrichment Analysis")
#' @param maintitlesize main title size (default, 20px)
#' @param subtitle subtitle (default, "color corresponds to p-value")
#' @param subtitlesize subitle size (default 15px)
#' @param xlabelsize size of xlabel (default, 10px)
#' @param ylabelsize size of ylabel (default, 10px)
#' @param xlabel label for x-axis (default, Experiment-specific, shared, Reference-specific )
#' @return heatmap
#'
#' @examples
#' \dontrun{
#' csvfile <- loadCSVFile("DNAseEncodeExample.csv")
#' samplePeaks <- loadBedFiles(csvfile)
#' consensusPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps = 2)
#' TSSannot <- getTSS()
#' consensusPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consensusPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           distancefromTSSdist = 1500,
#'                                           distancefromTSSprox = 1000)
#' consensusPeaksCounts <- getCounts(annotpeaks = consensusPeaksAnnotated,
#'                               reference = 'SAEC',
#'                               sampleinfo = csvfile,
#'                               chrom = 'chr21')
#' alteredPeaks <- countanalysis(counts=consensusPeaksCounts,
#'                              pval=0.01,
#'                              lfcvalue=1)
#' alteredPeaksCategorized <- categAltrePeaks(analysisresults = alteredPeaks,
#'                              lfctypespecific = 1.5,
#'                              lfcshared = 1.2,
#'                              pvaltypespecific = 0.01,
#'                              pvalshared = 0.05)
#' callPaths <- runGREAT(peaks = alteredPeaksCategorized)
#' pathResults <- processPathways(callPaths, pathway_category = "GO",
#' enrichcutoff = 2, adjpvalcutoff = 0.05)
#' plotGREATenrich(pathResults, maintitle = "GREAT Enrichment Analysis",
#' pathwaycateg = "GO_Molecular_Function")
#'}
#' @export
plotGREATenrich <- function(input,
                            maintitle = "GREAT Enrichment Analysis",
                            pathwaycateg = NULL,
                            test = "Binom",
                            numshow = 10,
                            maintitlesize = "20px",
                            ylabelsize = "10px",
                            xlabelsize = "10px",
                            xlabel = NULL,
                            subtitle = "(color corresponds to p-value)",
                            subtitlesize = "13px") {


  variable = value = Experiment_specific = Reference_specific = Shared = c()

  if (is.null(pathwaycateg)) {
    stop("Please designate a pathway with the parameter pathwaycateg")
  }

  if (is.list(input) == FALSE) {
    stop(
      "The input is not a list! Please make sure you are
      using the output from the enrichment analysis"
    )
  }

  if (is.na(match(test, c("Hyper", "Binom")))) {
    stop("test must be either 'Hyper' or 'Binom'")
  }

  mycols = c("name",
             paste0(test, "_Fold_Enrichment"),
             paste0(test, "_adj_PValue"))

  if (is.list(input$ExperimentSpecificByIntensity$Sig_Pathways) == FALSE |
      is.list(input$ReferenceSpecificByIntensity$Sig_Pathways) == FALSE |
      is.list(input$Shared$Sig_Pathways) == FALSE |
      length(input) != 3 |
      length(which(!is.na(match(
        mycols,
        colnames(input$ExperimentSpecificByIntensity$Sig_Pathways[[pathwaycateg]])
      )))) !=
      length(mycols) |
      length(which(!is.na(match(
        mycols,
        colnames(input$ReferenceSpecificByIntensity$Sig_Pathways[[pathwaycateg]])
      )))) !=
      length(mycols) |
      length(match(mycols,
                   colnames(input$Shared$Sig_Pathways[[pathwaycateg]]))) !=
      length(mycols) |
      all(
        names(input) != c(
          "ExperimentSpecificByIntensity",
          "ReferenceSpecificByIntensity",
          "Shared"
        )
      )) {
    stop(
      "The input is not a list of three dataframes or there are no enriched pathways to plot.
      Be sure the input is the output from running processPathways(()"
      )
  }

  up <-
    input$ExperimentSpecificByIntensity$Sig_Pathways[[pathwaycateg]][, mycols]

  if (is.null(nrow(up))) {
    up$name <- NA
  } else {
    if (nrow(up) > numshow) {
      # order by last row, which is always adjusted p-value
      up <- up[order(up[, 3])[1:numshow], ]
    }
  }

  reference <-
    input$ReferenceSpecificByIntensity$Sig_Pathways[[pathwaycateg]][,mycols]
  if (is.null(nrow(reference))) {
    reference$name <- NA
  } else {
    if (nrow(reference) > numshow) {
      # order by last row, which is always adjusted p-value
      reference <- reference[order(reference[, 3])[1:numshow], ]
    }
  }

  shared <- input$Shared$Sig_Pathways[[pathwaycateg]][, mycols]
  if (is.null(nrow(shared))) {
    shared$name <- NA
  } else {
    if (nrow(shared) > numshow) {
      # order by last row, which is always adjusted p-value
      shared <- shared[order(shared[, 3])[1:numshow], ]
    }
  }

  # make a list of all the pathways in up, down, and shared
  pathways <- unique(c(up$name,
                       reference$name,
                       shared$name))
  pathways <- pathways[!is.na(pathways)]
  if (is.na(pathways) || length(pathways) == 0) {
    stop("No pathways are significant!")
  }

  # make a matrix with as many row as there are pathways
  heatmapmatrix <- matrix(data = NA,
                          nrow = length(pathways),
                          ncol = 3)
  # name the rows with the pathway names
  row.names(heatmapmatrix) <- pathways

  colnames(heatmapmatrix) <-
    c("Experiment_specific", "Reference_specific", "Shared")

  # places the adjusted p-value in the matrix is there is one
  for (i in 1:nrow(heatmapmatrix)) {
    if (row.names(heatmapmatrix)[i] %in% up$name) {
      num1 <- which(up$name == row.names(heatmapmatrix)[i])
      heatmapmatrix[i, 1] <- up[num1, mycols[3]]
    }

    if (row.names(heatmapmatrix)[i] %in% reference$name) {
      num2 <- which(reference$name == row.names(heatmapmatrix)[i])
      heatmapmatrix[i, 2] <- reference[num2, mycols[3]]
    }

    if (row.names(heatmapmatrix)[i] %in% shared$name) {
      num3 <- which(shared$name == row.names(heatmapmatrix)[i])
      heatmapmatrix[i, 3] <- shared[num3, mycols[3]]
    }
  }

  # Create a data.frame of the heatmapmatrix and sort
  heatmapdata <- as.data.frame(heatmapmatrix)
  heatmapdata <- heatmapdata[order(
    heatmapdata$Reference_specific,
    heatmapdata$Experiment_specific,
    heatmapdata$Shared,
    decreasing = TRUE
  ),]
  # Create ids:
  heatmapdata$id <- rownames(heatmapdata)
  rownames(heatmapdata) <- c(1:nrow(heatmapdata))

  #suppressMessages(meltedheatmapdata <- reshape2::melt(heatmapdata))
  suppressMessages(
    meltedheatmapdata <- tidyr::gather(
      heatmapdata,
      variable,
      value,
      Experiment_specific,
      Reference_specific,
      Shared
    )
  )

  meltedheatmapdata$newid <-
    stringr::str_wrap(meltedheatmapdata$id, width = 80)

  meltedheatmapdata$id <- factor(meltedheatmapdata$id,
                                 levels = unique(meltedheatmapdata$id))
  #all possible values of X (type) and Y (pathways)
  theXAxis <- as.character(meltedheatmapdata[, 2])
  theYAxis <- meltedheatmapdata[, 4]

  #unique values of X and Y
  theUniqueY <- unique(meltedheatmapdata$newid)
  theUniqueX <-
    c("Experiment_specific", "Shared", "Reference_specific")

  # Substitute words with position on the meatrix
  for (i in 0:(length(theUniqueY) - 1))
  {
    num <- which(theYAxis == theUniqueY[i + 1])
    theYAxis[num] <- i
  }
  for (i in 0:(length(theUniqueX) - 1))
  {
    num <- which(theXAxis == theUniqueX[i + 1])
    theXAxis[num] <- i
  }

  #create final formatting
  dataforHeatmap <- as.data.frame(cbind(
    as.numeric(theXAxis),
    as.numeric(theYAxis),
    as.numeric(meltedheatmapdata$value)
#as.numeric(format(meltedheatmapdata$value,scientific=T,digits=2))
  ))

  formattedHeatmapData <- list_parse2(dataforHeatmap)

  fntltp <- JS(
    "function(){
    return 'pval='+this.point.value;
    }"
   )

  if (is.null(xlabel))
  {categ = c("Experiment-specific", "Shared", "Reference-specific")}
  else
  {categ = xlabel}

  p <- highchart(width = 800, height = 700) %>%
    hc_chart(type = "heatmap", spacingRight = 160) %>%
    hc_title(text = maintitle,
             style = list(color = '#2E1717',fontSize = maintitlesize,
                          fontWeight = 'bold')) %>%
    hc_subtitle(text = subtitle,
                style = list(fontSize = subtitlesize)) %>%
    hc_xAxis(categories = categ,
             labels = list(style = list(fontSize = xlabelsize))) %>%
    hc_yAxis(categories = theUniqueY, labels = list(style = list(fontSize = ylabelsize))) %>%
    hc_add_series(name = "matrix location, p-value",
                  data = formattedHeatmapData) %>%
    hc_tooltip(formatter = fntltp, valueDecimals = 2) %>%
    hc_colorAxis(stops = color_stops(2, colors = c("#5097D1", "#DEEFF5")),
                 min = min(as.numeric(dataforHeatmap[ , 3]), na.rm = T),
                 max = max(as.numeric(dataforHeatmap[ , 3]), na.rm = T)) %>%
    hc_legend(
      enabled = TRUE,
      layout = "vertical",
      align = "right",
      verticalAlign = "top",
      floating = FALSE,
      maxWidth = 200,
      x = -10, # 90
      y = 100, # 70
      padding = 2
      #title = list(text="p-value")
    ) %>%
    hc_exporting(enabled = TRUE)
  #create final formatting
  return(p)
  } # end plotGREATenrich


