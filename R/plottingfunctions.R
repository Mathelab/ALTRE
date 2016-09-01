#' Given the output from getConsensusPeaks, generate a barplot
#' of countstatistics
#'
#' @param samplepeaks output generated from getConsensusPeaks
#' @param palette RColorBrewer palette to change graph colors
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
plotConsensusPeaks <- function(samplepeaks, palette = "Set1") {

  cols <- RColorBrewer::brewer.pal(3, palette)


  CellType <- NULL
  #without this R CMD check throws no visible binding for global variable error

  consPeaksStats <- samplepeaks$consPeaksStats
  row.names(consPeaksStats) <- NULL
  consPeaksStats[ , 2] <-  as.numeric(as.character(consPeaksStats[[2]]))
  consPeaksStats[ , 3] <-  as.numeric(as.character(consPeaksStats[[3]]))
  statsFormated <- tidyr::gather(consPeaksStats, "CellType", "Count", 2:3)

  plottingData <- statsFormated %>%
    split(levels(as.factor(statsFormated$CellType)))

  col1 <- cols[1]
  col2 <- cols[2]

  p <- highchart(height = 700) %>%
    hc_title(text = "Peak Counts by Cell Type",
             style = list(color = '#2E1717',
                          fontWeight = 'bold')) %>%
    hc_subtitle(text = "For bioreplicates and their merged consensus track") %>%
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
    hc_yAxis(title = list(text = "Peak Counts"),
             labels = list(format = "{value}")) %>%
    hc_xAxis(categories = plottingData[[1]]$PeakType) %>%
    hc_legend(
      enabled = TRUE,
      layout = "horizonal",
      align = "center",
      verticalAlign = "bottom",
      floating = FALSE,
      maxHeight = 100,
      x = 15,
      y = 16
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
#' plotCombineAnnotatePeaks(consensusPeaksAnnotated)
#' }
#' @export
#'
plotCombineAnnotatePeaks <- function(conspeaks, viewer = TRUE,
                                     palette = "Set1") {

    cols <- RColorBrewer::brewer.pal(3, palette)
    CellType <- NULL
    #R CMD check throws no visible binding for global variable error

    mergeStats <- conspeaks$mergestats
    row.names(mergeStats) <- NULL
    mergeStats[ , 2] <-  as.numeric(as.character(mergeStats[[2]]))
    mergeStats[ , 3] <-  as.numeric(as.character(mergeStats[[3]]))
    mergeStatsFormatted <- tidyr::gather(mergeStats, "CellType", "Count", 2:3)


    if ( nrow(mergeStatsFormatted) == 1 ) {
        stop("No plot to show since merging was not performed
             when calling combineAnnotatePeaks function")
    }


    feature <- "TotalNumber"
    if ( feature == "TotalNumber") {
        mergeStatsTotal <- dplyr::filter(mergeStatsFormatted,
                                         CellType == "TotalNumber")
        thecondition <- matrix(unlist(strsplit(mergeStatsTotal$Condition, "_")),
                               nrow = 3, ncol = 4)[2,]
        mergeStatsBefore <- dplyr::filter(mergeStatsTotal,
                                          thecondition == "before")
        mergeStatsAfter <- dplyr::filter(mergeStatsTotal,
                                         thecondition == "after")

        p1 <- highchart(height = 400) %>%
            hc_title(text = "Number of REs",
                     style = list(color = '#2E1717',
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
                )) %>%
            hc_add_series(
                data = mergeStatsAfter$Count,
                name = c("TSS-proximal"),
                type = "column",
                dataLabels = list(
                    enabled = TRUE,
                    rotation = 270,
                    color = '#FFFFFF',
                    y = 40
                )) %>%
            hc_yAxis(title = list(text = "Number of REs"),
                     labels = list(format = "{value}")) %>%
            hc_xAxis(categories = c("Before Merging", "After Merging")) %>%
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
    if ( feature == "MeanLength" ) {
      mergeStatsMean <- dplyr::filter(mergeStatsFormatted,
                                      CellType == "MeanLength")
      thecondition <- matrix(unlist(strsplit(mergeStatsMean$Condition, "_")),
                             nrow = 3,
                             ncol = 4)[2, ]
      mergeStatsBefore <- dplyr::filter(mergeStatsMean,
                                        thecondition == "before")
      mergeStatsAfter <- dplyr::filter(mergeStatsMean,
                                       thecondition == "after")


        p2 <- highchart(height = 400) %>%
            hc_title(text = "Mean length of REs",
                     style = list(color = '#2E1717',
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
                )) %>%
            hc_add_series(
                data = mergeStatsAfter$Count,
                name = c("TSS-proximal"),
                type = "column",
                dataLabels = list(
                    enabled = TRUE,
                    rotation = 270,
                    color = '#FFFFFF',
                    y = 40
                )) %>%
            hc_yAxis(title = list(text = "Mean Length of REs"),
                     labels = list(format = "{value}")) %>%
            hc_xAxis(categories = c("Before Merging", "After Merging")) %>%
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
    p <- htmltools::browsable(hw_grid(p1, p2, ncol = 2, rowheight = 550))
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

plotGetCounts <- function(countsConsPeaks, palette = "Set1") {
  region <- NULL
  variable <- NULL
  #set to null forR CMD Check error: Undefined global functions/variables

  cols <- RColorBrewer::brewer.pal(4, palette)

  mydf <- countsConsPeaks$regioncountsforplot
  varstack <- suppressMessages(reshape2::melt(mydf))
  TSSdistal_A549 <- dplyr::filter(varstack, region == "TSS-distal"
                                  & variable == "A549")
  TSSproximal_A549 <- dplyr::filter(varstack, region == "TSS-proximal"
                                    & variable == "A549")
  TSSdistal_SAEC <- dplyr::filter(varstack, region == "TSS-distal"
                                  & variable == "SAEC")
  TSSproximal_SAEC <- dplyr::filter(varstack, region == "TSS-proximal"
                                    & variable == "SAEC")

  p <- hchart(stats::density(TSSdistal_A549$value), area = TRUE,
              name = "A549 TSS-distal") %>%
    hc_title(text = "Density of log2 read counts
             (normalized by library and region sizes)",
             style = list(color = '#2E1717',
                          fontWeight = 'bold')) %>%
    hc_yAxis(title = "density") %>%
    hc_xAxis(title = "log2 read counts") %>%
    hc_add_series_density(stats::density(TSSproximal_A549$value),
                          area = TRUE, name = "A549 TSS-proximal") %>%
    hc_add_series_density(stats::density(TSSdistal_SAEC$value),
                          area = TRUE, name = "SAEC TSS-distal") %>%
    hc_add_series_density(stats::density(TSSproximal_SAEC$value),
                          area = TRUE, name = "SAEC TSS-proximal") %>%
    hc_colors(cols) %>%
    hc_exporting(enabled = TRUE)
  return(p)
}

#'
#' #' Create a volcano plot from the output of categAltrePeaks
#' #'
#' #' @param altrepeakscateg output generated from countanalysis() then
#' #' categAltrePeaks()
#' #' @param viewer whether the plot should be displayed in the RStudio viewer or
#' #'        in Shiny/Knittr
#' #' @param palette RColorBrewer palette to change graph colors
#'
#' #'
#' #' @return a highcharter object
#' #'
#' #' @examples
#' #'
#' #' csvfile <- loadCSVFile("DNAseEncodeExample.csv")
#' #' samplePeaks <- loadBedFiles(csvfile)
#' #' consensusPeaks <- getConsensusPeaks(samplepeaks = samplePeaks,
#' #' minreps = 2)
#' #' TSSannot <- getTSS()
#' #' consensusPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consensusPeaks,
#' #' TSS = TSSannot,
#' #' merge = TRUE,
#' #' regionspecific = TRUE,
#' #' distancefromTSSdist = 1500,
#' #' distancefromTSSprox = 1000)
#' #' consensusPeaksCounts <- getCounts(annotpeaks = consensusPeaksAnnotated,
#' #'                                  sampleinfo = csvfile,
#' #'                                  reference = 'SAEC',
#' #'                                  chrom = 'chr21')
#' #' alteredPeaks <- countanalysis(counts = consensusPeaksCounts,
#' #' pval = 0.01,
#' #' lfcvalue = 1)
#' #' alteredPeaksCategorized <- categAltrePeaks(alteredPeaks,
#' #'                                           lfctypespecific = 1.5,
#' #'                                           lfcshared = 1.2,
#' #'                                           pvaltypespecific = 0.01,
#' #'                                           pvalshared = 0.05)
#' #' plotCountAnalysis(alteredPeaksCategorized)
#' #'
#' #' @export

# plotCountAnalysis <- function(altrepeakscateg, viewer = TRUE, palette = NULL ) {
#
#
#   if ( !is.null(palette) ) {
#     cols <- RColorBrewer::brewer.pal(4, palette) }
#   else{cols <- c("#d3d3d3", "#C71585", "#00E5EE","#000080")}
#                         #grey (ambiguous)
#                         #magenta (experiment-specific)
#                         #blue (reference specific)
#   #blue (shared)
#
#   log2FoldChange <- NULL
#   padj <- NULL
#   REaltrecateg <- NULL
#   #To prevent R CMD check error
#
#   colnames(altrepeakscateg$analysisresults)
#
#   toplot <- altrepeakscateg$analysisresults[ ,c("region",
#                                                 "log2FoldChange",
#                                                 "padj",
#                                                 "REaltrecateg")]
#   tssdist <- toplot[which(toplot$region == "TSS-distal"), ]
#   tssdist$padj <- round(-log10(tssdist$padj), 2)
#   tssdist$log2FoldChange <- round(tssdist$log2FoldChange, 2)
#   tssprox <- toplot[which(toplot$region == "TSS-proximal"), ]
#   tssprox$padj <- round(-log10(tssprox$padj), 2)
#   tssprox$log2FoldChange <- round(tssprox$log2FoldChange, 2)
#   lengthRE <- rep("", length(tssdist$REaltrecateg))
#
#   p1 <- highchart() %>%
#     hc_chart(type = "scatter") %>%
#     hc_title(text = "TSS-distal",
#              style = list(color = '#2E1717',
#                           fontWeight = 'bold')) %>%
#     hc_add_series_df(data = tssdist, x = log2FoldChange, y = padj,
#                      type = "scatter", group = REaltrecateg)  %>%
#     hc_xAxis(title = list(text = "log2fold change")) %>%
#     hc_yAxis(title = list(text = "-log10 pvalue")) %>%
#     hc_tooltip(headerFormat = "",
#                pointFormat  = "<b>log2FC</b> = {point.x}<br> <b>-log10pvalue</b>
#                = {point.y}<br>") %>%
#     hc_colors(cols) %>%
#     hc_exporting(enabled = TRUE)
#
#   lengthRE <- rep("", length(tssprox$REaltrecateg))
#
#   num1 <- min(which(tssprox$REaltrecateg == "Experiment Specific"))
#   num2 <- min(which(tssprox$REaltrecateg == "Reference Specific"))
#   num3 <- min(which(tssprox$REaltrecateg == "Shared"))
#   num4 <- min(which(tssprox$REaltrecateg == "Ambiguous"))
#
#   lengthRE[num1] <- "Experiment Specific"
#   lengthRE[num2] <- "Reference Specific"
#   lengthRE[num3] <- "Shared"
#   lengthRE[num4] <- "Ambiguous"
#
#
#   p2 <- highchart() %>%
#     hc_chart(type = "scatter") %>%
#     hc_title(text = "TSS-proximal",
#              style = list(color = '#2E1717',
#                           fontWeight = 'bold')) %>%
#     hc_add_series_df(data = tssprox, x = log2FoldChange, y = padj,
#                      type = "scatter", group = REaltrecateg)  %>%
#     hc_xAxis(title = list(text = "log2fold change")) %>%
#     hc_yAxis(title = list(text = "-log10 pvalue")) %>%
#     hc_tooltip(headerFormat = "",
#                pointFormat  = "<b>log2FC</b> = {point.x}<br> <b>-log10pvalue</b>
#                = {point.y}<br>") %>%
#     hc_colors(cols) %>%
#     hc_exporting(enabled = TRUE)
#
#
#   if (viewer == TRUE) {
#     p <- htmltools::browsable(hw_grid(p1, p2, ncol = 2, rowheight = 700))
#   }
#   else {
#     p <- hw_grid(p1, p2, ncol = 2)
#   }
#   return(p)
# }
#


#' Create a volcano plot from the output of categAltrePeaks
#'
#' @param altrepeakscateg output generated from countanalysis() then
#' categAltrePeaks()
#' @param viewer whether the plot should be displayed in the RStudio viewer or
#'         in Shiny/Knittr
#' @param palette RColorBrewer palette to change graph colors
#'
#' @return a highcharter object
#'
#' @export
#'
plotCountAnalysisTemp <- function(altrepeakscateg = NULL, viewer = TRUE, palette = NULL) {

    if ( !is.null(palette) ) {
      cols <- RColorBrewer::brewer.pal(4, palette) }
    else{cols <- c("#d3d3d3", "#C71585", "#00E5EE","#000080")}
                          #grey (ambiguous)
                          #magenta (experiment-specific)
                          #blue (reference specific)

  cat(altrepeakscateg)

  dataraw <-
  ' DOK6 0.51 1.861e-08 0.0003053
  TBX5 -2.129 5.655e-08 0.0004191
  SLC32A1 0.9003 7.664e-08 0.0004191
  IFITM1 -1.687 3.735e-06 0.006809
  NUP93 0.3659 3.373e-06 0.006809
  EMILIN2 1.534 2.976e-06 0.006809
  TPX2 -0.9974 2.097e-06 0.006809
  LAMA2 -1.425 2.39e-06 0.006809
  CAV2 -1.052 3.213e-06 0.006809
  TNN -1.658 8.973e-06 0.01472
  POU3F4 1.181 1.062e-05 0.01584
  COL13A1 -1.647 1.394e-05 0.01592
  IFITM3 -1.61 1.202e-05 0.01592
  SHISA3 -1.477 1.31e-05 0.01592
  LOC285954 1.05 1.456e-05 0.01592
  VEPH1 1.137 2.211e-05 0.02267
  ARHGAP29 -1.526 3.675e-05 0.03547
  KIAA1755 -1.562 3.972e-05 0.0362
  LAMC3 -1.563 4.29e-05 0.03704
  ITM2A -1.398 4.972e-05 0.04078
  DTHD1 1.54 5.594e-05 0.04371
  RBMS1 -0.9139 6.688e-05 0.04988
  CEBPD -1.202 7.859e-05 0.05606
  DMBX1 -1.425 9.529e-05 0.06486
  PAPLN -1.253 9.883e-05 0.06486
  ADM -1.357 0.0001089 0.06872
  COL2A1 -1.187 0.0001424 0.07794
  HS3ST3A1 -1.004 0.0001388 0.07794
  DYSF -1.03 0.0001425 0.07794
  PI16 1.495 0.0001297 0.07794
  CDC42EP5 -1.355 0.0001581 0.08146
  SLC12A8 -0.9425 0.0001589 0.08146
  ZNF391 -1.024 0.0001913 0.09512
  GALNTL2 1.075 0.0002298 0.1109
  C4orf45 1.288 0.0002472 0.1159
  KIF18B -0.8849 0.0002551 0.1162
  KIF20A -0.9505 0.0002972 0.1318
  PDE1B 1.053 0.0003356 0.1449
  BCAN 1.117 0.0003698 0.1477
  APLNR -1.365 0.000378 0.1477
  CILP -1.11 0.0003582 0.1477
  TEC -1.373 0.0003701 0.1477
  KLF5 -0.8177 0.0004159 0.1578
  ACSS2 -0.5578 0.0004232 0.1578
  RAPGEF2 0.3371 0.0004513 0.1645
  C1orf51 -0.6451 0.0005237 0.1665
  IGF2-AS1 -1.235 0.0005835 0.1665
  RPLP0P2 -1.096 0.0005689 0.1665
  COTL1 -0.7376 0.0005886 0.1665
  MYO1D -0.8454 0.0005529 0.1665
  CIAO1 0.2695 0.000522 0.1665
  POU3F3 0.6857 0.000588 0.1665
  CFLAR -0.9694 0.0005598 0.1665
  BHLHE40 -1.127 0.0004785 0.1665
  PLSCR4 -1.317 0.0004978 0.1665
  HECW1 0.5135 0.0005373 0.1665
  KCNQ3 1.147 0.000483 0.1665
  TIMP1 -1.15 0.0005267 0.1665
  CAV1 -1.115 0.0006722 0.1869
  LTBP4 -0.8186 0.0006991 0.1912
  HDAC1 -0.4542 0.0007203 0.1937
  HSPG2 -1.218 0.0007963 0.1947
  CYR61 -1.102 0.0008071 0.1947
  SAP30 -0.02968 0.904 0.9994
  FBXO8 -0.07336 0.6573 0.9994
  CEP44 -0.1132 0.6782 0.9994
  HPGD -0.3122 0.4238 0.9994
  GLRA3 0.3489 0.3221 0.9994
  ADAM29 0.4294 0.2735 0.9994
  WDR17 0.1681 0.4038 0.9994
  SPATA4 0.4892 0.2052 0.9994
  SPCS3 -0.06255 0.7011 0.9994
  NEIL3 -0.4648 0.1156 0.9994
  AGA -0.4363 0.1089 0.9994
  MGC45800 -0.2157 0.4652 0.9994
  ODZ3 -0.3824 0.1282 0.9994
  FAM92A3 -0.4185 0.282 0.9994
  C4orf38 -0.05564 0.8868 0.9994
  CDKN2AIP 0.07392 0.7821 0.9994
  LOC389247 -0.06532 0.8677 0.9994
  ING2 -0.1931 0.1921 0.9994
  RWDD4 0.08427 0.6963 0.9994
  PDLIM3 0.195 0.5488 0.9994
  GLP1R -0.001237 0.9971 0.9996
  TCP1 -0.00206 0.9931 0.9996
  CCM2 -0.000938 0.9955 0.9996
  WBSCR16 0.001843 0.9928 0.9996
  LOC155060 -0.003811 0.9916 0.9996
  C8orf42 0.001161 0.9975 0.9996
  CSGALNACT1 -0.001458 0.9967 0.9996
  PPAPDC1B 0.0008906 0.9974 0.9996
  LINC00535 0.003013 0.9929 0.9996
  VPS28 0.0005436 0.9965 0.9996
  DNAJC25 -0.000678 0.9974 0.9996
  TSC1 0.001054 0.9962 0.9996
  RLIM -0.001962 0.9932 0.9996
  TAF9B -0.002379 0.9915 0.9996
  DDX3Y -0.001194 0.996 0.9996
  ZNRF2 0.0002994 0.9988 0.9998
  ZNF559 0.0002839 0.9994 0.9998'

  datmat <- t(matrix(unlist(strsplit(unlist(strsplit(dataraw,"\n " ))," ")), 5, 100))[,-1]

  dat <- data.frame(log2FoldChange = as.numeric(datmat[ , 2]),
                       pvalue = as.numeric(datmat[ , 3]),
                       padj = as.numeric(datmat[ , 4]),
                       categ = rep(c("A", "B"), c(50, 50)))
  rownames(dat) <- datmat[ , 1]

  ind <- as.logical(stats::rbinom(100, 1, 0.5))
  dat1 <- dat[ind, ]
  dat2 <- dat[!ind, ]
  p1 <- highchart() %>%
         hc_title(text = "TSS-distal",
                  style = list(color = '#2E1717',
                               fontWeight = 'bold')) %>%
    hc_add_series_scatter(x = -log10(dat1$pvalue),
                          y = dat1$log2FoldChange,
                          color = dat1$categ,
                          label = rownames(dat1))  %>%
         hc_xAxis(title = list(text = "log2fold change")) %>%
         hc_yAxis(title = list(text = "-log10 pvalue")) %>%
    hc_tooltip(headerFormat = "",
               pointFormat  = "<b>log2FC</b> = {point.x}<br> <b>-log10pvalue</b> = {point.y}<br>") %>%
    hc_colors(cols) %>%
    hc_exporting(enabled = TRUE)


  p2 <- highchart() %>%
         hc_title(text = "TSS-proximal",
                  style = list(color = '#2E1717',
                               fontWeight = 'bold')) %>%
    hc_add_series_scatter(x = -log10(dat2$pvalue),
                          y = dat2$log2FoldChange,
                          color = dat2$categ,
                          label = rownames(dat2))  %>%
         hc_xAxis(title = list(text = "log2fold change")) %>%
         hc_yAxis(title = list(text = "-log10 pvalue")) %>%
    hc_tooltip(headerFormat = "",
               pointFormat  = "<b>log2FC</b> = {point.x}<br> <b>-log10pvalue</b> = {point.y}<br>") %>%
    hc_colors(cols) %>%
    hc_exporting(enabled = TRUE)


    if (viewer == TRUE) {
      p <- htmltools::browsable(hw_grid(p1, p2, ncol = 2, rowheight = 700))
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
plotDistCountAnalysis <- function(analysisresults, counts, palette = NULL){

  if ( !is.null(palette) ) {
    cols <- RColorBrewer::brewer.pal(4, palette) }
  else{cols <- c("#C71585", "#d3d3d3", "#000080", "#00E5EE")}
           #magenta (experiment-specific) #grey (ambiguous) #blue (shared))
  #blue (reference specific)

  readcounts <- counts$regioncounts
  analysisresults <- analysisresults[[1]]
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
  countsinfo <- as.data.frame(SummarizedExperiment::rowRanges(readcounts))
  countcoord <- paste0(countsinfo$seqnames, countsinfo$start, countsinfo$end)
  analcoord <- paste0(analysisresults$chr,
                     analysisresults$start,
                     analysisresults$stop)

  if (!all.equal(analcoord, countcoord)) {
    stop("The peaks in the analysisresults and counts are not the same")
  }

  PEcateg <- analysisresults$region
  altrecateg <- analysisresults$REaltrecateg

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
                    altrecateg = altrecateg)
  TSSdistal <- dplyr::filter(mydf, PEcateg == "TSS-distal")
  distal1 <- dplyr::filter(mydf, altrecateg == "Experiment Specific")
  distal2 <- dplyr::filter(mydf, altrecateg == "Ambiguous")
  distal3 <- dplyr::filter(mydf, altrecateg == "Shared")
  distal4 <- dplyr::filter(mydf, altrecateg == "Reference Specific")

  TSSproximal <- dplyr::filter(mydf, PEcateg == "TSS-proximal")
  proximal1 <- dplyr::filter(mydf, altrecateg == "Experiment Specific")
  proximal2 <- dplyr::filter(mydf, altrecateg == "Ambiguous")
  proximal3 <- dplyr::filter(mydf, altrecateg == "Shared")
  proximal4 <- dplyr::filter(mydf, altrecateg == "Reference Specific")

  distal1_5num_A549 <- stats::fivenum(distal1$meanlog2FPM.A549)
  proximal1_5num_A549 <- stats::fivenum(proximal1$meanlog2FPM.A549)
  distal1_5num_SAEC <- stats::fivenum(distal1$meanlog2FPM.SAEC)
  proximal1_5num_SAEC <- stats::fivenum(proximal1$meanlog2FPM.SAEC)

  distal2_5num_A549 <- stats::fivenum(distal2$meanlog2FPM.A549)
  proximal2_5num_A549 <- stats::fivenum(proximal2$meanlog2FPM.A549)
  distal2_5num_SAEC <- stats::fivenum(distal2$meanlog2FPM.SAEC)
  proximal2_5num_SAEC <- stats::fivenum(proximal2$meanlog2FPM.SAEC)

  distal3_5num_A549 <- stats::fivenum(distal3$meanlog2FPM.A549)
  proximal3_5num_A549 <- stats::fivenum(proximal3$meanlog2FPM.A549)
  distal3_5num_SAEC <- stats::fivenum(distal3$meanlog2FPM.SAEC)
  proximal3_5num_SAEC <- stats::fivenum(proximal3$meanlog2FPM.SAEC)

  distal4_5num_A549 <- stats::fivenum(distal4$meanlog2FPM.A549)
  proximal4_5num_A549 <- stats::fivenum(proximal4$meanlog2FPM.A549)
  distal4_5num_SAEC <- stats::fivenum(distal4$meanlog2FPM.SAEC)
  proximal4_5num_SAEC <- stats::fivenum(proximal4$meanlog2FPM.SAEC)

  Experimentspecific_list <- list(round(distal1_5num_A549,3),
                                  round(proximal1_5num_A549, 3),
                                  round(distal1_5num_SAEC,3),
                                  round(proximal1_5num_SAEC,3))
  Ambiguous_list <- list(round(distal2_5num_A549,3),
                         round(proximal2_5num_A549,3),
                         round(distal2_5num_SAEC,3),
                         round(proximal2_5num_SAEC,3))
  Shared_list <- list(round(distal3_5num_A549,3),
                      round(proximal3_5num_A549,3),
                      round(distal3_5num_SAEC,3),
                      round(proximal3_5num_SAEC,3))
  Referencespecific_list <- list(round(distal4_5num_A549,3),
                                 round(proximal4_5num_A549,3),
                                 round(distal4_5num_SAEC,3),
                                 round(proximal4_5num_SAEC,3))

  categ <- c('A549-specific TSS-distal', 'A549-specific TSS-proximal',
             'SAEC-specific TSS-distal', 'SAEC-specific TSS-proximal')
  p <- highchart() %>%
    hc_title(text = "Distribution of Normalized Counts",
             style = list(color = '#2E1717',
                          fontWeight = 'bold')) %>%
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
                  name = 'Experiment Specific',
                  type = "boxplot") %>%
    hc_add_series(data = Ambiguous_list,
                  name = 'Ambiguous',
                  type = "boxplot") %>%
    hc_add_series(data = Shared_list,
                  name = 'Shared',
                  type = "boxplot") %>%
    hc_add_series(data = Referencespecific_list,
                  name = 'Reference Specific',
                  type = "boxplot") %>%
    hc_yAxis(title = list(text = "Observations"),
             labels = list(format = "{value}")) %>%
    hc_xAxis(categories = categ, title = "Experiment No.") %>%
    hc_tooltip(headerFormat = "<b>{point.key}</b><br>",
               pointFormat = "{point.y}") %>%
    hc_colors(cols) %>%
    hc_exporting(enabled = TRUE)
    return(p)
}

##############################################################################

#' Given the output from enrichment(), creates a heatmap from
#' the ouput of the enrichment analysis. Presence or absence of
#' the pathway in enrichment of both type-specific (increased or decreased
#' log2fold change, low p-value) and shared (no change, higher p-value)
#' regulatory regions is plotted.
#'
#' @param input results from enrichment analysis
#' @param title title of the heatmap
#' @param pvalfilt p-value cut-off for inclusion in heatmap
#' @param removeonlyshared removes regions that come up signifigant only shared
#' regulatory regions when set to TRUE. Default is FALSE.
#' @param numshow number of top pathways (ranked according to p-value) of each
#' type (expt, reference, shared) to show in the plot (default=10)
#'
#' @return heatmap
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
#' MFenrich <- pathenrich(analysisresults = alteredPeaksCategorized,
#' ontoltype = 'MF',
#' enrichpvalfilt = 0.99)
#' enrichHeatmap(MFenrich, title = "TSS-distal GO:MF", pvalfilt = 0.99)
#' }
#' @export

enrichHeatmap <- function(input,
                          title,
                          pvalfilt = 0.01,
                          removeonlyshared = FALSE,
                          numshow=10) {
  # input=input[[1]]

  if (is.list(input) == FALSE) {
    stop("The input is not a list! Please make sure you are
         using the output from the enrichment analysis")
  }

  if (is.data.frame(input$expt) == FALSE |
      is.data.frame(input$reference) == FALSE |
      is.data.frame(input$shared) == FALSE |
      length(input) != 3 |
      all(names(input) != c("expt", "reference", "shared"))) {
    stop("The input is not a list of three dataframes or
         there are no enriched pathways to plot")
  }

  up <- input$expt
  if (length(up) <= 1) {
    up$Description <- NA
  } else {
    up <- up[up$p.adjust < pvalfilt, ]
    if ( nrow(up) > numshow ) {
      up <- up[order(up$p.adjust)[1:numshow],]
    }
  }
  reference <- input$reference
  if (length(reference) <= 1) {
    reference$Description <- NA
  } else {
    reference <- reference[reference$p.adjust < pvalfilt, ]
    if ( nrow(reference) > numshow) {
      reference <- reference[order(reference$p.adjust)[1:numshow],]
    }
  }
  shared <- input$shared
  if (length(shared) <= 1) {
    shared$Description <- NA
  } else {
    shared <- shared[shared$p.adjust < pvalfilt, ]
    if ( nrow(shared) > numshow ) {
      shared <- shared[order(shared$p.adjust)[1:numshow],]
    }
  }

  pathways <- unique(c(up$Description,
                       reference$Description,
                       shared$Description))
  #print(paste("Pathways", pathways))
  pathways <- pathways[!is.na(pathways)]
  if (is.na(pathways) || length(pathways) == 0) {
    stop("No pathways are significant
         (with adjusted pvalues < user input cutoff)")
  }
  # make a list of all the pathways in up, down, and shared
  heatmapmatrix <- matrix(data = NA,
                          nrow = length(pathways),
                          ncol = 3)
  # make a matrix with as many row as there are pathways
  row.names(heatmapmatrix) <- pathways
  # name the rows with the pathway names

  colnames(heatmapmatrix) <- c("up", "down", "shared")
  # put up, down, and shared as the pathway names

  #print(paste("Dim heatmapmatrix", dim(heatmapmatrix)))

  for (i in 1:length(row.names(heatmapmatrix))) {
    #print(row.names(heatmapmatrix)[i])
    if (row.names(heatmapmatrix)[i] %in% up$Description) {
      num1 <- which(up$Description == row.names(heatmapmatrix)[i])
      heatmapmatrix[i, 1] <- up[num1, 6]
    }

    if (row.names(heatmapmatrix)[i] %in% reference$Description) {
      num2 <- which(reference$Description == row.names(heatmapmatrix)[i])
      heatmapmatrix[i, 2] <- reference[num2, 6]
    }

    if (row.names(heatmapmatrix)[i] %in% shared$Description) {
      num3 <- which(shared$Description == row.names(heatmapmatrix)[i])
      heatmapmatrix[i, 3] <- shared[num3, 6]
    }
  }

  # places the adjusted p-value in the matrix is there is one

  if (removeonlyshared == TRUE) {
    # finds the shared pathways the are not present in up or down
    mycounts <- as.numeric(apply(heatmapmatrix,
                                 1,
                                 function(x) is.na(x[1]) & is.na(x[2])))
    # keeps those that are not only shared
    heatmapinput <- heatmapmatrix[mycounts == 0, ]
  }
  if (removeonlyshared == FALSE) {
    heatmapinput <- heatmapmatrix
  }


  heatmapdata <- as.data.frame(heatmapinput)
  heatmapdata <- heatmapdata[order(heatmapdata$down,
                                   heatmapdata$up,
                                   heatmapdata$shared,
                                   decreasing = TRUE), ]
  # sorts matrix
  heatmapdata$id <- rownames(heatmapdata)
  # makes id
  rownames(heatmapdata) <- c(1:nrow(heatmapdata))


  suppressMessages(meltedheatmapdata <- reshape2::melt(heatmapdata))

  meltedheatmapdata$newid <- stringr::str_wrap(meltedheatmapdata$id, width = 80)

  meltedheatmapdata$id <- factor(meltedheatmapdata$id,
                                 levels = unique(meltedheatmapdata$id))

  theXAxis <- as.character(meltedheatmapdata[,2])
  theYAxis <- meltedheatmapdata[,4]
  #all possible values of X (type) and Y (pathways)

  theUniqueY <- unique(meltedheatmapdata$newid)
  theUniqueX <- c("up", "down", "shared")
  #unique values of X and Y


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
  #Subsitute words with position on the matrix

  dataforHeatmap <- as.data.frame(cbind(as.numeric(theXAxis),
                                      as.numeric(theYAxis),
                                      round(as.numeric(meltedheatmapdata$value)
                                              ,3)))
  formattedHeatmapData <- list_parse2(dataforHeatmap)
  #create final formatting

  hc <- highchart() %>%
    hc_chart(type = "heatmap") %>%
    hc_title(text = title) %>%
    hc_xAxis(categories = theUniqueX) %>%
    hc_yAxis(categories = theUniqueY) %>%
    hc_add_series(name = "matrix location, p-value",
                  data = formattedHeatmapData) %>%
    hc_legend(
      title = "p-value",
      enabled = TRUE
    ) %>%
    hc_exporting(enabled = TRUE)
  p <- hc_colorAxis(hc, minColor = "#000080", maxColor = "#FFFFFF")
  #create final formatting

  return(p)
  }


##############################################################################

#' Given the output from processPathways(), creates a heatmap from
#' the ouput of the GREAT enrichment analysis. Presence or absence of
#' the pathway in enrichment of both type-specific (increased or decreased
#' log2fold change, low p-value) and shared (no change, higher p-value)
#' regulatory regions is plotted.
#'
#' #@param input results from GREAT enrichment analysis
#' #@param title title of the heatmap
#' #@param pathwaycateg ontology
#' #@param numshow number of top pathways (ranked according to p-value) of each type (expt, reference, shared) to show in the plot (default=10)

#' @return heatmap
#'
#' @examples
#' \dontrun{
#' csvfile <- file.path(dir="yourfilepath", 'sampleinfo.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps = 2)
#' plotConsensusPeaks(samplepeaks = consPeaks)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           distancefromTSSdist = 1500,
#'                                           distancefromTSSprox = 1000)
#' counts_consPeaks <- getCounts(annotpeaks = consPeaksAnnotated,
#'                               sampleinfo = sampleinfo,
#'                               reference = 'SAEC',
#'                               chrom = 'chr21')
#' altre_peaks <- countanalysis(counts=counts_consPeaks,
#'                              pval=0.01,
#'                              lfcvalue=1)
#' categaltre_peaks <- categAltrePeaks(altre_peaks,
#'                              lfctypespecific = 1.5,
#'                              lfcshared = 1.2,
#'                              pvaltypespecific = 0.01,
#'                              pvalshared = 0.05)
#' GREAToutput <- runGREAT(peaks = categaltre_peaks)
#' GREATpathways <- processPathways(temp)
#' names(GREATpathways$Sig_Pathways)
# plotGREATenrich(GREATpathways,
#                 title = "GREAT Enrichment Analysis",
#                 pathwaycateg ="GO Molecular Function")
#' }
#'
#' @export
plotGREATenrich <- function(){}
###plotGREATenrich <- function(input,
###                          title,
###                          pathwaycateg,
###                          numshow=10) {
###
###
###  if (is.list(input) == FALSE) {
###    stop("The input is not a list! Please make sure you are
###         using the output from the enrichment analysis (function processPathways())")
###
###  if (is.data.frame(input$expt) == FALSE |
###      is.data.frame(input$reference) == FALSE |
###      is.data.frame(input$shared) == FALSE |
###      length(input) != 3 |
###      all(names(input) != c("expt", "reference", "shared"))) {
###    stop("The input is not a list of three dataframes or
###         there are no enriched pathways to plot")
###  }






#' Plots a venn diagram that compares altered regions as determined by peak presence or by
#' differential counts.  The type of regulatory region (TSS-proximal, TSS-distal, or both)
#' and type of peak comparison (intensity or peak) must be specified.
#' Plots a venn diagram that compares altered regions as determined by peak
#' presence or by differential counts.  The type of regulatory region
#' (TSS-proximal, TSS-distal, or both) and type of peak comparison
#' (intensity or peak) must be specified.
#' @param analysisresultsmatrix analysisresults of Intensity analysis place into
#' analysisresults matrix by the analyzeanalysisresults function
#' @param region pick a region, regions can be 'TSS-distal', 'TSS-proximal',
#' or 'both' -- INCLUDE quotes
#' @param method pick a method, methods can be 'Intensity' or 'Peak'
#' include quotes
#' @param palette RColorBrewer palette to change graph colors
#' @return venn diagram
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
#'plotCompareMethods(comparePeaksAnalysisResults)
#'}

plotCompareMethods <- function(analysisresultsmatrix,
                     region = "both", method = "Intensity", palette = NULL) {


  if ( !is.null(palette) ) {
    cols <- RColorBrewer::brewer.pal(3, palette) }
  else{cols <- c("#00E5EE", "#C71585","#000080")}

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
  string <- paste(rownames(analysisresultsmatrix)[1],
                  rownames(analysisresultsmatrix)[2],
                  rownames(analysisresultsmatrix)[3],
                  rownames(analysisresultsmatrix)[4],
                  rownames(analysisresultsmatrix)[5],
                  rownames(analysisresultsmatrix)[6],
                  rownames(analysisresultsmatrix)[7],
                  rownames(analysisresultsmatrix)[8],
                  rownames(analysisresultsmatrix)[9])

  stringsplit <- strsplit(string, " ")
  uniquestringsplit <- unique(stringsplit[[1]])
  split <- unlist(strsplit(rownames(analysisresultsmatrix)[1], split = " "))
  names <- split[!(split %in% c("TSS-dists"))]
  names <- paste(names, collapse = " ")
  casename <- names

  split <- unlist(strsplit(rownames(analysisresultsmatrix)[4], split = " "))
  names <- split[!(split %in% c("TSS-dists"))]
  names <- paste(names, collapse = " ")
  referencename <- names

  # this is a way to the name of the 'case'
  # from the analysisresults matrix

  p <- highchart() %>%
    hc_chart(type = "pie") %>%
    hc_title(text = paste(region, method),
             style = list(color = '#2E1717',
                          fontWeight = 'bold')) %>%
    hc_plotOptions(
      series = list(showInLegend = TRUE)
    ) %>%
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
    hc_add_series(data = list(
      list(y = case, name = casename, dataLabels = FALSE),
      list(y = reference, name = referencename, dataLabels = FALSE),
      list(y = shared, name = "Shared", dataLabels = FALSE)
    ), name = paste(region, method)) %>%
    hc_colors(cols)  %>%
    hc_exporting(enabled = TRUE)
  return(p)
}

#' Plots venn diagrams for comparison of two methods of identifying altered
#' regulatory regions Makes venn diagrams for TSS-proximal, TSS-distal, and
#' combined for both intensity-based peaks and for peaks identified by hotspot
#' calling algorithms.  There is no return value. Six venn diagrams will be
#' plotted
#' @param analysisresultsmatrix analysisresults of countanalysis function
#' place into a a analysisresults matrix by the analyzeanalysisresults function
#' @param viewer whether the plot should be displayed in the RStudio viewer or
#' in Shiny/Knittr
#' @param palette RColorBrewer palette to change graph colors
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
#'plotCompareMethodsAll(comparePeaksAnalysisResults)
#' }
#' @export
#'
plotCompareMethodsAll <- function(analysisresultsmatrix, viewer = TRUE,
                                  palette = NULL) {


  if ( !is.null(palette) ) {
    cols <- RColorBrewer::brewer.pal(3, palette) }
  else{cols <- c("#00E5EE", "#C71585","#000080")}

  analysisresultsmatrix <- analysisresultsmatrix[[1]]

  if (is.matrix(analysisresultsmatrix) ==
      FALSE) {
    stop("The input is not a matrix!")
  }

  p1 <- plotCompareMethods(analysisresultsmatrix,
                    "TSS-proximal", "Intensity", palette = palette)
  p2 <- plotCompareMethods(analysisresultsmatrix,
                    "TSS-distal", "Intensity", palette = palette)
  p3 <- plotCompareMethods(analysisresultsmatrix,
                    "both", "Intensity", palette = palette)
  p4 <- plotCompareMethods(analysisresultsmatrix,
                    "TSS-proximal", "Peak", palette = palette)
  p5 <- plotCompareMethods(analysisresultsmatrix,
                    "TSS-distal", "Peak", palette = palette)
  p6 <- plotCompareMethods(analysisresultsmatrix,
                    "both", "Peak", palette = palette)

  if (viewer == TRUE) {
    p <- htmltools::browsable(hw_grid(p1, p2, p3, p4, p5, p6, ncol = 3,
                                      rowheight = 300))
  }
  else {
    p <- hw_grid(p1, p2, p3, p4, p5, p6, ncol = 3)
  }
  return(p)

return(p)
}

