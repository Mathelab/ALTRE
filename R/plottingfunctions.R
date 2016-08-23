#' Given the output from getConsensusPeaks, generate a barplot
#' of countstatistics
#'
#' @param samplepeaks output generated from getConsensusPeaks
#'
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' csvfile <- file.path(dir="yourfilepath", 'sampleinfo.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks, minreps=2)
#' plotConsensusPeaks(samplepeaks = consPeaks)
#' }
#' @export
#'
plotConsensusPeaks <- function(samplepeaks) {

    CellType <- NULL
    #without this R CMD check throws no visible binding for global variable error

    consPeaksStats <- samplepeaks$consPeaksStats
    row.names(consPeaksStats) <- NULL
    consPeaksStats[ , 2] <-  as.numeric(as.character(consPeaksStats[[2]]))
    consPeaksStats[ , 3] <-  as.numeric(as.character(consPeaksStats[[3]]))
    statsFormated <- tidyr::gather(consPeaksStats, "CellType", "Count", 2:3)

    plottingData <- statsFormated %>% split(levels(as.factor(statsFormated$CellType)))

    p <- highchart() %>%
        hc_title(text = "Peak Counts by Cell Type",
                 style = list(color = '#2E1717',
                              fontWeight = 'bold')) %>%
        hc_subtitle(text = "For bioreplicates and their merged consensus track") %>%
        hc_add_series(
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
            layout = "vertical",
            align = "right",
            verticalAlign = "top",
            floating = TRUE,
            x = -5,
            y = 60
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
#' @param viewer whether the plot should be displayed in the RStudio viewer or
#'        in Shiny/Knittr
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' csvfile <- file.path(dir="yourfilepath", 'sampleinfo.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           distancefromTSSdist = 1500,
#'                                           distancefromTSSprox = 1000)
#' plotCombineAnnotatePeaks(consPeaksAnnotated)
#' }
#' @export
#'
plotCombineAnnotatePeaks <- function(conspeaks, viewer = TRUE) {

    CellType <- NULL
    #without this R CMD check throws no visible binding for global variable error

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
        mergeStatsTotal <- dplyr::filter(mergeStatsFormatted, CellType == "TotalNumber")
        thecondition <- matrix(unlist(strsplit(mergeStatsTotal$Condition, "_")), nrow = 3, ncol = 4)[2,]
        mergeStatsBefore <- dplyr::filter(mergeStatsTotal, thecondition == "before")
        mergeStatsAfter <- dplyr::filter(mergeStatsTotal, thecondition == "after")

        p1 <- highchart() %>%
            hc_title(text = "Number of REs before/after merging",
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
            hc_xAxis(categories = c("beforemerging", "aftermerging")) %>%
            hc_legend(
                enabled = TRUE,
                layout = "vertical",
                align = "right",
                verticalAlign = "top",
                floating = TRUE,
                x = -5,
                y = 60
            ) %>%
            hc_tooltip(
                headerFormat = "<b>{series.name}_{point.key}</b><br>",
                pointFormat = "{point.y}",
                valueSuffix = ' peaks'
            ) %>%
            hc_exporting(enabled = TRUE)
    }

    feature <- "MeanLength"
    if ( feature == "MeanLength" ) {
      mergeStatsMean <- dplyr::filter(mergeStatsFormatted, CellType == "MeanLength")
      thecondition <- matrix(unlist(strsplit(mergeStatsMean$Condition, "_")),
                             nrow = 3,
                             ncol = 4)[2, ]
      mergeStatsBefore <- dplyr::filter(mergeStatsMean, thecondition == "before")
      mergeStatsAfter <- dplyr::filter(mergeStatsMean, thecondition == "after")


        p2 <- highchart() %>%
            hc_title(text = "Mean length of REs before/after merging",
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
            hc_xAxis(categories = c("beforemerging", "aftermerging")) %>%
            hc_legend(
                enabled = TRUE,
                layout = "vertical",
                align = "right",
                verticalAlign = "top",
                floating = TRUE,
                x = -5,
                y = 60
            ) %>%
            hc_tooltip(
                headerFormat = "<b>{series.name}_{point.key}</b><br>",
                pointFormat = "{point.y}",
                valueSuffix = ' peaks'
            ) %>%
            hc_exporting(enabled = TRUE)
    }


    if (viewer == TRUE) {
    p <- htmltools::browsable(hw_grid(p1, p2, ncol = 1, rowheight = 300))
    }
    else {
    p <- hw_grid(p1, p2)
    }
    return(p)
}


#' Given the output from getCounts, plot a density plot
#'  of log2 RPKM values of regulation regions
#'
#' @param altrepeakscateg output generated from countanalysis() then
#' categAltrePeaks()
#'
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' csvfile <- file.path(dir="yourfilepath", 'sampleinfo.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot<- getTSS()
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
#' altre_peaks <- countanalysis(counts = counts_consPeaks,
#'                              pval = 0.01,
#'                              lfcvalue = 1)
#' categaltre_peaks <- categAltrePeaks(altre_peaks,
#'                                     lfctypespecific = 1.5,
#'                                     lfcshared = 1.2,
#'                                     pvaltypespecific = 0.01,
#'                                     pvalshared = 0.05)
#' plotCountAnalysis(categaltre_peaks)
#' }
#' @export

plotCountAnalysis <- function(altrepeakscateg) {

  toplot <- altrepeakscateg$analysisresults[ ,c("region",
                                                "log2FoldChange",
                                                "padj",
                                                "REaltrecateg")]
  tssdist <- toplot[which(toplot$region == "TSS-distal"), ]
  tssprox <- toplot[which(toplot$region == "TSS-proximal"), ]
  lengthRE <- rep("", length(tssdist$REaltrecateg))

  num1 <- min(which(tssdist$REaltrecateg == "Experiment Specific"))
  num2 <- min(which(tssdist$REaltrecateg == "Reference Specific"))
  num3 <- min(which(tssdist$REaltrecateg == "Shared"))
  num4 <- min(which(tssdist$REaltrecateg == "Ambiguous"))

  lengthRE[num1] <- "Experiment Specific"
  lengthRE[num2] <- "Reference Specific"
  lengthRE[num3] <- "Shared"
  lengthRE[num4] <- "Ambiguous"

  p1 <- highchart() %>%
    hc_title(text = "TSS-distal",
             style = list(color = '#2E1717',
                          fontWeight = 'bold')) %>%
    hc_add_series_scatter(y = -log10(tssdist$padj),
                          x = tssdist$log2FoldChange,
                          color = tssdist$REaltrecateg,
                          label = lengthRE)  %>%
    hc_tooltip(headerFormat = "",
               pointFormat  = "<b>log2FC</b> = {point.x}<br> <b>-log10pvalue</b> = {point.y}<br>") %>%
    hc_exporting(enabled = TRUE)

  lengthRE <- rep("", length(tssprox$REaltrecateg))

  num1 <- min(which(tssprox$REaltrecateg == "Experiment Specific"))
  num2 <- min(which(tssprox$REaltrecateg == "Reference Specific"))
  num3 <- min(which(tssprox$REaltrecateg == "Shared"))
  num4 <- min(which(tssprox$REaltrecateg == "Ambiguous"))

  lengthRE[num1] <- "Experiment Specific"
  lengthRE[num2] <- "Reference Specific"
  lengthRE[num3] <- "Shared"
  lengthRE[num4] <- "Ambiguous"


  p2 <- highchart() %>%
    hc_title(text = "TSS-proximal",
             style = list(color = '#2E1717',
                          fontWeight = 'bold')) %>%
    hc_add_series_scatter(y = -log10(tssprox$padj),
                          x = tssprox$log2FoldChange,
                          color = tssprox$REaltrecateg,
                          label = lengthRE)  %>%
    hc_tooltip(headerFormat = "",
               pointFormat  = "<b>log2FC</b> = {point.x}<br> <b>-log10pvalue</b> = {point.y}<br>") %>%
    hc_exporting(enabled = TRUE)

  plot <- htmltools::browsable(hw_grid(p1, p2, ncol = 2, rowheight = 550))
  htmlcode <- hw_grid(p1, p2)

  return(plot)
}


###############################################################################
#' Creates a boxplot to see the distribution of read counts in type-specific and
#' shared TSS-proximal and TSS-distal regions.
#'
#' Takes the rlog transformation of the RRKM (Reads Per Kilobase of transcript
#' per Million) of the read counts of type-specific and shared regulatory regions
#' and plots the distribution of those read counts in all sample types analyzed
#' in the workflow.
#'
#' @param analysisresults output generated from countanalysis() then categAltrePeaks()
#' @param counts output generated from getCounts()
#'
#' @return a highcharter object
#'
#' @examples
#' \dontrun{
#' csvfile <- file.path(dir="yourfilepath", 'sampleinfo.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot<- getTSS()
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
#' altre_peaks <- countanalysis(counts = counts_consPeaks,
#'                              pval = 0.01,
#'                              lfcvalue = 1)
#' categaltre_peaks <- categAltrePeaks(altre_peaks,
#'                                     lfctypespecific = 1.5,
#'                                     lfcshared = 1.2,
#'                                     pvaltypespecific = 0.01,
#'                                     pvalshared = 0.05)
#' plotDistCountAnalysis(categaltre_peaks, counts_consPeaks)
#' }
#' @export
#'
plotDistCountAnalysis <- function(analysisresults, counts) {

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

  Experimentspecific_list <- list(distal1_5num_A549, proximal1_5num_A549, distal1_5num_SAEC, proximal1_5num_SAEC)
  Ambiguous_list <- list(distal2_5num_A549, proximal2_5num_A549, distal2_5num_SAEC, proximal2_5num_SAEC)
  Shared_list <- list(distal3_5num_A549, proximal3_5num_A549, distal3_5num_SAEC, proximal3_5num_SAEC)
  Referencespecific_list <- list(distal4_5num_A549, proximal4_5num_A549, distal4_5num_SAEC, proximal4_5num_SAEC)

  categ <- c('A549-specific TSS-distal', 'A549-specific TSS-proximal', 'SAEC-specific TSS-distal', 'SAEC-specific TSS-proximal')
  p <- highchart() %>%
    hc_title(text = "Distribution of Normalized Counts",
             style = list(color = '#2E1717',
                          fontWeight = 'bold')) %>%
    hc_plotOptions(
      boxplot = list(
        fillColor = '#F0F0E0',
        lineWidth = 2,
        medianColor = '#0C5DA5',
        medianWidth = 3,
        stemColor = '#A63400',
        stemDashStyle = 'dot',
        stemWidth = 1,
        whiskerColor = '#3D9200',
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
    hc_exporting(enabled = TRUE)
    return(p)
}

###############################################################################


#' Given the output from getCounts, plot a density plot of
#' log2 RPKM values of regulation regions
#'
#' @param countsconspeaks output generated from getCounts
#'
#' @return a ggplot
#'
#' @examples
#' \dontrun{
#' csvfile <- file.path(dir="yourfilepath", 'sampleinfo.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps=2)
#' plotConsensusPeaks(samplepeaks = consPeaks)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           distancefromTSSdist = 1500,
#'                                           distancefromTSSprox = 1000 )
#'
#' counts_consPeaks <- getCounts(annotpeaks = consPeaksAnnotated,
#'                              sampleinfo = sampleinfo,
#'                              reference = 'SAEC',
#'                              chrom = 'chr21')
#' plotgetcounts(counts_consPeaks)
#' }
#' @export

plotgetcounts <- function(countsconspeaks) {
  mydf <- countsconspeaks$regioncountsforplot
  varstack <- suppressMessages(reshape2::melt(mydf))
  varstack$concat <- paste(varstack$region,
                           varstack$variable,
                           sep = ": ")
  varstack$concat <- sub("librarysize.*", "", varstack$concat)
  densityplot <- ggplot(varstack, aes(x = varstack$value)) +
    geom_density(aes(group = varstack$concat,
                     color = varstack$concat,
                     fill = varstack$concat),
                 alpha = 0.3) +
    theme_bw(base_size = 15, base_family = "") +
    theme(legend.title = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "log2 read counts \n(normalized by library and region sizes)")

  return(densityplot)
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
#' @param numshow number of top pathways (ranked according to p-value) of each type (expt, reference, shared) to show in the plot (default=10)
#'
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
#'				lfctypespecific = 1.5,
#'				lfcshared = 1.2,
#'				pvaltypespecific = 0.01,
#'				pvalshared = 0.05)
#' MFenrich <- pathenrich(analysisresults = categaltre_peaks,
#'                        ontoltype = 'MF',
#'                        enrichpvalfilt = 0.01)
#' BPenrich <- pathenrich(analysisresults=categaltre_peaks,
#'                        ontoltype='BP',
#'                        enrichpvalfilt=0.01)
#' plot1 <- enrichHeatmap(MFenrich, title='GO:MF, p<0.01')
#' plot2 <- enrichHeatmap(BPenrich, title='GO:BP, p<0.01')
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

  dataforHeatmap <- as.data.frame(cbind(as.numeric(theXAxis),as.numeric(theYAxis),as.numeric(meltedheatmapdata$value)))
  formattedHeatmapData <- list_parse2(dataforHeatmap)
  #create final formatting

  hc <- highchart() %>%
    hc_chart(type = "heatmap") %>%
    hc_title(text = title) %>%
    hc_xAxis(categories = theUniqueX) %>%
    hc_yAxis(categories = theUniqueY) %>%
    hc_add_series(name = "matrix location, p-value", data = formattedHeatmapData)
  p <- hc_colorAxis(hc, minColor = "#FFFFFF", maxColor = "#000080")
  #create final formatting

  return(p)
  }


#' Multiple plot function
#'
#' Plots multiple ggplot objects in one window.
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' This function was not written by the authours of this package. Link:
#' http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#'
#' @param ... list of plots
#' @param cols number of columns in layout
#' @param layout matrix specifying layout, if present cols is ignored
#' @param plotlist plotlist
#' @param file file
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
#'                                           distancefromTSSprox = 1000 )
#'counts_consPeaks <- getCounts(annotpeaks = consPeaksAnnotated,
#'                              sampleinfo = sampleinfo,
#'                              reference = 'SAEC',
#'                              chrom = 'chr21')
#' altre_peaks <- countanalysis(counts = counts_consPeaks,
#'                              pval = 0.01,
#'                              lfcvalue = 1)
#' categaltre_peaks <- categAltrePeaks(altre_peaks,
#'                                     lfctypespecific = 1.5,
#'                                     lfcshared = 1.2,
#'                                     pvaltypespecific = 0.01,
#'                                     pvalshared = 0.05)
#'MFenrich <- pathenrich(analysisresults = categaltre_peaks,
#'                       ontoltype = 'MF',
#'                       enrichpvalfilt = 0.01)
#'BPenrich <- pathenrich(analysisresults= categaltre_peaks,
#'                       ontoltype='BP',
#'                       enrichpvalfilt=0.01)
#' plot1 <- enrichHeatmap(MFenrich, title='GO:MF, p<0.01')
#' plot2 <- enrichHeatmap(BPenrich, title='GO:BP, p<0.01')
#' multiplot(plot1,plot2,cols=1)
#' }
#'
#' @export
#'
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots <- length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel ncol: Number of columns of plots nrow: Number of rows
    # needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols,
                     nrow = ceiling(numPlots/cols))
  }

  if (numPlots == 1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' Plots a venn diagram that compares altered regions as determined by peak presence or by
#' differential counts.  The type of regulatory region (TSS-proximal, TSS-distal, or both)
#' and type of peak comparison (intensity or peak) must be specified.
#' @param analysisresultsmatrix analysisresults of Intensity analysis place into
#' analysisresults matrix by the analyzeanalysisresults function
#' @param region pick a region, regions can be 'TSS-distal', 'TSS-proximal', or 'both'
#' INCLUDE quotes
#' @param method pick a method, methods can be 'intensity' or 'peak'
#' include quotes
#' @param color include the colors you want in your venn diagram
#' @return venn diagram
#' @examples
#' \dontrun{
#' csvfile <- file.path(dir="yourfilepath", 'sampleinfo.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps = 2)
#' plotConsensusPeaks(samplepeaks = consPeaks)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                            TSS = TSSannot,
#'                                            merge = TRUE,
#'                                            regionspecific = TRUE,
#'                                            distancefromTSSdist = 1500,
#'                                            distancefromTSSprox = 1000)
#' counts_consPeaks <- getCounts(annotpeaks = consPeaksAnnotated,
#'                               sampleinfo = sampleinfo,
#'                               reference = 'SAEC')
#' altre_peaks <- countanalysis(counts=counts_consPeaks,
#'                              pval=0.01,
#'                              lfcvalue=1)
#' categaltre_peaks=categAltrePeaks(altre_peaks,
#'	lfctypespecific=1.5,
#' 	lfcshared=1.2,
#' 	pvaltypespecific=0.01,
#' 	pvalshared=0.05)
#' analysisresults <- comparePeaksAltre(categaltre_peaks, reference= "SAEC")
#' plot1 <- plotvenn(analysisresults,
#'                   region='TSS-distal',
#'                   method='intensity')
#'}

plotvenn <- function(analysisresultsmatrix,
                     region = "both", method = "intensity",
                     color = "redorange") {

  if (region == "TSS-proximal") {
    feature <- c("TSS-proxs")
    coordinates <- c(2, 5, 8)
  }
  if (region == "TSS-distal") {
    feature <- c("TSS-dists")
    coordinates <- c(1, 4, 7)
  }
  if (region == "both") {
    region <- c("TSS-distal/TSS-proximal")
    feature <- c("TSS-dists", "TSS-proxs")
    coordinates <- c(3, 6, 9)
  }
  # identifies the correct numbers from the
  # analysisresults matrix based on the
  # regulatory region of interest
  if (method == "intensity") {
    case <- analysisresultsmatrix[coordinates[1], 1]
    reference <- analysisresultsmatrix[coordinates[2], 1]
    shared <- analysisresultsmatrix[coordinates[3], 1]
  }

  if (method == "peak") {
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
  names <- split[!(split %in% c("TSS-distal"))]
  names <- paste(names, collapse = " ")
  casename <- names

  split <- unlist(strsplit(rownames(analysisresultsmatrix)[4], split = " "))
  names <- split[!(split %in% c("TSS-distal"))]
  names <- paste(names, collapse = " ")
  referencename <- names

  # this is a way to the name of the 'case'
  # from the analysisresults matrix


  p <- highchart() %>%
    hc_chart(type = "pie") %>%
    hc_title(text = paste(method, region),
             style = list(color = '#2E1717',
                          fontWeight = 'bold')) %>%
    hc_plotOptions(
      series = list(showInLegend = TRUE)
    ) %>%
    hc_legend(
      enabled = FALSE,
      layout = "vertical",
      align = "right",
      verticalAlign = "top",
      floating = TRUE,
      x = -5,
      y = 60
    ) %>%

    hc_add_series(data = list(
      list(y = case, name = casename),
      list(y = reference, name = referencename),
      list(y = shared, name = "Shared")
    )
    )

  return(p)
}

#' Plots venn diagrams for comparison of two methods of identifying altered
#' regulatory regions Makes venn diagrams for TSS-proximal, TSS-distal, and combined
#' for both intensity-based peaks and for peaks identified by hotspot calling
#' algorithms.  There is no return value. Six venn diagrams will be plotted
#' @param analysisresultsmatrix analysisresults of countanalysis function
#' place into a a analysisresults matrix by the analyzeanalysisresults function
#' @examples
#' \dontrun{
#' csvfile <- file.path(dir="yourfilepath", 'sampleinfo.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps = 2)
#' plotConsensusPeaks(samplepeaks = consPeaks)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                            TSS = TSSannot,
#'                                            merge = TRUE,
#'                                            regionspecific = TRUE,
#'                                            distancefromTSSdist = 1500,
#'                                            distancefromTSSprox = 1000 )
#' counts_consPeaks <- getCounts(annotpeaks = consPeaksAnnotated,
#'                               sampleinfo = sampleinfo,
#'                               reference = 'SAEC')
#' altre_peaks <- countanalysis(counts=counts_consPeaks,
#'                              pval=0.01,
#'                              lfcvalue=1)
#' categaltre_peaks=categAltrePeaks(altre_peaks,
#'      lfctypespecific=1.5,
#'      lfcshared=1.2,
#'      pvaltypespecific=0.01,
#'      pvalshared=0.05)
#' analysisresults <- resultsComparison(altre_peaks, reference= "SAEC")
#' plotallvenn(analysisresults)
#' }
#' @export


plotallvenn <- function(analysisresultsmatrix) {
analysisresultsmatrix <- analysisresultsmatrix[[1]]

  if (is.matrix(analysisresultsmatrix) ==
      FALSE) {
    stop("The input is not a matrix!")
  }

  p1 <- plotvenn(analysisresultsmatrix,
                    "TSS-proximal", "intensity", "bluegreen")
  p2 <- plotvenn(analysisresultsmatrix,
                    "TSS-distal", "intensity", "bluegreen")
  p3 <- plotvenn(analysisresultsmatrix,
                    "both", "intensity", "bluegreen")
  p4 <- plotvenn(analysisresultsmatrix,
                    "TSS-proximal", "peak", "redorange")
  p5 <- plotvenn(analysisresultsmatrix,
                    "TSS-distal", "peak", "redorange")
  p6 <- plotvenn(analysisresultsmatrix,
                    "both", "peak", "redorange")

  plot <- htmltools::browsable(hw_grid(p1, p2, p3, p4, p5, p6, ncol = 3, rowheight = 300))

return(plot)
}

