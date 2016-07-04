#' Given the output from getConsensusPeaks, generate a barplot
#' of countstatistics
#'
#' @param samplepeaks output generated from getConsensusPeaks
#'
#' @return a ggplot
#'
#' @examples
#' \dontrun{
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks, minreps=2)
#' plotConsensusPeaks(samplepeaks = consPeaks)
#' }
#' @export
#'
plotConsensusPeaks <- function(samplepeaks) {
  dfstats <- samplepeaks$consPeaksStats
  # quick fix: change to numeric
  row.names(dfstats) <- NULL
  dfstats[ , 2] <-  as.numeric(as.character(dfstats[[2]]))
  dfstats[ , 3] <-  as.numeric(as.character(dfstats[[3]]))
  ##
  mydf <- tidyr::gather(dfstats, "CellType", "Count", 2:3)

  p <- ggplot(data = mydf, aes_string(x = "CellType",
                                      y = "Count",
                                      fill = "PeakType")) +
    geom_bar(stat = "identity",
             position = position_dodge()) +
    geom_text(aes_string(label = "Count",
                         x = "CellType"),
              position = position_dodge(width = 1),
              size = 3,
              hjust = 0.5,
              vjust = -.3) +
    theme_bw(base_size = 12, base_family = "Helvetica") +
    theme(aspect.ratio = 0.8,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "Sample Type",
         y = "Number of Peaks") +
    ggtitle("Number of Peaks for Bioreplicates  \n and their Merged Consensus") +
    scale_fill_brewer(palette = "Set2")
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
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           mergedistenh = 1500,
#'                                           mergedistprom = 1000)
#' plotCombineAnnotatePeaks(consPeaksAnnotated)
#' }
#' @export
#'
plotCombineAnnotatePeaks <- function(conspeaks) {

  # quick fix
  mydf <- conspeaks$mergestats[ , c(2, 3)]
  row.names(mydf) <- conspeaks$mergestats[ , 1]
  ##

  if (nrow(mydf) == 1) {
    stop("No plot to show since merging was not performed
         when calling combineAnnotatePeaks function")
  } else {
    numreg <- as.data.frame(mydf$TotalNumber)
    colnames(numreg) <- "TotalNumber"
    numreg$REtype <- gsub("_.*", "", rownames(mydf))
    numreg$REmerge <- rownames(mydf)
    numreg$REmerge <- gsub("enhancers_", "", numreg$REmerge)
    numreg$REmerge <- gsub("promoters_", "", numreg$REmerge)
    numreg$REmerge <- factor(numreg$REmerge,
                             levels = c("before_merging", "after_merging"))

    plot1 <- ggplot(numreg, aes_string(x = "REtype",
                                       y = "TotalNumber",
                                       fill = "REmerge")) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      geom_text(aes_string(label = "TotalNumber",
                           x = "REtype",
                           y = "TotalNumber",
                           ymax = "TotalNumber"),
                position = position_dodge(width = 1),
                size = 3,
                hjust = 0.5,
                vjust = -.5) +
      scale_fill_brewer(palette = "Set2") +
      theme_bw(base_size = 12) +
      theme(aspect.ratio = 1.5,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      labs(x = "", y = "Number of Regulatory Regions") +
      ggtitle("Number of REs\nBefore/After merging")


    meanlength <- as.data.frame(mydf$MeanLength)
    colnames(meanlength) <- "MeanLength"
    meanlength$REtype <- gsub("_.*", "", rownames(mydf))
    meanlength$REmerge <- rownames(mydf)
    meanlength$REmerge <- gsub("enhancers_", "", meanlength$REmerge)
    meanlength$REmerge <- gsub("promoters_", "", meanlength$REmerge)
    meanlength$REmerge <- factor(meanlength$REmerge,
                                 levels = c("before_merging", "after_merging"))

    plot2 <- ggplot(meanlength, aes_string(x = "REtype",
                                           y = "MeanLength",
                                           fill = "REmerge")) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      geom_text(aes_string(label = "MeanLength",
                           x = "REtype",
                           y = "MeanLength",
                           ymax = "MeanLength"),
                position = position_dodge(width = 1),
                size = 3,
                hjust = 0.5,
                vjust = -.5) +
      scale_fill_brewer(palette = "Set2") +
      theme_bw(base_size = 12) +
      theme(aspect.ratio = 1.5,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      labs(x = "", y = "Mean length of Regulatory Regions") +
      ggtitle("Mean length of REs\nBefore/After merging")

    multiplot(plot1, plot2, cols = 2)
  }
}


#' Given the output from getcounts, plot a density plot
#'  of log2 RPKM values of regulation regions
#'
#' @param altrepeakscateg output generated from countanalysis() then categAltrePeaks()
#'
#' @return a ggplot
#'
#' @examples
#' \dontrun{
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot<- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           mergedistenh = 1500,
#'                                           mergedistprom = 1000)
#' counts_consPeaks <- getcounts(annotpeaks = consPeaksAnnotated,
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
  enh <- toplot[which(toplot$region == "enhancer"), ]
  prom <- toplot[which(toplot$region == "promoter"), ]

  plot1 <- ggplot(enh,
                 aes(enh$log2FoldChange,
                     -log2(enh$padj))) +
    geom_point(aes(col = factor(enh$REaltrecateg))) +
    scale_colour_manual(values = c("dark grey","salmon","dark green","blue")) +
    theme_bw(base_size = 15) +
    theme(legend.title = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "log2FC", y = "-log2(pvalue)") +
    ggtitle("Enhancers")

  plot2 <- ggplot(prom,
                 aes(prom$log2FoldChange,
                     -log2(prom$padj))) +
    geom_point(aes(col = factor(prom$REaltrecateg))) +
      scale_fill_brewer(palette = "Set2") +
    theme_bw(base_size = 15) +
    theme(legend.title = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "log2FC", y = "-log2(pvalue)") +
    ggtitle("Promoters")

#    geom_hline(aes(yintercept = -log2(pval)), linetype = "dashed") +
#    geom_vline(aes(xintercept = (-lfcvalue)), linetype = "dashed") +
#    geom_vline(aes(xintercept =  (lfcvalue)), linetype = "dashed")

  return(multiplot(plot1,plot2))
}

###############################################################################
#' Creates a boxplot to see the distribution of read counts in type-specific and shared enhancers
#'
#' Takes the rlog transformation of the RRKM (Reads Per Kilobase of transcript per Million)
#' of the read counts of type-specific and shared regulatory regions and plots
#' the distribution of those read counts in all sample types analyzed in the workflow.
#'
#' @param analysisresults output generated from countanalysis() then categAltrePeaks()
#' @param counts output generated from getcounts()
#'
#' @return a ggplot
#'
#' @examples
#' \dontrun{
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot<- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           mergedistenh = 1500,
#'                                           mergedistprom = 1000)
#' counts_consPeaks <- getcounts(annotpeaks = consPeaksAnnotated,
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
	readcounts=counts$regioncounts
	analysisresults=analysisresults[[1]]
	errortest=try(SummarizedExperiment::assay(readcounts), silent=TRUE)
	if (inherits(errortest, 'try-error')==TRUE){
	  stop("The input for the readcounts arguement is not a summerized experiment object!")
	}

	if (is.data.frame(analysisresults)==FALSE)
	{stop("The input for the analysisresults arguement is not a dataframe!")

	}

	# Check that counts and analysisresults are in the same order
	countsinfo=as.data.frame(SummarizedExperiment::rowRanges(readcounts))
	countcoord=paste0(countsinfo$seqnames,countsinfo$start,countsinfo$end)
	analcoord=paste0(analysisresults$chr,analysisresults$start,analysisresults$stop)

	if(!all.equal(analcoord,countcoord)) {
	    stop("The peaks in the analysisresults and counts are not the same")
	}

	PEcateg=analysisresults$region
	altrecateg=analysisresults$REaltrecateg
	samplecateg=

	# Get log2FPM values:
	log2FPM=log2(DESeq2::fpkm(readcounts,robust=TRUE)+0.001)
	# Average log2FPM values over replicats:
	sampletypes= SummarizedExperiment::colData(readcounts)$sample
	meanlog2FPM=c()
	for (i in unique(sampletypes)) {
		samp=which(sampletypes==i)
		meanlog2FPM=cbind(meanlog2FPM,
			as.numeric(apply(log2FPM[,samp],1,mean)))
	}
	colnames(meanlog2FPM)=unique(sampletypes)

	mydf=data.frame(meanlog2FPM=meanlog2FPM,
	    PEcateg=PEcateg,
	    altrecateg=altrecateg)

        meltdf=reshape2::melt(mydf)
	meltdf$variable=gsub("meanlog2FPM.","",meltdf$variable)

	boxplot=ggplot(meltdf, aes_string(x="PEcateg", y="value")) +
             geom_boxplot(aes_string(fill="altrecateg"),position = position_dodge(width = .8)) +
	     facet_grid(.~variable) +
	     scale_fill_manual(values=c("grey","salmon","darkgreen","blue")) +
             theme_bw() + ggtitle("Distribution of Normalized Counts") +
             xlab("") + ylab("log2(FPKM)") +
        	theme(axis.line = element_line(colour = "black"),
                        axis.title=element_text(size=12,face="bold"),
                        plot.title=element_text(size=14,face="bold"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        legend.key=element_blank())
    return(boxplot)
}

###############################################################################


#' Given the output from getcounts, plot a density plot of
#' log2 RPKM values of regulation regions
#'
#' @param countsconspeaks output generated from getcounts
#'
#' @return a ggplot
#'
#' @examples
#' \dontrun{
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps=2)
#' plotConsensusPeaks(samplepeaks = consPeaks)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           mergedistenh = 1500,
#'                                           mergedistprom = 1000 )
#'
#' counts_consPeaks <- getcounts(annotpeaks = consPeaksAnnotated,
#'                              csvfile = csvfile,
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
#'
#' @return heatmap
#'
#' @examples
#' \dontrun{
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps = 2)
#' plotConsensusPeaks(samplepeaks = consPeaks)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           mergedistenh = 1500,
#'                                           mergedistprom = 1000)
#' counts_consPeaks <- getcounts(annotpeaks = consPeaksAnnotated,
#'                               csvfile = csvfile,
#'                               reference = 'SAEC',
#'                               chrom = 'chr21')
#' altre_peaks <- countanalysis(counts=counts_consPeaks,
#'                              pval=0.01,
#'                              lfcvalue=1)
#' MFenrich <- pathenrich(analysisresults = altre_peaks,
#'                        ontoltype = 'MF',
#'                        enrichpvalfilt = 0.01)
#' BPenrich <- pathenrich(analysisresults=altre_peaks,
#'                        ontoltype='BP',
#'                        enrichpvalfilt=0.01)
#' plot1 <- enrichHeatmap(MFenrich, title='GO:MF, p<0.01')
#' plot2 <- enrichHeatmap(BPenrich, title='GO:BP, p<0.01')
#' }
#' @export

enrichHeatmap <- function(input,
                          title,
                          pvalfilt = 0.01,
                          removeonlyshared = FALSE) {
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
  }
  reference <- input$reference
  if (length(reference) <= 1) {
    reference$Description <- NA
  } else {
    reference <- reference[reference$p.adjust < pvalfilt, ]
  }
  shared <- input$shared
  if (length(shared) <= 1) {
    shared$Description <- NA
  } else {
    shared <- shared[shared$p.adjust < pvalfilt, ]
  }

  pathways <- unique(c(up$Description,
                       reference$Description,
                       shared$Description))
  print(paste("Pathways", pathways))
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

  print(paste("Dim heatmapmatrix", dim(heatmapmatrix)))

  for (i in 1:length(row.names(heatmapmatrix))) {
    print(row.names(heatmapmatrix)[i])
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
    mycounts <- as.numeric(apply(heatmapmatrix,
                                 1,
                                 function(x) is.na(x[1]) & is.na(x[2])))
    # finds the shared pathways the are not present in up or down
    heatmapinput <- heatmapmatrix[mycounts == 0, ]
    # keeps those that are not only shared
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
  meltedheatmapdata <- reshape2::melt(heatmapdata)

  meltedheatmapdata$id <- factor(meltedheatmapdata$id,
                                 levels = unique(meltedheatmapdata$id))

  p1 <- ggplot(meltedheatmapdata,
               aes(y = meltedheatmapdata$id,
                   x = meltedheatmapdata$variable)) +
    geom_tile(aes(fill = meltedheatmapdata$value),
              colour = "black") +
    scale_fill_continuous(low = "navy",
                          high = "turquoise1",
                          na.value = "white",
                          guide = guide_legend(title = "Pvalue")) +
    theme(text = element_text(size = 13)) +
    theme(axis.text.x = element_text(angle=45,hjust=1)) +
    theme(strip.text.y = element_text(colour="red")) +
    ggtitle(title)
  p2 <- p1 +
    scale_x_discrete(expand = c(0, 0),
                     labels = c("High in Expt","High in Reference", "Shared")) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())

  return(p2)
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
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps = 2)
#' plotConsensusPeaks(samplepeaks = consPeaks)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           mergedistenh = 1500,
#'                                           mergedistprom = 1000 )
#'counts_consPeaks <- getcounts(annotpeaks = consPeaksAnnotated,
#'                              csvfile = csvfile,
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
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' Plots a venn digram to compare sample types Plots the number of type-specific
#' and shared regulatory regions in a venn diagram Type of regulatroy region
#' (enhancer, promoter, or both) and type  of peak comparison (intensity or peak)
#' must be specified.
#' @param analysisresultsmatrix analysisresults of Intensity analysis place into
#' analysisresults matrix by the analyzeanalysisresults function
#' @param region pick a region, regions can be 'enhancer', 'promoter', or 'both'
#' INCLUDE quotes
#' @param method pick a method, methods can be 'intensity' or 'peak'
#' include quotes
#' @param color include the colors you want in your venn diagram
#' @return venn diagram
#' @examples
#' \dontrun{
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps = 2)
#' plotConsensusPeaks(samplepeaks = consPeaks)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                            TSS = TSSannot,
#'                                            merge = TRUE,
#'                                            regionspecific = TRUE,
#'                                            mergedistenh = 1500,
#'                                            mergedistprom = 1000)
#' counts_consPeaks <- getcounts(annotpeaks = consPeaksAnnotated,
#'                               csvfile = csvfile,
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
#'                   region='enhancer',
#'                   method='intensity')
#'}

plotvenn <- function(analysisresultsmatrix,
                     region = "both", method = "intensity",
                     color = "redorange") {
  analysisresultsmatrix <- analysisresultsmatrix[[1]]
  if (is.matrix(analysisresultsmatrix) ==
      FALSE) {
    stop("The input is not a matrix!")
  }

  if (color == "redorange") {
    color2 <- c(0.1, 0.9)
  }
  if (color == "bluegreen") {
    color2 <- c(0.3, 0.6)
  }

  if (region == "promoter") {
    feature <- c("promoter")
    coordinates <- c(2, 5, 8)
  }
  if (region == "enhancer") {
    feature <- c("enhancer")
    coordinates <- c(1, 4, 7)
  }
  if (region == "both") {
    region <- c("enhancer/promoter")
    feature <- c("enhancer", "promoter")
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
  names <- split[!(split %in% c("enhancers"))]
  names <- paste(names, collapse = " ")
  casename <- names

  split <- unlist(strsplit(rownames(analysisresultsmatrix)[4], split = " "))
  names <- split[!(split %in% c("enhancers"))]
  names <- paste(names, collapse = " ")
  referencename <- names

  # this is a way to the name of the 'case'
  # from the analysisresults matrix

  caselabel <- paste(casename, "\n", case)
  controllabel <- paste(referencename,  "\n", reference)
  # case=case+shared
  # reference=reference+shared
  p <- venneuler::venneuler(c(A = c(case),
                              B = c(reference),
                              `A&B` = shared))
  p$labels <- c("", "")
  p$colors <- color2
  graphics::plot(p)

  graphics::text(0.15, 0.6, controllabel, cex = 1.1)
  graphics::text(0.75, 0.4, caselabel, cex = 1.1)
  graphics::text(0.5, 0.5, shared, cex = 1.1)
  title <- paste(method, region)

  graphics::title(title, cex = 1.3)

  return(p)
}



#' Plots venn diagrams for comparison of two methods of identifying altered
#' regulatory regions Makes venn diagrams for enhancers, promoters, and combined
#' for both intensity-based peaks and for peaks identified by hotspot calling
#' algorithms.  There is no return value. Six venn diagrams will be plotted
#' @param analysisresultsmatrix analysisresults of countanalysis function
#' place into a a analysisresults matrix by the analyzeanalysisresults function
#' @examples
#' \dontrun{
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks = samplePeaks, minreps = 2)
#' plotConsensusPeaks(samplepeaks = consPeaks)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                            TSS = TSSannot,
#'                                            merge = TRUE,
#'                                            regionspecific = TRUE,
#'                                            mergedistenh = 1500,
#'                                            mergedistprom = 1000 )
#' counts_consPeaks <- getcounts(annotpeaks = consPeaksAnnotated,
#'                               csvfile = csvfile,
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
#  analysisresultsmatrix <- analysisresultsmatrix[[1]]

  if (is.matrix(analysisresultsmatrix[[1]]) ==
      FALSE) {
    stop("The input is not a matrix!")
  }
  graphics::par(mfrow = c(2, 3), oma = c(0,0,2,0))
  plot1 <- plotvenn(analysisresultsmatrix,
                    "promoter", "intensity", "bluegreen")
  plot2 <- plotvenn(analysisresultsmatrix,
                    "enhancer", "intensity", "bluegreen")
  plot3 <- plotvenn(analysisresultsmatrix,
                    "both", "intensity", "bluegreen")
  plot4 <- plotvenn(analysisresultsmatrix,
                    "promoter", "peak", "redorange")
  plot5 <- plotvenn(analysisresultsmatrix,
                    "enhancer", "peak", "redorange")
  plot6 <- plotvenn(analysisresultsmatrix,
                    "both", "peak", "redorange")

  graphics::title("Venn Diagrams Comparing the Two Methods",
                  outer = TRUE,
                  cex.main = 2)

}

