


#' Given the output from getConsensusPeaks, generate a barplot of count statistics
#'
#' @param samplepeaks output generated from getConsensusPeaks
#'
#' @return a ggplot
#'
#' @examples
#' \dontrun{
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' }
#' @export

plotConsensusPeaks <- function(samplepeaks) {
    dfstats = samplepeaks$consPeaksStats
    # dfstats$Replicate=dfstats$PeakType dfstats=dfstats[,-which(colnames(dfstats)=='PeakType')] mydf=na.omit(reshape2::melt(dfstats))
    mydf = tidyr::gather(dfstats, "CellType", "count", 2:3)

    p <- ggplot(data = mydf, aes_string(x = "CellType", y = "count", fill = "PeakType")) + geom_bar(stat = "identity", position = position_dodge()) +
        geom_text(aes_string(label = "count", x = "CellType", y = "count", ymax = "count"), position = position_dodge(width = 1), size = 3, hjust = 0.5,
            vjust = -1.5) + scale_colour_manual(values = c("red", "dark grey")) + theme_bw(base_size = 15) + theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + labs(x = "Sample Type", y = "Number of Peaks (Regulatory Regions)")

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
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot=getTSS()
#'  consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks,
#'  TSS=TSSannot, merge=TRUE, regionspecific=TRUE, mergedistenh=1500,
#'   mergedistprom=1000)
#' plotCombineAnnotatePeaks(consPeaksAnnotated)
#' }
#' @export

plotCombineAnnotatePeaks <- function(conspeaks) {
    mydf = conspeaks$mergestats
    if (nrow(mydf) == 1) {
        stop("No plot to show since merging was not performed when calling combineAnnotatePeaks function")
    } else {
        graphics::par(mfrow = c(1, 2), oma = c(4, 1, 1, 1), mar = c(1, 3, 3, 3))
        toplot = data.frame(Enhancers = mydf$total_number[1:2], Promoters = mydf$total_number[3:4])
        rownames(toplot) = c("NoMerge", "Merge")

        b = graphics::barplot(as.matrix(toplot), ylab = "Number of Regulatory Regions", beside = TRUE, main = "Number of Regulatory Regions\nBefore/After merging",
            ylim = c(0, max(toplot) + max(toplot) * 0.05), col = c("blue4", "coral2"))
        graphics::text(x = b, y = as.matrix(toplot), label = as.matrix(toplot), pos = 3, cex = 0.8)

        toplot = data.frame(Enhancers = mydf$mean_length[1:2], Promoters = mydf$mean_length[3:4])
        rownames(toplot) = c("NoMerge", "Merge")

        b = graphics::barplot(as.matrix(toplot), ylab = "Mean Length of Regulatory Regions", beside = TRUE, main = "Mean length of Regulatory Regions\nBefore/After merging",
            ylim = c(0, max(toplot) + max(toplot) * 0.05), col = c("blue4", "coral2"))
        graphics::text(x = b, y = as.matrix(toplot), label = as.matrix(toplot), pos = 3, cex = 0.8)

        graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        graphics::legend("bottom", c("No Merging", "After Merging"), fill = c("blue4", "coral2"), horiz = TRUE, inset = c(0, 0), xpd = TRUE, bty = "n")
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
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot=getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot,
#' regionspecific=TRUE,mergedistenh=1500,mergedistprom=1000 )
#' counts_consPeaks <-getcounts(annotpeaks=consPeaksAnnotated, csvfile=csvfile,
#' reference='SAEC', chrom='chr21')
#' altre_peaks=countanalysis(counts=counts_consPeaks, pval=0.01, lfcvalue=1)
#' plotCountAnalysis(altre_peaks)
#' }
#' @export

plotCountAnalysis <- function(altrepeaks) {
    toplot = altrepeaks$dftoplot$toplot
    pval = altrepeaks$dftoplot$pval
    lfcvalue = altrepeaks$dftoplot$lfcvalue
    plot = ggplot(toplot, aes(toplot$log2FoldChange, -log2(toplot$padj))) + geom_point(aes(col = factor(toplot$col))) + scale_colour_manual(values = c("red",
        "dark grey")) + theme_bw(base_size = 15) + theme(legend.title = element_blank()) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0,
        0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "log2FC", y = "-log2(pvalue)") + geom_hline(aes(yintercept = -log2(pval)),
        linetype = "dashed") + geom_vline(aes(xintercept = (-lfcvalue)), linetype = "dashed") + geom_vline(aes(xintercept = (lfcvalue)), linetype = "dashed")

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
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot=getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot,merge=TRUE,
#' regionspecific=TRUE,mergedistenh=1500,mergedistprom=1000 )
#' counts_consPeaks <- getcounts(annotpeaks=consPeaksAnnotated, csvfile=csvfile,
#' reference='SAEC', chrom='chr21')
#' plotgetcounts(counts_consPeaks)
#' }
#' @export

plotgetcounts <- function(countsconspeaks) {
    mydf = countsconspeaks$regioncountsforplot
    varstack = suppressMessages(melt(mydf))
    varstack$concat = paste(varstack$region, varstack$variable, sep = ": ")
    varstack$concat = sub("librarysize.*", "", varstack$concat)
    densityplot = ggplot(varstack, aes(x = varstack$value)) + geom_density(aes(group = varstack$concat, color = varstack$concat, fill = varstack$concat),
        alpha = 0.3) + theme_bw(base_size = 15, base_family = "") + theme(legend.title = element_blank()) + scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "log2 read counts \n(normalized by library and region sizes)")

    return(densityplot)
}

##############################################################################

#' Given the output from enrichment(), creates a heatmap from
#' the ouput of the enrichment analysis. Presence or
#' absence of the pathway in enrichment of both type-specific (increased
#' or decreased log2fold change, low p-value) and shared (no change, higher p-value)
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
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot=getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot,merge=TRUE,
#' regionspecific=TRUE,mergedistenh=1500,mergedistprom=1000 )
#' counts_consPeaks <- getcounts(annotpeaks=consPeaksAnnotated, csvfile=csvfile,
#' reference='SAEC')
#'  altre_peaks=countanalysis(counts=counts_consPeaks, pval=0.01, lfcvalue=1)
#' MFenrich=pathenrich(analysisresults=altre_peaks, ontoltype='MF', enrichpvalfilt=0.01)
#' BPenrich=pathenrich(analysisresults=altre_peaks, ontoltype='BP', enrichpvalfilt=0.01)
#' plot1=enrichHeatmap(MFenrich, title='GO:MF, p<0.01')
#' plot2=enrichHeatmap(BPenrich, title='GO:BP, p<0.01')
#' }
#' @export

enrichHeatmap <- function(input, title, pvalfilt = 0.001, removeonlyshared = FALSE) {
    # input=input[[1]]

    if (is.list(input) == FALSE) {
        stop("The input is not a list! Please make sure you are using the output from the enrichment analysis")
    }

    if (is.data.frame(input$expt) == FALSE | is.data.frame(input$reference) == FALSE | is.data.frame(input$shared) == FALSE | length(input) != 3 |
        all(names(input) != c("expt", "reference", "shared"))) {
        stop("The input is not a list of three dataframes! Please make sure you are using the output from the enrichment analysis")
    }

    up = input$expt
    if (length(up) <= 1) {
        up$Description = NA
    } else {
        up = up[up$p.adjust < pvalfilt, ]
    }
    reference = input$reference
    if (length(reference) <= 1) {
        reference$Description = NA
    } else {
        reference = reference[reference$p.adjust < pvalfilt, ]
    }
    shared = input$shared
    if (length(shared) <= 1) {
        shared$Description = NA
    } else {
        shared = shared[shared$p.adjust < pvalfilt, ]
    }

    pathways = unique(c(up$Description, reference$Description, shared$Description))
    print(paste("Pathways", pathways))
    pathways = pathways[!is.na(pathways)]
    if (is.na(pathways) || length(pathways) == 0) {
        stop("No pathways are significant (with adjusted pvalues < user input cutoff)")
    }
    # make a list of all the pathways in up, down, and shared
    heatmapmatrix = matrix(data = NA, nrow = length(pathways), ncol = 3)
    # make a matrix with as many row as there are pathways
    row.names(heatmapmatrix) = pathways
    # name the rows with the pathway names

    colnames(heatmapmatrix) = c("up", "down", "shared")
    # put up, down, and shared as the pathway names

    print(paste("Dim heatmapmatrix", dim(heatmapmatrix)))

    for (i in 1:length(row.names(heatmapmatrix))) {
        print(row.names(heatmapmatrix)[i])
        if (row.names(heatmapmatrix)[i] %in% up$Description) {
            num1 = which(up$Description == row.names(heatmapmatrix)[i])
            heatmapmatrix[i, 1] = up[num1, 6]
        }

        if (row.names(heatmapmatrix)[i] %in% reference$Description) {
            num2 = which(reference$Description == row.names(heatmapmatrix)[i])
            heatmapmatrix[i, 2] = reference[num2, 6]
        }

        if (row.names(heatmapmatrix)[i] %in% shared$Description) {
            num3 = which(shared$Description == row.names(heatmapmatrix)[i])
            heatmapmatrix[i, 3] = shared[num3, 6]
        }
    }

    # places the adjusted p-value in the matrix is there is one

    if (removeonlyshared == TRUE) {
        mycounts = as.numeric(apply(heatmapmatrix, 1, function(x) is.na(x[1]) & is.na(x[2])))
        # finds the shared pathways the are not present in up or down
        heatmapinput = heatmapmatrix[mycounts == 0, ]
        # keeps those that are not only shared
    }
    if (removeonlyshared == FALSE) {
        heatmapinput = heatmapmatrix
    }


    heatmapdata = as.data.frame(heatmapinput)
    heatmapdata = heatmapdata[order(heatmapdata$down, heatmapdata$up, heatmapdata$shared, decreasing = TRUE), ]
    # sorts matrix
    heatmapdata$id = rownames(heatmapdata)
    # makes id
    rownames(heatmapdata) = c(1:nrow(heatmapdata))
    meltedheatmapdata = melt(heatmapdata)

    meltedheatmapdata$id = factor(meltedheatmapdata$id, levels = unique(meltedheatmapdata$id))

    p1 = ggplot(meltedheatmapdata, aes(y = meltedheatmapdata$id, x = meltedheatmapdata$variable)) + geom_tile(aes(fill = meltedheatmapdata$value),
        colour = "black") + scale_fill_continuous(low = "turquoise1", high = "navy", na.value = "white", guide = guide_legend(title = "Pvalue")) +
        theme(text = element_text(size = 13)) + ggtitle(title)
    p2 = p1 + scale_x_discrete(expand = c(0, 0), labels = c("High in Expt", "High in Reference", "Shared")) + scale_y_discrete(expand = c(0, 0)) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank())

    return(p2)
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
#'
#' @examples
#' \dontrun{
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot=getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot,merge=TRUE,
#' regionspecific=TRUE,mergedistenh=1500,mergedistprom=1000 )
#' counts_consPeaks <- getcounts(annotpeaks=consPeaksAnnotated, csvfile=csvfile,
#' reference='SAEC')
#' altre_peaks=countanalysis(counts=counts_consPeaks, pval=0.01, lfcvalue=1)
#' MFenrich=pathenrich(analysisresults=altre_peaks, ontoltype='MF', enrichpvalfilt=0.01)
#' BPenrich=pathenrich(analysisresults=altre_peaks, ontoltype='BP', enrichpvalfilt=0.01)
#' plot1=enrichHeatmap(MFenrich, title='GO:MF, p<0.01')
#' plot2=enrichHeatmap(BPenrich, title='GO:BP, p<0.01')
#' multiplot(plot1,plot2,cols=1)
#' }
#'
#' @export
#'
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel ncol: Number of columns of plots nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
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

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
        }
    }
}
