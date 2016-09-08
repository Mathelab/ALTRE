#' Internal function: converts granges to dataframe
#'
#' Takes a Grange input and converts it into a dataframe
#' Not for use by package user
#'
#' @param grange Grange
#'
#' @return dataframe
#'
#' Dataframe1 <- grangestodataframe(Grange1)

grangestodataframe <- function(grange) {
  chr <- seqnames(grange)
  dataframe <- data.frame(chr)
  dataframe <- data.frame(lapply(dataframe, as.character),
                          stringsAsFactors = FALSE)
  dataframe$start <- start(grange)
  dataframe$stop <- end(grange)
  # create the beginning of the dataframe
  # (chr, start, stop, etc...)
  if (length(mcols(grange)) != 0) {
    for (i in 1:length(mcols(grange))) {
      dataframe[, 3 + i] <- lapply(mcols(grange)[i],
                                   as.character,
                                   stringsAsFactors = FALSE)
    }
  }
  # if the additional columns of the grange
  # are not 0, do this....
  colnames(dataframe) <- c("chr",
                           "start",
                           "stop",
                           colnames(mcols(grange)))
  # name columns off grange
  return(dataframe)
}

################################################################################
#' Internal function: annotates peaks as TSS-proximals or TSS-distals
#'
#' Takes a Grange input
#' Not for use by package user
#'
#' @param grange GRanges object
#' @param TSS GRanges of Transcription Start Site Locations
#' @param distancefromTSS distancefromTSS
#' @return GRanges object
#'

tssannotgrange <- function(grange, TSS, distancefromTSS) {
  # read in TSS, make column header, change to
  # Granges Find distance to transcription
  # start site
  distancetoTSS <- distanceToNearest(grange, TSS)
  # make dataframe from grange
  newdataframe <- grangestodataframe(grange)
  newdataframe$distance <- mcols(distancetoTSS)$distance
  newdataframe <- within(newdataframe, {
    region <- ifelse(distance <= distancefromTSS,
                     "TSS-proximal",
                     "TSS-distal")
  })
  # annotate anything <=1500 bp away as
  # TSS-proximal, otherwise TSS-distal
  chr <- c()
  annotatedgrange <- with(newdataframe,
                          GRanges(chr,
                                  IRanges(start, stop),
                                  meta = newdataframe[, 5]))
  # create a grange
  colnames(mcols(annotatedgrange)) <- c("region")
  return(annotatedgrange)
}

################################################################################

#' Internal function: creates a stat matrix for combined and annotated peaks
#'
#' Takes a data frame of GRanges
#' Not for use by package user
#'
#' @param tableofinfo tableinfo
#' @param input input
#' @param row1 row1
#' @param row2 row2
#'
#' @return matrix
#'

statscombineannotate <- function(tableofinfo,
                                 input, row1, row2) {
  inputdata <- input
  inputdata$size <- inputdata$stop - inputdata$start
  tableofinfo[row1, 2] <- round(mean(inputdata[inputdata[, 4] ==
                                                 "TSS-distal", ]$size),
                                digits = 0)
  tableofinfo[row2, 2] <- round(mean(inputdata[inputdata[, 4] ==
                                                 "TSS-proximal", ]$size),
                                digits = 0)
  tableofinfo[row1, 1] <- nrow(inputdata[inputdata[, 4] == "TSS-distal", ])
  tableofinfo[row2, 1] <- nrow(inputdata[inputdata[, 4] == "TSS-proximal", ])
  return(tableofinfo)
}
################################################################################

#' Internal function: merges peaks within user input distance of each other
#'
#' Takes a data frame of GRanges and distance as input
#' Not for use by package user
#'
#' @param grange data frame of GRanges
#' @param mergedist numeric
#' @param peaklist peaklist
#' @param distancefromTSS distancefromTSS
#' @param TSS TSS
#'
#' @return GRanges object
#'

mergeclosepeaks <- function(peaklist, grange,
                            mergedist, TSS, distancefromTSS) {
  chr <- region <- c()
  halfdistance <- mergedist/2

  # how to deal with odd distances entered by
  # users
  if (halfdistance%%1 != 0) {
    halfdistup <- halfdistance + 0.5
    halfdistdown <- halfdistance - 0.5
  }
  if (halfdistance%%1 == 0) {
    halfdistup <- halfdistance
    halfdistdown <- halfdistance
  }

  # determine the size of each of the regions
  grange$sizeofregionbeforemerge <- grange$stop - grange$start

  # add (approximately +/- 0.5) the same
  # amount to both the beginning and end of
  # each region these will be subtracted later
  # so it doesn't really matter what they are
  # INDIVIDUALLY It only matters what number
  # they are TOGETHER ('merged distance')
  grange$start <- grange$start - halfdistup
  grange$stop <- grange$stop + halfdistdown

  # reduce the merged regions
  beforemerge <- GRanges(grange$chr,
                         IRanges(grange$start,
                                             grange$stop),
                         meta = region)
  aftermerge <- reduce(beforemerge)

  # annotate with transcription start site
  # (again)
  aftermerge <- tssannotgrange(aftermerge,
                               TSS,
                               distancefromTSS)

  # get the names of the input files and place
  # them in a vector
  namesvector <- c()
  for (i in 1:length(peaklist)) {
    newgranges <- peaklist[[i]]
    name <- as.character(unique(mcols(newgranges)[1])[1, 1])
    namesvector <- c(namesvector, name)
  }

  # this will annotate the regions with
  # type-specificity
  for (i in 1:length(peaklist)) {
    typespecific <- findOverlaps(aftermerge, peaklist[[i]])
    newdataframe <- data.frame(matrix(nrow = length(aftermerge)))
    newdataframe[queryHits(typespecific),
                 1] <- namesvector[i]
    values(aftermerge) <- cbind(values(aftermerge), newdataframe)
    colnames(mcols(aftermerge))[i + 1] <- c(namesvector[i])
  }

  aftermergedata <- grangestodataframe(aftermerge)
  aftermergedata$start <- aftermergedata$start + halfdistup
  aftermergedata$stop <- aftermergedata$stop - halfdistdown
  # remove the extra length added to each
  # region for merging purposes
  lengthdata <- ncol(aftermergedata)

  # convert back to grange
  aftermerge <- GRanges(aftermergedata$chr,
                        IRanges(aftermergedata$start, aftermergedata$stop),
                        meta = aftermergedata[, c(4:lengthdata)])

  colnames(mcols(aftermerge)) <- c("region",  namesvector)

  return(aftermerge)
}  # end merging function


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

