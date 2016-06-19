#' Combine and annotate peaks from different sample types. Optionally merge nearby regions.
#'
#' This function accomplishes three tasks:
#' (1) Combines peaks from different sample types into one master list, and annotates
#' each peak with it's sample type-specificity (which cell or tissue types the peak can
#' be found in)
#' (2) Categorizes peaks as either a promoter or enhancer (promoter distance
#' is defined as 1500bp away from a transcription start site (TSS) by default, but the distance
#' can be changed with the distancefromTSS argument)
#' (3) Optionally, regulatory regions that are within a certain distance of each other
#' can be merged to form a larger regulatory region.
#'
#'
#' @param conspeaks list of GRanges objects for each sample type (output by getConsensusPeaks() function)
#' @param TSS file of transcription start sites
#' @param distancefromTSS in bp; peaks within distFromTSS of an annotated Transcription Start Site (TSS) will be annotated as a promoter (default=1500bp)
#' @param merge whether or not regions should be merged if they are within a user set distance to each other (default is FALSE)
#' @param regionspecific logical to if TRUE, merging occurs within same type peaks (e.g. merge promoters, then merge enhancers)
#' @param mergedist merge promoters and enhancers if they are < mergedist apart (set when regionspecific is FALSE)
#' @param mergedistprom merge promoters if they are < mergedistprom apart (set when regionspecific is TRUE)
#' @param mergedistenh merge enhancers peaks if they are < mergedistenh apart (set when regionspecific is TRUE)
#'
#' @return List containing two items:
#' (1) GRanges object:  all sample types combined, regions annotated as type-specific
#' and enhancer/promoter specific.
#' (2) Matrix: number and size of enhancers and promoters before and after merging
#' nearby regulatory regions
#'
#'
#' @examples
#' TSSannot <- getTSS()
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' consPeaksAnnotated=combineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot,merge=TRUE,
#'	regionspecific=TRUE,mergedistenh=1500,mergedistprom=1000)
#'
#' @export

combineAnnotatePeaks<-function(conspeaks, TSS, merge=FALSE, mergedistenh=NA, mergedistprom=NA, mergedist=NA,
	regionspecific=NA, distancefromTSS=1500){

  if (class(conspeaks[[1]])[1]!="GRangesList")
  	{stop("The input conspeaks is not in the correct format!  (try to rerun getConsensusPeaks())")}

  peaklist=conspeaks[[1]]

  #################################
  # Aggregate and combine consensus peaks from all sample types
  allregregions=c(peaklist[[1]])
  for (i in 2:length(peaklist)){
    allregregions=c(allregregions, peaklist[[i]])
  }
  reducedallregregions=reduce(allregregions)

  TSSgranges=tssannotgrange(reducedallregregions,TSS,distancefromTSS)

  #################################
  # Annotate combined peaks (TSSgranges) by the cell type they came from

  # If merge is set to false, don't merge
  if(merge == FALSE){
    namesvector=c()
    for (i in 1:length(peaklist)){
      newgranges=peaklist[[i]]
      name=as.character(unique(mcols(newgranges)[1])[1,1])
      namesvector=c(namesvector, name)
    }

    for (i in 1:length(peaklist)){
      typespecific=findOverlaps(TSSgranges, peaklist[[i]])
      newdataframe=data.frame(matrix(nrow=length(TSSgranges)))
      newdataframe[queryHits(typespecific),1]=namesvector[i]
      values(TSSgranges)=cbind(values(TSSgranges),newdataframe)
      colnames(mcols(TSSgranges))[i+1]=c(namesvector[i])
    }

    #this will annotate the regions with type-specificity
    notmerging=grangestodataframe(TSSgranges)
    listtoreturn = list(consPeaksAnnotated=GRanges(notmerging,meta=notmerging[,4:ncol(notmerging)]),
	mergestats=as.data.frame("No merging because mergedistenh and mergedistprom were set to zero"))
  } # end if no merging

  else { # Do the merging
    if (is.na(regionspecific)) {
	stop("If merging, then the regionspecific paramater must be set to TRUE or FALSE")
    }
    # if merging enhancer and promoter regions seperately,
    # then run the merging function on them
    # seperately and then combine them (WITHOUT reducing or you will lose the annotation)
    if (regionspecific == TRUE){
      if(is.na(mergedistprom) || is.na(mergedistenh)) {
	stop("If regionspecific is true, then mergedistprom and mergedistenh must be set")
      }
      dataframeformerge=grangestodataframe(TSSgranges)
      enhancerbeforemergedata=dataframeformerge[dataframeformerge$region=="enhancer",]
      promoterbeforemergedata=dataframeformerge[dataframeformerge$region=="promoter",]

      # Merge enhancers and promoters independently if they're within user defined distances
      enhancerafter=mergeclosepeaks(peaklist,enhancerbeforemergedata, mergedist=mergedistenh, TSS, distancefromTSS)
      promoterafter=mergeclosepeaks(peaklist,promoterbeforemergedata, mergedist=mergedistprom, TSS, distancefromTSS)

      bothafter=sort(sortSeqlevels(c(enhancerafter,promoterafter)))
    }

    # if merging enhancer and promoter regions at the same time,
    # then you just need to run the function once
    if (regionspecific == FALSE){
      if(is.na(mergedist)) {
	stop("If regionspecific is FALSE, then mergedist must be set")
      }
      dataframeformerge=grangestodataframe(TSSgranges)
      #create grange from dataframe

      bothafter=mergeclosepeaks(peaklist,dataframeformerge,mergedist,TSS, distancefromTSS)
    }

    # Create matrix for comparison:
    resultuserinput=grangestodataframe(bothafter)
    result0=dataframeformerge
    tableofinfo=matrix(nrow=4, ncol=2)
    rownames(tableofinfo)=c("enhancers_before_merging","enhancers_after_merging",
	"promoters_before_merging", "promoters_after_merging")
    colnames(tableofinfo)=c("total_number", "mean_length")

    # Callstatscombineannotate internal function to create table:
    tableofinfo=statscombineannotate(tableofinfo,result0,1,3)
    tableofinfo=statscombineannotate(tableofinfo,resultuserinput,2,4)

    listtoreturn=list(consPeaksAnnotated=GRanges(resultuserinput,meta=resultuserinput[,4:ncol(resultuserinput)]),
	mergestats=as.data.frame(tableofinfo))
  }
return(listtoreturn)
}


