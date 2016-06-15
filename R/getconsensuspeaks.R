#' Given a GRangeslist object with peaks for each samples, determine the consensus peaks (found in at least N replicates, where N is input by the user) for each sample type
#'
#' @param samplePeaks A GRangesList object comprising one GRanges object (peaks) for each sample (output of loadPeaks() function)

#' @param minreps vector with minimum number of replicates that a peak should be contained in to be called as a consensus peak for each sample.  "minreps" for each sample type should be lsited in alphabetical order (e.g. A549 would be in pos 1, SAEC in pos2...)
#'
#' @return a list comprising:
#' 	1) a GRangeslist with one GRange for each sample type which contains consensus peaks
#'	2) a summary statistic table
#'
#' @examples
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#'  
#' @export

getConsensusPeaks <- function(samplepeaks, minreps) {
  if(class(hotspots)!="GRangesList")
	stop("Peaks must be a GRangesList Object")

  sampnames=names(samplePeaks)
  sampletypes=sort(unique(gsub("_.*","",sampnames)) )

  conspeaks=GRangesList()
  hotspotsinreps=c()
  conspeaks_stats=list()

  for (mytype in 1:length(sampletypes)) {
	mytypepeaks <- samplePeaks[grep(sampletypes[mytype],sampnames)]

	 # Concatenate all peaks pertaining to the same sample type and merge peaks:
	 allregregions=c(mytypepeaks[[1]])
	 for (i in 2:length(mytypepeaks)){
     		 allregregions=c(allregregions, mytypepeaks[[i]])
    	 }
    	 reducedallregregions=reduce(allregregions)
	
	# For each reduced peak, determine whether it was present in each sample type
	for (i in 1:length(mytypepeaks)){
	   typespecific=findOverlaps(reducedallregregions, mytypepeaks[[i]])
	   newdataframe=data.frame(i=matrix(nrow=length(reducedallregregions)))
	   newdataframe[queryHits(typespecific),1]="present"
	   values(reducedallregregions)=cbind(values(reducedallregregions),newdataframe)
        }
	colnames(values(reducedallregregions))=names(mytypepeaks)

	# Find regions that are present in at least N replicates (from user in put minreps)
	reducedallregionsdata=grangestodataframe(reducedallregregions)
	applymatrix=as.matrix(reducedallregionsdata[4:ncol(reducedallregionsdata)])
	keepers=which(apply(applymatrix, 1, function(x) length(which(x=="present")) >= minreps[i]))

	reducedallregionsdatakeepers=reducedallregionsdata[keepers,]

	# Convert back to GRanges object
	finalgranges=with(reducedallregionsdatakeepers, GRanges(chr, IRanges(start, stop)))
	mcols(finalgranges)[1]=sampletypes[mytype]
	colnames(mcols(finalgranges))="sampletype"
	
	# Construct output
	conspeaks[[mytype]]=finalgranges

	#these are the metrics for the hotspots
	hotspotsinreps=c(hotspotsinreps, nrow(reducedallregionsdatakeepers))
	names(hotspotsinreps)=sampletypes[mytype]
	conspeaks_stats[[mytype]]=hotspotsinreps		
  } # end looping through sampletypes

  return(list(conspeaks,conspeaks_stats))
}

