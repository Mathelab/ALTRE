#' Creates a custom track for visualization on genome browser
#'
#' Creates a color-coded BED file for visualization of peaks and their ALTRE categories in a genome
#' browser: red indicates increased log2fold change/low p-value, blue indicates
#' decreased lod2fold change/low p-value, and purple indicates regions
#' with little to no change and insignificant p-values. Based on the log2fold
#' change and p-value inputs, there is a possibility that some regions will not
#' fulfill any of the supplied criteria ("in-between" shared and type-specific) –
#' they are colored grey.
#'
#' @param analysisresults output generated from countanalysis() then categAltrePeaks()
#' @param colors colors, in RGB, used for tracks for ambigious, experiment-specific, reference-specific, 
#' and shared peaks (defaults colors are c(grey, salmon, green, blue), respectively)
#' @param trackname name of track
#' @param trackdescript description of track
#'
#' @examples
#' \dontrun{
#' dir <- system.file('extdata', package='ALTRE', mustWork=TRUE)
#' csvfile <- file.path(dir, 'lung.csv')
#' sampleinfo <- loadCSVFile(csvfile)
#' samplePeaks <- loadBedFiles(sampleinfo)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consPeaks,
#'                                           TSS = TSSannot,
#'                                           merge = TRUE,
#'                                           regionspecific = TRUE,
#'                                           mergedistenh = 1500,
#'                                           mergedistprom = 1000 )
#' #Need to run getcounts on all chromosomes
#' counts_consPeaks <- getcounts(annotpeaks = consPeaksAnnotated,
#'                              csvfile = csvfile,
#'                              reference = 'SAEC')
#' altre_peaks <- countanalysis(counts=counts_consPeaks,
#'                              pval=0.01,
#'                              lfcvalue=1)
#' categaltre_peaks=categAltrePeaks(altre_peaks, lfctypespecific=1.5,
#'	lfcshared=1.2,pvaltypespecific=0.01,pvalshared=0.05)
#' genomebrowvisual(categaltre_peaks)
#' }
#' @export

genomebrowvisual<-function(analysisresults, colors=c("grey","salmon","green","blue"),
	trackname="ALTRE_categories", trackdescript="ALTRE categories") {

	analysisresults=analysisresults[[1]]
	if (is.data.frame(analysisresults)==FALSE) {
		stop("The input for the analysisresults arguement is not a dataframe!")
	}

bedfile=data.frame(chr=analysisresults$chr,
	start=analysisresults$start,
	end=analysisresults$stop,
	name=paste0("peak",1:nrow(analysisresults)),
	score="0",
	strand="+",
	thickStart=analysisresults$start,
	thickEnd=analysisresults$stop)
mycol=analysisresults$REaltrecateg
mycol[which(mycol=="Ambiguous")]=paste(as.numeric(grDevices::col2rgb(colors[1])),sep=",",collapse=",")
mycol[which(mycol=="Shared")]=paste(as.numeric(grDevices::col2rgb(colors[2])),sep=",",collapse=",")
mycol[which(mycol=="Reference Specific")]=paste(as.numeric(grDevices::col2rgb(colors[3])),sep=",",collapse=",")
mycol[which(mycol=="Experiment Specific")]=paste(as.numeric(grDevices::col2rgb(colors[4])),sep=",",collapse=",")
bedfile$itemRgb=mycol

trackline=c(paste0('track name="', substitute(trackname), '" description="', 
	substitute(trackdescript), 'demonstration" visibility=2 itemRgb="On"'))

  outfile="ALTREtrack.bed"
  utils::write.table(trackline, file=outfile, row.names=FALSE, col.names= FALSE, quote=FALSE)
  utils::write.table(bedfile, file=outfile, row.names=FALSE, col.names= FALSE, append = TRUE, quote=FALSE,sep="\t")

} # end of function
