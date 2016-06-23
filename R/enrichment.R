#' Enrichment analysis to identify putative pathways of interest for further investigation
#'
#' Determine which pathways are overrepresented in the
#' altered promoters and enhancers. Pathways are determined by linking
#' the enhancer/promoter to the nearest gene, and then linking genes to pathways using
#' the gene ontology database. The "gene" argument limits how few genes a pathway can contain,
#' while the "offspring" argument limits how many offspring a pathway can contain.
#' Pathways with low gene counts are less reliable (often false positives), while pathways
#' with many offspring are vague and unlikely to be of much use -- the enrichment of
#' their more precise offspring is the more interesting question.
#' ************CHANGE ANALYSISRESULTS TO ALTREPEAKS!!!!******************
#' ****************CHANGE REGIONSUBSET TO REGIONTYPE!!!!*****************
#' @param analysisresults results from analysis of counts
#' @param enrichpvalfilt adjusted pval for enrichment to filter on (adjusted for multiple testing)
#' @param ontoltype one of three categories: "MF" (molecular function), "CC" (cellular component), "BP" (biological process)
#' @param lfctypespecific log2fold change (of chromatin accessibility) for type specific enhancers/promoters
#' @param lfcshared log2fold change (of chromatin accessibility) for shared enhancers/promoters
#' @param pvaltypespecific p-value (of chromatin accessibility) for type specific enhancers/promoters
#' @param pvalshared p-value (of chromatin accessibility) for shared enhancers/promoters
#' @param regionsubset "promoter" or "enhancer"
#' @param genes minimum number of genes allowable in a pathway
#' @param offspring maximum number of offspring allowable in a pathway
#' @examples
#' \dontrun{
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' TSSannot <- getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot,merge=TRUE,
#'	regionspecific=TRUE,mergedistenh=1500,mergedistprom=1000 )
#' # Need to run getcounts on all chromosomes
#' counts_consPeaks=getcounts(annotpeaks=consPeaksAnnotated, csvfile=csvfile, reference="SAEC")
#' altre_peaks=countanalysis(counts=counts_consPeaks, pval=0.01, lfcvalue=1)
#' MFenrich=pathenrich(analysisresults=altre_peaks, ontoltype="MF", enrichpvalfilt=0.01)
#' BPenrich=pathenrich(analysisresults=altre_peaks, ontoltype="BP", enrichpvalfilt=0.01)
#'}
#' @return dataframe identifying p-values for enriched pathways -- pathways also annotated with additional information
#'
#'
#' @export

pathenrich<-function(analysisresults,
	ontoltype="MF",
	enrichpvalfilt=0.01,
	lfctypespecific=1.5,
	lfcshared=1.2,
	pvaltypespecific=0.01,
	pvalshared=0.05,
	genes=20,
	offspring=300,
	regionsubset="promoter"){

  analysisresults=analysisresults[[1]]

  if (is.data.frame(analysisresults)==FALSE) {
	stop("analysisresults parameter is not in the correct format, make sure you are using the output from countanalysis()")
  }
  analysisresultsdata=as.data.frame(analysisresults)


  if ( regionsubset == "all" ){
	newanalysisresults=analysisresultsdata
  }
  else {
	newanalysisresults=analysisresultsdata[analysisresultsdata$meta.region==regionsubset,]
  }

  # Define regions that are more open, less open, or shared
  up=newanalysisresults[!(is.na(newanalysisresults$padj)) &
	newanalysisresults$log2FoldChange > lfctypespecific &
	newanalysisresults$padj < pvaltypespecific,]
  down=newanalysisresults[!(is.na(newanalysisresults$padj)) &
	newanalysisresults$log2FoldChange < -lfctypespecific &
	newanalysisresults$padj < pvaltypespecific,]
  shared=newanalysisresults[(newanalysisresults$log2FoldChange <= lfcshared &
	newanalysisresults$log2FoldChange >= -lfcshared) &
	(newanalysisresults$padj >= pvalshared | is.na(newanalysisresults$padj)), ]
  all=rbind(up,down,shared)
  subsets=list(up,down,shared,newanalysisresults)
  names(subsets)=c("up","down","shared","all")

  message("finding expt-specific...")
  if(nrow(subsets[["up"]])==0) {expt=as.data.frame("No REs higher in experiment group")}  else {
	expt=rundose(set=subsets[["up"]], background=subsets[["all"]],log2FoldChange=lfctypespecific,
	ontoltype=ontoltype, pvalfilt=enrichpvalfilt, genes=genes, offspring=offspring)
        if(nrow(shared)==0) {
               expt=as.data.frame("No enrichment found for experiment REs")
	}
  }

  message("finding reference-specific...")
  if(nrow(subsets[["down"]])==0) {reference=as.data.frame("No REs higher in reference group")} else {
	reference=rundose(set=subsets[["down"]], background=subsets[["all"]],log2FoldChange=lfctypespecific,
        ontoltype=ontoltype, pvalfilt=enrichpvalfilt, genes=genes, offspring=offspring)
        if(nrow(expt)==0) {
                reference=as.data.frame("No enrichment found for reference REs")
	}
  }

  message("finding shared...")
  if(nrow(subsets[["shared"]])==0) {shared=as.data.frame("No shared REs")} else {
	shared=rundose(set=subsets[["shared"]], background=subsets[["all"]],log2FoldChange=lfcshared,
	        ontoltype=ontoltype, pvalfilt=enrichpvalfilt, genes=genes, offspring=offspring)
	print(paste("Number of rows", nrow(shared)))
	if(nrow(shared)==0) {
		shared=as.data.frame("No enrichment found for shared REs")
	}
  }

#  enrichstats=data.frame(

  allthree=list(expt=expt, reference=reference, shared=shared)

  return(allthree)
}

