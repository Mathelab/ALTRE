#' Comparison of methods for identifying altered regulatory regions
#'
#' Creates a table to compare two methods of identifying altered regulatory regions
#' – one based on peak intensity, the other on peak presence as determined by
#' hotspot calling algorithms.
#'
#' @param analysisresults analysisresults of countanalysis.
#' @param samplenames vector of sample types
#' @param reference cell type to be considered "reference" or "reference" to which
#' other cell types will be compared
#' @param lfctypespecific log2fold change for type specific enhancers/promoters
#' @param lfcshared log2fold chance for shared enhancers/promoters
#' @param pvaltypespecific p-value for type specific enhancers/promoters
#' @param pvalshared p-value for shared enhancers/promoters
#'
#' @return matrix comparing the two methods of identifying altered regulatory
#' regions – one based on peak intensity, the other on peak presence as
#' determined by hotspot calling algorithms.
#'
#' ************* analysisresults is categaltre_peaks ***************
#' @examples
#' liver_summary=resultscomparison(analysisresults_liver, reference="Stellate")
#' @export
#'


comparePeaksAltre <- function(analysisresults, 
	samplenames,
	reference,
	lfctypespecific=1.5, 
	lfcshared=1.2, 
	pvaltypespecific=0.01, 
	pvalshared=0.05){
  analysisresults=analysisresults[[1]]

  if (!is.data.frame(analysisresults) ||
	!all.equal(colnames(analysisresults),
		c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj",
		"chr","start","stop","region","A549","SAEC","REaltrecateg"))) {
	stop("Make sure the output of the analysis is from categAltrePeaks() function")
   }

  # removed this bc prone to errors: samplenames=colnames(analysisresults[,11:ncol(analysisresults)])
  analysisresultsmatrix=matrix(nrow=9, ncol=2)
  nonreference=samplenames[!(samplenames %in% reference)]

  if(length(nonreference)==1)
  {string=nonreference
  nonreferencestr=nonreference}
  else{
    string=c()
    for (i in 1:length(nonreference)){
      or="and/or"
      if(i==length(nonreference))
      {string=paste(string,or,nonreference[i], sep=" ")}
      else
      {string=paste(string,nonreference[i], sep=", ")}}
    nonreferencestr=substr(string, 3, nchar(string))
  }



  one=c(rep(nonreferencestr, 3), rep(reference,3), rep("Shared",3))
  two=rep(c("enhancers", "promoters", "all"), 3)

  matrixrownames=paste(one,two)
  rownames(analysisresultsmatrix)=matrixrownames
  colnames(analysisresultsmatrix)=c("intensity", "peak")

  #create a matrix for my analysisresults which says how mant of one-specific total, two-specific promoters, etc. there are

  analysisresultsmatrix[1,1]=length(which(analysisresults$log2FoldChange > lfctypespecific & analysisresults$padj < pvaltypespecific & analysisresults$region == "enhancer"))
  analysisresultsmatrix[2,1]=length(which(analysisresults$log2FoldChange > lfctypespecific & analysisresults$padj < pvaltypespecific & analysisresults$region == "promoter"))
  analysisresultsmatrix[3,1]=length(which(analysisresults$log2FoldChange > lfctypespecific & analysisresults$padj < pvaltypespecific))
  analysisresultsmatrix[4,1]=length(which(analysisresults$log2FoldChange < -lfctypespecific & analysisresults$padj < pvaltypespecific & analysisresults$region == "enhancer"))
  analysisresultsmatrix[5,1]=length(which(analysisresults$log2FoldChange < -lfctypespecific & analysisresults$padj < pvaltypespecific & analysisresults$region == "promoter"))
  analysisresultsmatrix[6,1]=length(which(analysisresults$log2FoldChange < -lfctypespecific & analysisresults$padj < pvaltypespecific))
  analysisresultsmatrix[7,1]= length(which((analysisresults$log2FoldChange >= -lfcshared & analysisresults$log2FoldChange <= lfcshared) & (analysisresults$padj >= pvalshared | is.na(analysisresults$padj)) & analysisresults$region == "enhancer"))
  analysisresultsmatrix[8,1]=length(which((analysisresults$log2FoldChange >= -lfcshared & analysisresults$log2FoldChange <= lfcshared) & (analysisresults$padj >= pvalshared | is.na(analysisresults$padj)) & analysisresults$region == "promoter"))
  analysisresultsmatrix[9,1]=length(which((analysisresults$log2FoldChange >= -lfcshared & analysisresults$log2FoldChange <= lfcshared) & analysisresults$padj >= pvalshared))

  typespecsub=analysisresults[,11:ncol(analysisresults)]

  referencesub=as.matrix(analysisresults[,reference])
  nonreferencesub=as.matrix(analysisresults[,nonreference])

  analysisresultsmatrix[1,2]=enhancernonreference=length(which(apply(nonreferencesub,1,function(x){any(is.na(x)==FALSE)}) & apply(referencesub,1,function(x){any(is.na(x))}) & analysisresults$region == "enhancer" ))
  analysisresultsmatrix[2,2]=promoternonreference=length(which(apply(nonreferencesub,1,function(x){any(is.na(x)==FALSE)}) & apply(referencesub,1,function(x){any(is.na(x))}) & analysisresults$region == "promoter" ))
  analysisresultsmatrix[3,2]=allnonreference=length(which(apply(nonreferencesub,1,function(x){any(is.na(x)==FALSE)}) & apply(referencesub,1,function(x){any(is.na(x))})))
  analysisresultsmatrix[4,2]=enhancerreference=length(which(apply(referencesub,1,function(x){any(is.na(x)==FALSE)}) & apply(nonreferencesub,1,function(x){all(is.na(x))}) & analysisresults$region == "enhancer" ))
  analysisresultsmatrix[5,2]=promoterreference=length(which(apply(referencesub,1,function(x){any(is.na(x)==FALSE)}) & apply(nonreferencesub,1,function(x){all(is.na(x))}) & analysisresults$region == "promoter" ))
  analysisresultsmatrix[6,2]=allreference=length(which(apply(referencesub,1,function(x){any(is.na(x)==FALSE)}) & apply(nonreferencesub,1,function(x){all(is.na(x))})))
  analysisresultsmatrix[7,2]=enhancershared=length(which(apply(referencesub,1,function(x){any(is.na(x)==FALSE)}) & apply(nonreferencesub,1,function(x){any(is.na(x)==FALSE)}) & analysisresults$region == "enhancer" ))
  analysisresultsmatrix[8,2]=promotershared=length(which(apply(referencesub,1,function(x){any(is.na(x)==FALSE)}) & apply(nonreferencesub,1,function(x){any(is.na(x)==FALSE)}) & analysisresults$region == "promoter" ))
  analysisresultsmatrix[9,2]=allshared=length(which(apply(referencesub,1,function(x){any(is.na(x)==FALSE)}) & apply(nonreferencesub,1,function(x){any(is.na(x)==FALSE)})))

  return(list(analysisresultsmatrix))
}




#' Plots a venn digram to compare sample types
#
# Plots the number of type-specific and shared regulatory regions in a venn diagram
# Type of regulatroy region (enhancer, promoter, or both) and type of peak comparison
# (intensity or peak) must be specified.
#
#
# @param analysisresultsmatrix analysisresults of Intensity analysis place into a a analysisresults matrix by the analyzeanalysisresults function
# @param region pick a region, regions can be "enhancer", "promoter", or "both" -- INCLUDE quotes
# @param method pick a method, methods can be "intensity" or "peak" -- include quotes
# @param color include the colors you want in your venn diagram
#
# @return venn diagram
#
# @examples
# plot1=plotvenn(lung_summary, region="enhancer", method="intensity")
#
# @export
#

plotvenn<-function(analysisresultsmatrix, region="both", method="intensity", color="redorange"){
  if (is.matrix(analysisresultsmatrix)==FALSE){
    stop("The input is not a matrix!")
  }

  if(color=="redorange"){
    color2=c(.1,.9)
  }
  if(color=="bluegreen"){
    color2=c(.3,.6)
  }


  if (region == "promoter")
  {feature = c("promoter")
  coordinates=c(2,5,8)}
  if (region == "enhancer")
  {feature = c("enhancer")
  coordinates=c(1,4,7)}
  if (region == "both")
  {   region = c("enhancer/promoter")
  feature = c("enhancer", "promoter")
  coordinates=c(3,6,9)}
  #identifies the correct numbers from the analysisresults matrix based on the regulatory region of interest
  if(method == "intensity"){
    case=analysisresultsmatrix[coordinates[1],1]
    reference=analysisresultsmatrix[coordinates[2],1]
    shared=analysisresultsmatrix[coordinates[3],1]}

  if(method == "peak"){
    case=analysisresultsmatrix[coordinates[1],2]
    reference=analysisresultsmatrix[coordinates[2],2]
    shared=analysisresultsmatrix[coordinates[3],2]}
  #identifies the correct numbers from the analysisresults matrix based on the method of region
  string=paste(rownames(analysisresultsmatrix)[1], rownames(analysisresultsmatrix)[2],
               rownames(analysisresultsmatrix)[3], rownames(analysisresultsmatrix)[4],
               rownames(analysisresultsmatrix)[5], rownames(analysisresultsmatrix)[6],
               rownames(analysisresultsmatrix)[7], rownames(analysisresultsmatrix)[8],
               rownames(analysisresultsmatrix)[9])
  stringsplit=strsplit(string, " ")
  uniquestringsplit=unique(stringsplit[[1]])

  split=unlist(strsplit(rownames(analysisresultsmatrix)[1], split=" "))
  names=split[!(split %in% c("enhancers"))]
  names=paste(names, collapse=" ")
  casename=names

  split=unlist(strsplit(rownames(analysisresultsmatrix)[4], split=" "))
  names=split[!(split %in% c("enhancers"))]
  names=paste(names, collapse=" ")
  referencename=names

  #this is a way to the name of the "case" from the analysisresults matrix

  caselabel=paste(casename, "\n", case)
  controllabel=paste(referencename, "\n", reference)
  #case=case+shared
  #reference=reference+shared
  p=venneuler(c(A=c(case), B=c(reference), "A&B"=shared))
  p$labels=c("","")
  p$colors=color2
  plot(p)

  text(0.15, 0.6, controllabel, cex=1.3)
  text(0.75, 0.4, caselabel, cex=1.3)
  text(0.5, 0.5, shared, cex=1.3)
  title=paste(method, region)

  text(0.5, 0.99, title, cex=1.9)


  return(p)
}



# Plots venn diagrams for comparison of two methods of identifying altered regulatory regions
#
# Makes venn diagrams for enhancers, promoters, and combined for both intensity-based peaks
# and for peaks identified by hotspot calling algorithms. There is no return value, six venn
# diagrams will be plotted
#
# @param analysisresultsmatrix analysisresults of countanalysis function place into a a analysisresults matrix by the analyzeanalysisresults function
#
# @examples
# plotallvenn(lung_summary)
#
# @export


plotallvenn<-function(analysisresultsmatrix){
  analysisresultsmatrix=analysisresultsmatrix[[1]]

  if (is.matrix(analysisresultsmatrix)==FALSE){
    stop("The input is not a matrix!")
  }
  par(mfrow=c(2,3))
  plot1=plotvenn(analysisresultsmatrix, "promoter", "intensity", "bluegreen")
  plot2=plotvenn(analysisresultsmatrix, "enhancer", "intensity", "bluegreen")
  plot3=plotvenn(analysisresultsmatrix, "both", "intensity", "bluegreen")
  plot4=plotvenn(analysisresultsmatrix, "promoter", "peak", "redorange")
  plot5=plotvenn(analysisresultsmatrix, "enhancer", "peak", "redorange")
  plot6=plotvenn(analysisresultsmatrix, "both", "peak", "redorange")

}

#' Creates a custom track for visualization on genome browser
#'
#' Creates a colors coded HOTSPOT (bed) file for visualization in UCSC genome
#' browser: red indicates increased log2fold change/low p-value, blue indicates
#' decreased lod2fold change/low p-value, and purple indicates regions
#' with little to no change and insignificant p-values. Based on the log2fold
#' change and p-value inputs, there is a possibility that some regions will not
#' fulfill any of the supplied criteria ("in-between" shared and type-specific) –
#' they are colored grey.
#'
#' @param analysisresults analysisresults of count analysis
#' @param lfctypespecific log2fold change for type specific enhancers/promoters
#' @param lfcshared log2fold chance for shared enhancers/promoters
#' @param pvaltypespecific p-value for type specific enhancers/promoters
#' @param pvalshared p-value for shared enhancers/promoters
#' @param remove removes filtered out regions
#' @param uncatagorized color of track for uncategorized regions
#' @param referencespecific color of track for reference specific regions
#' @param diseasespecific color of track for disease specific regions
#' @param shared color of track for shared regions
#' @param filteredout color of track for regions that are filtered out
#'
#'
#' @examples
#' genomebrowvisual(analysisresults_liver[[1]])
#'
#' @export
#'

genomebrowvisual<-function(analysisresults, lfctypespecific=1.5, lfcshared=1.2, pvaltypespecific=0.01, pvalshared=0.05, remove=FALSE, uncatagorized="turquoise", referencespecific="red", diseasespecific="blue", shared="purple", filteredout="grey"){
  #if (is.data.frame(analysisresults)==FALSE)
  #{stop("The input is not a dataframe!")

  #}

  colors=c()
  for (type in c(uncatagorized, referencespecific, diseasespecific, shared, filteredout)){
  if (type=="red")
  {type="255,0,0"}
  else if (type=="orange")
  {type="255,165,0"}
  else if (type=="yellow")
  type="255,255,0"
  else if (type=="green")
  type="0,255,0"
  else if (type=="blue")
  type="0,0,255"
  else if (type=="purple")
  type="160,32,240"
  else if (type=="pink")
  type="255,192,203"
  else if (type=="salmon")
  type="250,128,114"
  else if (type=="turquoise")
  type="64,224,208"
  else if (type=="grey")
  type="190,190,190"
  else if (type=="black")
  type="0,0,0"
  else if (type=="white")
  type="255,255,255"
  colors=c(colors, type)
  }

  uncatagorized=colors[1]
  referencespecific=colors[2]
  diseasespecific=colors[3]
  shared=colors[4]
  filteredout=colors[5]

  filename=deparse(substitute(analysisresults))
  bedfile=analysisresults[,c(7:9)]
  bedfile$name=paste(analysisresults[,10], ",l2fc:", round(analysisresults[,2], digits=4), ",pval:", sprintf("%0.7g", analysisresults[,5]), sep="")
  bedfile$score="0"
  bedfile$strand="+"
  bedfile$thickStart=analysisresults$start
  bedfile$thickEnd=analysisresults$start
  #create the columns necessary for bed files
  bedfile$itemRgb=uncatagorized
  bedfile$itemRgb=ifelse((!(is.na(analysisresults$padj))) & analysisresults$log2FoldChange > lfctypespecific & analysisresults$padj < pvaltypespecific, diseasespecific, bedfile$itemRgb)
  #mark cancer specific with blue
  bedfile$itemRgb=ifelse((!(is.na(analysisresults$padj))) & analysisresults$log2FoldChange < -lfctypespecific & analysisresults$padj < pvaltypespecific, referencespecific, bedfile$itemRgb)
  #mark reference specific with red
  bedfile$itemRgb=ifelse((analysisresults$padj > pvalshared) & analysisresults$log2FoldChange > -lfcshared & analysisresults$log2FoldChange < lfcshared & analysisresults$padj > pvalshared, shared, bedfile$itemRgb)
  #mark shared regions with purple
  bedfile$itemRgb=ifelse(is.na(bedfile$itemRgb), filteredout, bedfile$itemRgb)

  if(remove == TRUE){
  bedfile=bedfile[bedfile$itemRgb!=filteredout,]}

  #create the bed file in an object
  trackline=c(paste0('track name="', substitute(filename), '" description="', substitute(filename), 'demonstration" visibility=2 itemRgb="On"'))
  #create the trackline
  write.table(trackline, file=paste0(noquote(substitute(filename)),"genomebrowser.bed"), row.names=FALSE, col.names= FALSE, quote=FALSE)
  write.table(bedfile, file=paste0(noquote(substitute(filename)),"genomebrowser.bed"), row.names=FALSE, col.names= FALSE, append = TRUE, quote=FALSE)
  #create the bedfile in an ACTUAL file
}


