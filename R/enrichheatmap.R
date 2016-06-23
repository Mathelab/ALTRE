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
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot=getTSS()
#' consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot,merge=TRUE,
#'	regionspecific=TRUE,mergedistenh=1500,mergedistprom=1000 )
#' counts_consPeaks <- getcounts(annotpeaks=consPeaksAnnotated, csvfile=csvfile,
#'	reference="SAEC")
#'  altre_peaks=countanalysis(counts=counts_consPeaks, pval=0.01, lfcvalue=1)
#' MFenrich=pathenrich(analysisresults=altre_peaks, ontoltype="MF", pvalfilt=0.01)
#' BPenrich=pathenrich(analysisresults=altre_peaks, ontoltype="BP", pvalfilt=0.01)
#' plot1=enrichheatmap(MFenrich, title="GO:MF, p<0.01")
#' plot2=enrichheatmap(BPenrich, title="GO:BP, p<0.01")
#' multiplot(plot1,plot2,cols=1)
#' @export

enrichHeatmap <- function(input, title, pvalfilt=0.001, removeonlyshared=FALSE){
  #input=input[[1]]

  if (is.list(input)==FALSE)
  {stop("The input is not a list! Please make sure you are using the output from the enrichment analysis")}

  if (is.data.frame(input$expt)==FALSE | is.data.frame(input$reference)==FALSE | 
	is.data.frame(input$shared)==FALSE | length(input) != 3 |
	all(names(input) != c("expt","reference","shared"))){
		stop("The input is not a list of three dataframes! Please make sure you are using the output from the enrichment analysis")
	}

  up=input$expt
  if (length(up)<=1) {up$Description=NA} else 
	{up=up[up$p.adjust<pvalfilt,]}
  reference=input$reference
  if (length(reference)<=1) {reference$Description=NA} else 
	{reference=reference[reference$p.adjust<pvalfilt,]}
  shared=input$shared
  if(length(shared)<=1) {shared$Description=NA} else {
	shared=shared[shared$p.adjust<pvalfilt,]}

  pathways=unique(c(up$Description, reference$Description, shared$Description))
  print(paste("Pathways",pathways))
  pathways=pathways[!is.na(pathways)]
  if(is.na(pathways) || length(pathways)==0) {
	stop("No pathways are significant (with adjusted pvalues < user input cutoff)")
  }
  #make a list of all the pathways in up, down, and shared
  heatmapmatrix=matrix(data=NA, nrow=length(pathways), ncol=3)
  #make a matrix with as many row as there are pathways
  row.names(heatmapmatrix)=pathways
  #name the rows with the pathway names

  colnames(heatmapmatrix)=c("up","down","shared")
  #put up, down, and shared as the pathway names

  print(paste("Dim heatmapmatrix",dim(heatmapmatrix)))

  for (i in 1:length(row.names(heatmapmatrix))){
	print(row.names(heatmapmatrix)[i])
    if (row.names(heatmapmatrix)[i] %in% up$Description){
      num1=which(up$Description==row.names(heatmapmatrix)[i])
      heatmapmatrix[i,1]=up[num1,6]}

    if (row.names(heatmapmatrix)[i] %in% reference$Description){
      num2=which(reference$Description==row.names(heatmapmatrix)[i])
      heatmapmatrix[i,2]=reference[num2,6]}

    if (row.names(heatmapmatrix)[i] %in% shared$Description){
      num3=which(shared$Description==row.names(heatmapmatrix)[i])
      heatmapmatrix[i,3]=shared[num3,6]}
  }

  #places the adjusted p-value in the matrix is there is one

  if (removeonlyshared == TRUE){
    mycounts=as.numeric(apply(heatmapmatrix,1,function(x) is.na(x[1]) &  is.na(x[2]) ))
    #finds the shared pathways the are not present in up or down
    heatmapinput=heatmapmatrix[mycounts == 0,]
    #keeps those that are not only shared
  }
  if (removeonlyshared == FALSE){
    heatmapinput=heatmapmatrix
  }


  heatmapdata=as.data.frame(heatmapinput)
  heatmapdata=heatmapdata[order(heatmapdata$down, heatmapdata$up, heatmapdata$shared, decreasing = TRUE),]
  #sorts matrix
  heatmapdata$id=rownames(heatmapdata)
  #makes id
  rownames(heatmapdata)=c(1:nrow(heatmapdata))
  meltedheatmapdata=melt(heatmapdata)

  meltedheatmapdata$id=factor(meltedheatmapdata$id, levels=unique(meltedheatmapdata$id))

  p1=ggplot(meltedheatmapdata, aes(y=meltedheatmapdata$id, x=meltedheatmapdata$variable)) +
	geom_tile(aes(fill=meltedheatmapdata$value), colour = "black") + 
	scale_fill_continuous(low="turquoise1", high="navy", na.value = 'white',
	guide=guide_legend(title="Pvalue")) + 
	theme(text = element_text(size=13)) + ggtitle(title) 
  p2=p1 + scale_x_discrete(expand = c(0, 0), labels=c("High in Expt",
	"High in Reference","Shared")) + 
	scale_y_discrete(expand = c(0, 0)) + 
	theme(axis.title.x = element_blank(), axis.title.y = element_blank())

  return(p2)
}
