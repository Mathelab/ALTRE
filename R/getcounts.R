#' Generates count data for regulatory regions
#'
#' Counts the number of reads in each regulatory region for each sample type --
#' read count is derived from user-input BAM filex, and regions of interest
#' are supplied in a Granges object, ideally output of combineannotatepeaks.R.
#'
#' @param annotpeaks list output from combineannotatepeaks function
#' @param csvfile name of CSV file, with complete path, that contains the following 4 columns for each sample 1) complete filepath; 2) name of bamfiles; 3) name of peak files; 4) name or type of sample; 5) name of replicate
#' @param reference name of sample type to be considered "reference" in DESeq2 analysis
#' @param chrom optional, only chromosome chrom will be evaluated
#'
#' @return List containing three items:
#' (1) DESeqDataSet: contains count information for all replicates of all samples
#' (2) Matrix: contains number of enhancers and promoters before and after filtering (if applicable)
#' (3) Plot: contains density plot
#'
#'
#' @examples
#' dir=system.file("extdata", package="ALTRE", mustWork=TRUE)
#' TSSpath=file.path(dir,"Homosapiens_GRCh37.75_TSS.bed")
#' TSSannot=read.table(TSSpath, header=TRUE)
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' consPeaksAnnotated=CombineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot)
#' counts_consPeaks=getcounts(annotpeaks=consPeaksAnnotated, csvfile=csvfile, reference="SAEC", chrom="chr21")
#' @export

getcounts<-function(annotpeaks, csvfile, reference, chrom = NULL){
   sampleinfo <- loadCSVFile(csvfile)
   bamfileslist=loadBamFiles(sampleinfo)  

  if (is.null(chrom) == FALSE){
    inputgranges=annotpeaks[[1]][seqnames(annotpeaks[[1]]) == chrom,]
  }

  else {inputgranges=annotpeaks[[1]]}

  # Count number of reads overlapping each annotated peak
  countsse <- GenomicAlignments::summarizeOverlaps(features=inputgranges, reads=bamfileslist,
                                mode="Union",
                                singleEnd=TRUE,
                                ignore.strand=TRUE)
  #add column labels
  SummarizedExperiment::colData(countsse) <- DataFrame(sampleinfo[,c(1:4)])
  countsse$sample=as.factor(countsse$sample)

  countsse$status <- stats::relevel(countsse$sample, reference)
  countssedds <- DESeqDataSet(countsse, design = ~ sample)

 # Optional filtering out of lowcount regions
 # As part of the DESeq2 algorithm, more stringent filtering will be applied subsequently
  #countssedds[ rowSums(counts(countssedds)) > 1, ]

  #get counts referenceized by librarysize
  normcountssedds=SummarizedExperiment::assay(countssedds, norm = T)

  #get region/peak size
  originaldata=grangestodataframe(inputgranges)
  regionsize=originaldata$stop-originaldata$start

  # Calculate RPKM for plotting densities
  #multiply by 10^6 and divide by regions size to get rpkm
  myrpkm=as.data.frame(normcountssedds[,1]*10^6/regionsize)
  for (i in 2:ncol(normcountssedds)){
    myrpkm[,i]=normcountssedds[,i]*10^6/regionsize
  }
  #take the log2 so that it is a normalized distribution
  myrpkmlog2=log2(as.matrix(myrpkm)+1)
  colnames(myrpkmlog2)=unlist(lapply(sampleinfo$sample, as.character))

  ##########################################
  # Create stats matrix
  # originaldata is created ~ 10 lines lines above
  colnames(originaldata)=unlist(lapply(colnames(originaldata), gsub, pattern="meta.", replacement=""))
  enhancernum=length(which(originaldata$region=="enhancer"))
  promoternum=length(which(originaldata$region=="promoter"))
 
  statdf=data.frame(Num_Enhancers=enhancernum,
	Num_Promoters=promoternum)

  ##########################################
  # Create densityplot
  region=originaldata$region
  alldata=cbind(myrpkmlog2,as.data.frame(region))
  graphics::par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  varstack=suppressMessages(melt(alldata))
  varstack$concat=paste(varstack$region, varstack$variable, sep=": ")
  varstack$concat=sub("librarysize.*", "", varstack$concat)
  densityplot=ggplot(varstack, aes(x=varstack$value)) +
        geom_density(aes(group=varstack$concat, color=varstack$concat, fill=varstack$concat), alpha=0.3) +
        theme_bw(base_size = 15, base_family = "") +
        theme(legend.title=element_blank()) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank()) +
        labs(x = "log2 read counts \n(normalized by library and region sizes)" )



   return(list(regioncounts=countssedds, regioncountstats=statdf, regioncountsplot=densityplot))
}
