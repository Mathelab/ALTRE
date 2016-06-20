#' Internal function: runs enrichment analysis
#'
#' Not for use by package user
#'
#' @param set set
#' @param background background
#' @param log2FoldChange log2FoldChange
#' @param ontoltype ontoltype
#' @param pvalfilt pvalfilt
#' @param genes genes
#' @param offspring offspring
#'
#' @return dataframe

  rundose <- function(set, background,log2FoldChange,pvalfilt,ontoltype, genes, offspring){
    TSSgranges=getTSS()
    en2eg=as.list(org.Hs.egENSEMBL2EG)
    en2egfunc<-function(num){
      en2eg[[num]]}

    setgranges <- GRanges(set$chr, IRanges(set$start, set$stop), meta=set$log2FoldChange)
    setwithgene=TSSgranges[nearest(setgranges, TSSgranges)]
#    colnames(mcols(setwithgene))=c("genename")
#    mcols(setwithgene)[2]=mcols(setgranges)[1]
#    colnames(mcols(setwithgene))=c("gene", "log2foldchange")
#    setwithgenelist=unique(unlist(as.list(mcols(setwithgene)[1])))
#    setentrezlist=lapply(setwithgenelist,en2egfunc)
#    names(setentrezlist)=setwithgenelist
#    setentrezlist=unlist(setentrezlist)
     setentrezlist=unique(unlist(as.list(mcols(setwithgene)["gene_id"])))

    bggranges <- GRanges(background$chr, IRanges(background$start, background$stop), meta=background$log2FoldChange)
    bgwithgene=TSSgranges[nearest(bggranges, TSSgranges)]
    #colnames(mcols(bgwithgene))=c("genename")
    #bgwithgenelist=unique(unlist(as.list(mcols(bgwithgene)[1])))
    #bgentrezlist=lapply(bgwithgenelist,en2egfunc)
    #names(bgentrezlist)=bgwithgenelist
    #bgentrezlist=unlist(bgentrezlist)
    bgentrezlist=unique(unlist(as.list(mcols(bgwithgene)["gene_id"])))

    result=enrichGO(setentrezlist, OrgDb="org.Hs.eg.db",pvalueCutoff = 1,keytype = 'ENSEMBL', 
	qvalueCutoff = 1, pAdjustMethod = "fdr", universe=bgentrezlist, ont=ontoltype)

    newresult=summary(result)
    newresult=newresult[order(newresult$p.adjust, decreasing=TRUE),]
    newresultpval=newresult[newresult$p.adjust<pvalfilt,]

    if (nrow(newresultpval)==0) {
	return("No significant pathways")
    }

    else {
    GOtocheck=newresultpval$ID

    ####Gene Filtering###########
    GenesforGOterms <- function(Node){
      GOgenes=as.list(org.Hs.egGO2EG)
      genenum=length(GOgenes[[Node]])
      return(genenum)
    }

    genenum=lapply(GOtocheck, GenesforGOterms)
    compare=cbind(as.matrix(newresultpval$ID), as.matrix(newresultpval$Description), as.matrix(genenum))
    filternum=which(compare[,3]>genes)
    newresultpval2=newresultpval[filternum,]
    GOtocheck=newresultpval2$ID

    ####Offspring Filtering########
    GOfindtermsatlevel <- function(Node, subont) {
      subontlist=as.list(subont)
      result=subontlist[[Node]]
      if (length(result)==1){
        if(is.na(result)==TRUE)
        {result=NULL}
      }
      return(length(result))
     }

    if (ontoltype == "MF") {
	offspringnum=lapply(GOtocheck, GOfindtermsatlevel, subont=GO.db::GOMFOFFSPRING)
    }
    if (ontoltype == "BP") {
	offspringnum=lapply(GOtocheck, GOfindtermsatlevel, subont=GO.db::GOBPOFFSPRING)
    }
    if (ontoltype == "CC") {
	offspringnum=lapply(GOtocheck, GOfindtermsatlevel, subont=GO.db::GOCCOFFSPRING)
    }

    compare=cbind(as.matrix(newresultpval2$ID), as.matrix(newresultpval2$Description), as.matrix(offspringnum))
    filternum=which(compare[,3]<offspring)
    genefilt=newresultpval2[filternum,]

    return(genefilt)
    }
  }





