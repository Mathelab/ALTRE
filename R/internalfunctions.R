#' Internal function: converts granges to dataframe
#'
#' Takes a Grange input and converts it into a dataframe
#' Not for use by package user
#'
#' @param grange Grange
#'
#' @return dataframe
#'
#' Dataframe1=grangestodataframe(Grange1)

grangestodataframe <- function(grange) {
    chr = seqnames(grange)
    dataframe = data.frame(chr)
    dataframe = data.frame(lapply(dataframe, as.character), stringsAsFactors = FALSE)
    dataframe$start = start(grange)
    dataframe$stop = end(grange)
    # create the beginning of the dataframe (chr, start, stop, etc...)
    if (length(mcols(grange)) != 0) {
        for (i in 1:length(mcols(grange))) {
            dataframe[, 3 + i] = lapply(mcols(grange)[i], as.character, stringsAsFactors = FALSE)
        }
    }
    # if the additional columns of the grange are not 0, do this....
    colnames(dataframe) = c("chr", "start", "stop", colnames(mcols(grange)))
    # name columns off grange
    return(dataframe)
}

################################################################################ 
#' Internal function: annotates peaks as promoters or enhancers
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
    # read in TSS, make column header, change to Granges Find distance to transcription start site
    distancetoTSS <- distanceToNearest(grange, TSS)
    # make dataframe from grange
    newdataframe = grangestodataframe(grange)
    newdataframe$distance = mcols(distancetoTSS)$distance
    newdataframe = within(newdataframe, {
        region = ifelse(distance <= distancefromTSS, "promoter", "enhancer")
    })
    # annotate anything <=1500 bp away as promoter, otherwise enhancer
    chr = c()
    annotatedgrange <- with(newdataframe, GRanges(chr, IRanges(start, stop), meta = newdataframe[, 5]))
    # create a grange
    colnames(mcols(annotatedgrange)) = c("region")
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

statscombineannotate <- function(tableofinfo, input, row1, row2) {
    inputdata = input
    inputdata$size = inputdata$stop - inputdata$start
    tableofinfo[row1, 2] = round(mean(inputdata[inputdata[, 4] == "enhancer", ]$size), digits = 0)
    tableofinfo[row2, 2] = round(mean(inputdata[inputdata[, 4] == "promoter", ]$size), digits = 0)
    tableofinfo[row1, 1] = nrow(inputdata[inputdata[, 4] == "enhancer", ])
    tableofinfo[row2, 1] = nrow(inputdata[inputdata[, 4] == "promoter", ])
    return(tableofinfo)
}
################################################################################ 

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

rundose <- function(set, background, log2FoldChange, pvalfilt, ontoltype, genes, offspring) {
    TSSgranges = getTSS()
    en2eg = as.list(org.Hs.egENSEMBL2EG)
    en2egfunc <- function(num) {
        en2eg[[num]]
    }
    
    setgranges <- GRanges(set$chr, IRanges(set$start, set$stop), meta = set$log2FoldChange)
    setwithgene = TSSgranges[nearest(setgranges, TSSgranges)]
    # colnames(mcols(setwithgene))=c('genename') mcols(setwithgene)[2]=mcols(setgranges)[1] colnames(mcols(setwithgene))=c('gene', 'log2foldchange')
    # setwithgenelist=unique(unlist(as.list(mcols(setwithgene)[1]))) setentrezlist=lapply(setwithgenelist,en2egfunc)
    # names(setentrezlist)=setwithgenelist setentrezlist=unlist(setentrezlist)
    setentrezlist = unique(unlist(as.list(mcols(setwithgene)["gene_id"])))
    
    bggranges <- GRanges(background$chr, IRanges(background$start, background$stop), meta = background$log2FoldChange)
    bgwithgene = TSSgranges[nearest(bggranges, TSSgranges)]
    # colnames(mcols(bgwithgene))=c('genename') bgwithgenelist=unique(unlist(as.list(mcols(bgwithgene)[1])))
    # bgentrezlist=lapply(bgwithgenelist,en2egfunc) names(bgentrezlist)=bgwithgenelist bgentrezlist=unlist(bgentrezlist)
    bgentrezlist = unique(unlist(as.list(mcols(bgwithgene)["gene_id"])))
    
    result = enrichGO(setentrezlist, OrgDb = "org.Hs.eg.db", pvalueCutoff = 1, keytype = "ENSEMBL", qvalueCutoff = 1, pAdjustMethod = "fdr", universe = bgentrezlist, 
        ont = ontoltype)
    
    newresult = summary(result)
    newresult = newresult[order(newresult$p.adjust, decreasing = TRUE), ]
    newresultpval = newresult[newresult$p.adjust < pvalfilt, ]
    
    if (nrow(newresultpval) == 0) {
        return("No significant pathways")
    } else {
        GOtocheck = newresultpval$ID
        
        #### Gene Filtering###########
        GenesforGOterms <- function(Node) {
            GOgenes = as.list(org.Hs.egGO2EG)
            genenum = length(GOgenes[[Node]])
            return(genenum)
        }
        
        genenum = lapply(GOtocheck, GenesforGOterms)
        compare = cbind(as.matrix(newresultpval$ID), as.matrix(newresultpval$Description), as.matrix(genenum))
        filternum = which(compare[, 3] > genes)
        newresultpval2 = newresultpval[filternum, ]
        GOtocheck = newresultpval2$ID
        
        #### Offspring Filtering########
        GOfindtermsatlevel <- function(Node, subont) {
            subontlist = as.list(subont)
            result = subontlist[[Node]]
            if (length(result) == 1) {
                if (is.na(result) == TRUE) {
                  result = NULL
                }
            }
            return(length(result))
        }
        
        if (ontoltype == "MF") {
            offspringnum = lapply(GOtocheck, GOfindtermsatlevel, subont = GO.db::GOMFOFFSPRING)
        }
        if (ontoltype == "BP") {
            offspringnum = lapply(GOtocheck, GOfindtermsatlevel, subont = GO.db::GOBPOFFSPRING)
        }
        if (ontoltype == "CC") {
            offspringnum = lapply(GOtocheck, GOfindtermsatlevel, subont = GO.db::GOCCOFFSPRING)
        }
        
        compare = cbind(as.matrix(newresultpval2$ID), as.matrix(newresultpval2$Description), as.matrix(offspringnum))
        filternum = which(compare[, 3] < offspring)
        genefilt = newresultpval2[filternum, ]
        
        return(genefilt)
    }
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

mergeclosepeaks <- function(peaklist, grange, mergedist, TSS, distancefromTSS) {
    chr = region = c()
    halfdistance = mergedist/2
    
    # how to deal with odd distances entered by users
    if (halfdistance%%1 != 0) {
        halfdistup = halfdistance + 0.5
        halfdistdown = halfdistance - 0.5
    }
    if (halfdistance%%1 == 0) {
        halfdistup = halfdistance
        halfdistdown = halfdistance
    }
    
    # determine the size of each of the regions
    grange$sizeofregionbeforemerge = grange$stop - grange$start
    
    # add (approximately +/- 0.5) the same amount to both the beginning and end of each region these will be subtracted later so it doesn't really
    # matter what they are INDIVIDUALLY It only matters what number they are TOGETHER ('merged distance')
    grange$start = grange$start - halfdistup
    grange$stop = grange$stop + halfdistdown
    
    # reduce the merged regions
    beforemerge <- GRanges(grange$chr, IRanges(grange$start, grange$stop), meta = region)
    aftermerge = reduce(beforemerge)
    
    # annotate with transcription start site (again)
    aftermerge = tssannotgrange(aftermerge, TSS, distancefromTSS)
    
    # get the names of the input files and place them in a vector
    namesvector = c()
    for (i in 1:length(peaklist)) {
        newgranges = peaklist[[i]]
        name = as.character(unique(mcols(newgranges)[1])[1, 1])
        namesvector = c(namesvector, name)
    }
    
    # this will annotate the regions with type-specificity
    for (i in 1:length(peaklist)) {
        typespecific = findOverlaps(aftermerge, peaklist[[i]])
        newdataframe = data.frame(matrix(nrow = length(aftermerge)))
        newdataframe[queryHits(typespecific), 1] = namesvector[i]
        values(aftermerge) = cbind(values(aftermerge), newdataframe)
        colnames(mcols(aftermerge))[i + 1] = c(namesvector[i])
    }
    
    aftermergedata = grangestodataframe(aftermerge)
    aftermergedata$start = aftermergedata$start + halfdistup
    aftermergedata$stop = aftermergedata$stop - halfdistdown
    # remove the extra length added to each region for merging purposes
    lengthdata = ncol(aftermergedata)
    
    # convert back to grange
    aftermerge <- GRanges(aftermergedata$chr, IRanges(aftermergedata$start, aftermergedata$stop), meta = aftermergedata[, c(4:lengthdata)])
    colnames(mcols(aftermerge)) = c("region", namesvector)
    return(aftermerge)
}  # end merging function
