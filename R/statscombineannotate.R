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

statscombineannotate<-function(tableofinfo,input, row1, row2){
      inputdata=input
      inputdata$size=inputdata$stop-inputdata$start
      tableofinfo[row1,2]=round(mean(inputdata[inputdata[,4]=="enhancer",]$size), digits=0)
      tableofinfo[row2,2]=round(mean(inputdata[inputdata[,4]=="promoter",]$size), digits=0)
      tableofinfo[row1,1]=nrow(inputdata[inputdata[,4]=="enhancer",])
      tableofinfo[row2,1]=nrow(inputdata[inputdata[,4]=="promoter",])
      return(tableofinfo)
}
