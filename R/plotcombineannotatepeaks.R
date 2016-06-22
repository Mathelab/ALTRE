#' Given the output from combineAnnotatePeaks, plot a barplot showing number of peaks before/after merging 
#' (only works if peaks were merged)
#'
#' @param conspeaks output generated from combineAnnotatePeaks
#'
#' @return a ggplot
#'
#' @examples
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#' consPeaks <- getConsensusPeaks(samplepeaks=samplePeaks,minreps=2)
#' plotConsensusPeaks(samplepeaks=consPeaks)
#' TSSannot=getTSS()
#'  consPeaksAnnotated <- combineAnnotatePeaks(conspeaks=consPeaks, TSS=TSSannot,merge=TRUE,
#'	regionspecific=TRUE,mergedistenh=1500,mergedistprom=1000 )
#' plotCombineAnnotatePeaks(consPeaksAnnotated)
#'
#' @export

plotCombineAnnotatePeaks <- function(conspeaks) {
    mydf=conspeaks$mergestats
    if (nrow(mydf)==1) {
	stop("No plot to show since merging was not performed when calling combineAnnotatePeaks function")
    }
    else {
	graphics::par(mfrow=c(1,2),oma=c(4,1,1,1),mar=c(1,3,3,3))
	toplot=data.frame(Enhancers=mydf$total_number[1:2],
		Promoters=mydf$total_number[3:4])
	rownames(toplot)=c("NoMerge","Merge")

	b=graphics::barplot(as.matrix(toplot),ylab="Number of Regulatory Regions",beside=TRUE,
		main="Number of Regulatory Regions\nBefore/After merging",
		ylim=c(0,max(toplot)+max(toplot)*0.05),col=c("blue4","coral2"))
	graphics::text(x = b, y = as.matrix(toplot), label = as.matrix(toplot), pos = 3, cex = 0.8)

	toplot=data.frame(Enhancers=mydf$mean_length[1:2],
                Promoters=mydf$mean_length[3:4])
        rownames(toplot)=c("NoMerge","Merge")

        b=graphics::barplot(as.matrix(toplot),ylab="Mean Length of Regulatory Regions",beside=TRUE,
                main="Mean length of Regulatory Regions\nBefore/After merging",
                ylim=c(0,max(toplot)+max(toplot)*0.05),col=c("blue4","coral2"))
        graphics::text(x = b, y = as.matrix(toplot), label = as.matrix(toplot), pos = 3, cex = 0.8)  

	graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
	graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
 	graphics::legend("bottom",c("No Merging", "After Merging"),fill=c("blue4","coral2"),
		horiz=TRUE,inset=c(0,0),xpd=TRUE,bty="n")
     }
}
