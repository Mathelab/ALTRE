#' Read in open chromatin data peak files (bed format, 0-based)
#'
#' @param csvfile name of CSV file, with complete path, that contains the following 4 columns for each sample 1) complete filepath; 2) name of bamfiles; 3) name of peak files; 4) name or type of sample; 5) name of replicate
#'
#' @return a GRangesList object comprising one GRanges object (peaks) for each sample.
#'
#' @examples
#' dir <- system.file("extdata", package="ALTRE", mustWork=TRUE)
#' csvfile <- file.path(dir, "lung.csv")
#' samplePeaks <- loadPeaks(csvfile)
#'
#' @export

loadPeaks <- function(csvfile) {
  sampleinfo <- loadCSVFile(csvfile)
  samplepeaks <- loadBedFiles(sampleinfo)
  return(samplepeaks)
}
