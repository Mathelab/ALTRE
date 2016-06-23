#' Read in CSV (internal function)
#' @param csvPath csvPath
#' @export
loadCSVFile <- function(csvPath) {

  stopifnot(is.character(csvPath))
  if(!file.exists(csvPath)) {
	stop("CSV input file does not exist")
  }
  csvfile <- read_csv(
    csvPath,
    col_types = cols_only(
      datapath  = col_character(),
      bamfiles  = col_character(),
      peakfiles = col_character(),
      sample    = col_character(),
      replicate = col_character()
    )
  )

  csvfile=csvfile[order(csvfile$replicate,csvfile$sample),]
  return(csvfile)
}

#' Read in BED Files (internal function)
#' @param csvfile csvfile
#' @export
loadBedFiles <- function(csvfile) {

  if(!is(csvfile, "data.frame"))
    stop("csvfile must be a data.frame ")

  readBed <- function(bedPath, ind) {
    bed <- DataFrame(read_delim(
      bedPath,
      delim = "\t",
      col_names = FALSE,
      na = "."
    ))[, 1:3]
    colnames(bed) <- c("seqnames", "start", "end")
    bed           <- DataFrame(bed, csvfile[ind, c("sample","replicate")])
    bed <- within(bed, {
      start      <- start + 1L
      sample     <- factor(sample)
      replicate  <- factor(replicate)
    })
    return(bed)
  }

  bedFilesPath    <- file.path(csvfile$datapath, csvfile$peakfiles)
  bedFiles        <- mapply(readBed, bedFilesPath, seq_along(bedFilesPath))
  names(bedFiles) <- paste(csvfile$sample, csvfile$replicate,sep="_")
  hotspots        <- lapply(bedFiles, function(x) as(x, "GRanges"))

  return(GRangesList(hotspots))
}

#' Read in BAM files (internal function)
#' @param csvfile csvfile
#' @export
loadBamFiles <- function(csvfile) {
  if(!is(csvfile, "data.frame"))
    stop("csvfile must be a data.frame ")

  bamfiles <- file.path(csvfile$datapath, csvfile$bamfiles)
  if(!all(file.exists(bamfiles))) {
	stop("bamfiles with the specified paths do not exist; fix CSV file")
  }
  indexfiles <- file.path(csvfile$datapath, paste(csvfile$bamfiles,".bai",sep=""))
  if(!all(file.exists(indexfiles))) {
	bamFiles <- Rsamtools::BamFileList(bamfiles,yieldSize = 100000)
  }
  else {
  	bamFiles <- Rsamtools::BamFileList(bamfiles, index=indexfiles,yieldSize = 100000)
  }
return(bamFiles)
}



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
