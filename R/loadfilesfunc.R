#' Read in CSV (internal function)
#' @param csvPath csvPath
#' @export
loadCSVFile <- function(csvPath) {

  stopifnot(is.character(csvPath))
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

  bamfiles <- file.path(csvfile$datapath, csvfile$bamfiles,fsep="")
  indexfiles <- file.path(csvfile$datapath, paste(csvfile$bamfiles,".bai",sep=""),fsep="")
  bamFiles <- Rsamtools::BamFileList(bamfiles, index=indexfiles,yieldSize = 100000)

return(bamFiles)
}
