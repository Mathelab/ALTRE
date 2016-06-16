# Some functions to read in CSV and bedfiles

#####################################
# Read in CSV
#####################################
loadCSVFile <- function(csvPath) {
  stopifnot(is.character(csvPath))
  
  csvfile <- utils::read.csv(
    csvPath
    )
  return(csvfile)
}

#####################################
# Load Bed Files
#####################################
loadBedFiles <- function(csvfile) {
  if(!is(csvfile, "data.frame"))
    stop("csvfile must be a data.frame ")

  readBed <- function(bedPath, ind) {
    bed <- DataFrame(utils::read.delim(
      bedPath,
      header= FALSE,
      na.strings = "."
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

#####################################
# Read in bam files
#####################################
loadBamFiles <- function(csvfile) {

  if(!is(csvfile, "data.frame"))
    stop("csvfile must be a data.frame ")

  bamfiles <- file.path(csvfile$datapath, csvfile$bamfiles,fsep="")
  indexfiles <- file.path(csvfile$datapath, paste(csvfile$bamfiles,".bai",sep=""),fsep="")
  bamFiles <- Rsamtools::BamFileList(bamfiles, index=indexfiles,yieldSize = 100000)

  return(bamFiles)
}
