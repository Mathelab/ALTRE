#' Read in CSV file
#'
#' The metadata associated with data files to be analyzed in ALTRE is supplied
#' as a CSV file. The software will automatically retrieve the file path of
#' input CSV so it is important that all analysis files are in the same folder
#' as CSV file. If there are greater than two sample types in the file,
#' two samples must be selected for analysis using the paramaters sample1 and
#' sample2. This software can only compare two samples at a time.
#'
#' @param csvPath csvPath
#' @param sample1 optional, select first sample to use if >2 samples are in csv file
#' @param sample2 optional, select second sample to use if >2 samples are in csv file
#' @return dataframe of CSV file
#'
#' @examples
#' \dontrun{
#' csvfile <- loadCSVFile(csvPath = "DNAseEncodeExample.csv")
#' }
#'
#' @export
loadCSVFile <- function(csvPath, sample1 = NULL, sample2 = NULL) {

    stopifnot(is.character(csvPath))

    if (!file.exists(csvPath)) {
        stop("CSV input file does not exist")
    }
    csvfile <- readr::read_csv(csvPath,
                        col_types = readr::cols_only(
                          bamfiles = readr::col_character(),
                          peakfiles = readr::col_character(),
                          sample = readr::col_character(),
                          replicate = readr::col_character()
                          )
                        )


    if ( length(unique(csvfile$sample)) > 2 && (is.null(sample1) || is.null(sample2)) ) {
      stop("If there are greater than two sample types in the file, two samples must be selected for analysis using the sample1 and sample2 parameters. This software can only
           compare two samples at a time.")
    } else if ( length(unique(csvfile$sample)) > 2) {
      csvfile = csvfile[csvfile$sample %in% c(sample1, sample2), ]
    }

    if (ncol(csvfile) < 4) {
      stop("Columns are missing in the CSV file.  Check the format.")
    } else {
      csvfile <- csvfile[order(csvfile$replicate, csvfile$sample),]

      #If the current directory is used (no path given), then get the directory first:
      if (length(grep("/", csvPath)) + length(grep("\\*", csvPath)) == 0) {
        csvPath <- paste0(getwd(), "/", csvPath)
      }
      csvfile$datapath <- rep(gsub("(.*)\\/(.*)", "\\1", csvPath), nrow(csvfile))

      # Check to be sure that at least 2 replicates exist for each condition
      if ( min(table(csvfile$sample)) < 2 ) {
	stop("One of the conditions has less than 2 replicates.
	     ALTRE requires at least 2 replicates per condition")
      }

      # Check that peakfiles exist
      if (!all(file.exists(base::file.path(csvfile$datapath,csvfile$peakfiles)))) {
	   stop("One of the 'peakfiles' does not exist")
      }

      # Check that bamfiles exist
      if (!all(file.exists(base::file.path(csvfile$datapath,csvfile$bamfiles)))) {
           stop("One of the 'bamfiles' in your CSV does not exist")
      }

     # If "_" are present, remove them, or else downstream code will break
     csvfile$sample=gsub("_","",csvfile$sample)

      return(csvfile)
    }
}

#' Read in BED Files
#'
#' Read in the peak files (BED format) with the loadBedFiles() function. Only
#' the first three columns (chr, start, end) of the peak files are required an
#' read in. Additional columns are allowed but ignored.
#'
#' @param csvfile csvfile from loadCSVFile function
#' @return coordinates of peak files
#'
#' @examples
#' \dontrun{
#' csvfile <- loadCSVFile("DNAseEncodeExample.csv")
#' samplePeaks <- loadBedFiles(csvfile)
#' }
#'
#'
#' @export
loadBedFiles <- function(csvfile) {
    if (!is(csvfile, "data.frame"))
        stop("csvfile must be a data.frame ")

    readBed <- function(bedPath, ind) {

	trackline <- utils::read.table(bedPath, nrows = 1, sep = "\t")
	# if not track line:
	if (length(grep("track type", trackline$V1)) == 0) {
        	bed <- DataFrame(readr::read_delim(bedPath,
                                    delim = "\t",
                                    col_names = FALSE,
                                    na = "."))[, 1:3]
	}
	else {# if there is a trackline
		                bed <- DataFrame(readr::read_delim(bedPath,
                                delim = "\t",
                                col_names = FALSE,
                                na = ".", skip = 1))[, 1:3]
	}
        colnames(bed) <- c("seqnames", "start", "end")
        bed <- DataFrame(bed, csvfile[ind, c("sample", "replicate")])
        bed <- within(bed, {
            start <- start + 1L
            sample <- factor(sample)
            replicate <- factor(replicate)
        })
        return(bed)
    }

    bedFilesPath <- file.path(csvfile$datapath, csvfile$peakfiles)
    bedFiles <- mapply(readBed, bedFilesPath, seq_along(bedFilesPath))
    names(bedFiles) <- paste(csvfile$sample, csvfile$replicate, sep = "_")
    hotspots <- lapply(bedFiles, function(x) as(x, "GRanges"))

    return(GRangesList(hotspots))
}

#' Read in BAM files (internal function)
#' @param csvfile csvfile
#' @export
loadBamFiles <- function(csvfile) {
    if (!is(csvfile, "data.frame"))
        stop("csvfile must be a data.frame ")

    bamfiles <- file.path(csvfile$datapath, csvfile$bamfiles)
    if (!all(file.exists(bamfiles))) {
        stop("bamfiles with the specified paths do not exist; fix CSV file")
    }
    indexfiles <- file.path(csvfile$datapath,
                            paste(csvfile$bamfiles,
                                  ".bai",
                                  sep = ""))
    if (!all(file.exists(indexfiles))) {
        bamFiles <- Rsamtools::BamFileList(bamfiles,
                                          yieldSize = 1e+05)
    } else {
        bamFiles <- Rsamtools::BamFileList(bamfiles,
                                           index = indexfiles,
                                           yieldSize = 1e+05)
    }
    return(bamFiles)
}



#' get TSS file for use in combineAnnotatePeaks function.
#' @export
getTSS <- function() {
    edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
    ensembldb::seqlevelsStyle(edb) <- "UCSC"
    TSSdb <- ensembldb::promoters(edb,
                                  filter = list(
                                    ensembldb::SeqnameFilter(
                                      paste0("chr", c(1:22, "X", "Y"))),
                                    ensembldb::GeneidFilter("ENSG%", "like")),
        columns = c("tx_seq_start",
                    "tx_id",
                    "tx_biotype",
                    "gene_id",
                    "gene_name"),
        upstream = 0,
        downstream = 2)
    return(TSSdb)
}
