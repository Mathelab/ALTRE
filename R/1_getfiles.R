#' Read in CSV (internal function)
#' @param csvPath csvPath
#' @export
loadCSVFile <- function(csvPath) {
    stopifnot(is.character(csvPath))
    if (!file.exists(csvPath)) {
        stop("CSV input file does not exist")
    }
    csvfile <- readr::read_csv(csvPath,
                        col_types = readr::cols_only(datapath = readr::col_character(),
                                              bamfiles = readr::col_character(),
                                              peakfiles = readr::col_character(),
                                              sample = readr::col_character(),
                                              replicate = readr::col_character()))

    csvfile <- csvfile[order(csvfile$replicate, csvfile$sample), ]
    return(csvfile)
}

#' Read in BED Files (internal function)
#' @param csvfile csvfile
#' @export
loadBedFiles <- function(csvfile) {
    if (!is(csvfile, "data.frame"))
        stop("csvfile must be a data.frame ")

    readBed <- function(bedPath, ind) {
	trackline=utils::read.table(bedPath,nrows=1,sep="\t")
	# if not track line:
	if (length(grep("track type",trackline$V1)) == 0) { 
	print("No track line")
        	bed <- DataFrame(readr::read_delim(bedPath,
                                    delim = "\t",
                                    col_names = FALSE,
                                    na = "."))[, 1:3]
	}
	else { # if there is a trackline
		print("There's a track")
		                bed <- DataFrame(readr::read_delim(bedPath,
                                delim = "\t",
                                col_names = FALSE,
                                na = ".",skip=1))[, 1:3]
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



#' get TSS file
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
