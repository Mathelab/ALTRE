
#' get TSS file
#' @export
getTSS <- function() {
    edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
    ensembldb::seqlevelsStyle(edb) <- "UCSC"
    
    TSSdb <- ensembldb::promoters(edb, filter = list(ensembldb::SeqnameFilter(paste0("chr", c(1:22, "X", "Y"))), ensembldb::GeneidFilter("ENSG%", "like")), 
        columns = c("tx_seq_start", "tx_id", "tx_biotype", "gene_id", "gene_name"), upstream = 0, downstream = 2)
    
    return(TSSdb)
}
