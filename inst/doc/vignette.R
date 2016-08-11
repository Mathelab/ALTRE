## ------------------------------------------------------------------------
library(ALTRE)
csvfile <- loadCSVFile("/home/rick/Documents/AltreDataRepo/DNaseEncodeChr21.csv")
csvfile

## ------------------------------------------------------------------------
samplePeaks <- loadBedFiles(csvfile)
samplePeaks

## ------------------------------------------------------------------------
consensusPeaks <- getConsensusPeaks(samplepeaks = samplePeaks,
                                    minreps = 2)

## ------------------------------------------------------------------------
plotConsensusPeaks(samplepeaks = consensusPeaks)

## ------------------------------------------------------------------------
TSSannot <- getTSS()

## ------------------------------------------------------------------------
consensusPeaksAnnotated <- combineAnnotatePeaks(conspeaks = consensusPeaks,
                                           TSS = TSSannot,
                                           merge = TRUE,
                                           regionspecific = TRUE,
                                           mergedistenh = 1500,
                                           mergedistprom = 1000)

## ------------------------------------------------------------------------
plotCombineAnnotatePeaks(consensusPeaksAnnotated)

## ------------------------------------------------------------------------
consensusPeaksCounts <- getCounts(annotpeaks = consensusPeaksAnnotated,
                              sampleinfo = csvfile,
                              reference = 'SAEC',
                              chrom = 'chr21')
plotgetcounts(consensusPeaksCounts)

## ------------------------------------------------------------------------
alteredPeaks <- countanalysis(counts = consensusPeaksCounts,
                             pval = 0.01,
                             lfcvalue = 1)


## ------------------------------------------------------------------------
alteredPeaksCategorized <- categAltrePeaks(alteredPeaks,
                                    lfctypespecific = 1.5,
                                    lfcshared = 1.2,
                                    pvaltypespecific = 0.01,
                                    pvalshared = 0.05)
plotCountAnalysis(alteredPeaksCategorized)

## ------------------------------------------------------------------------
analysisresults <- comparePeaksAltre(alteredPeaksCategorized, reference = "SAEC")

## ------------------------------------------------------------------------
plotallvenn(analysisresults)

## ------------------------------------------------------------------------
plotDistCountAnalysis(alteredPeaksCategorized, consensusPeaksCounts)

## ---- message = FALSE----------------------------------------------------
MFenrich <- pathenrich(analysisresults = alteredPeaksCategorized,
                       ontoltype = 'MF',
                       enrichpvalfilt = 0.99)

## ------------------------------------------------------------------------
enrichHeatmap(MFenrich, title = "GO:MF", pvalfilt = 0.99)

## ---- eval = FALSE-------------------------------------------------------
#  writeBedFile(alteredPeaksCategorized, "output.txt")

## ------------------------------------------------------------------------
sessionInfo()

