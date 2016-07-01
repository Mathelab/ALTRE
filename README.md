# ALTRE R Package

## Install From Github

The  package in the private repository can be installed from Github using the personal access token provided in the *auth_token* field

```{R}
devtools::install_github("rfarouni/PackageALTRE", auth_token ="d640aa4cfb2f3eed24af7cd5ffc8d6e486e021c3")
```

### To Run

To launch the Shiny app inside R, run

```{R}
library(ALTRE)
runShinyApp()
```


## Data

A restricted subset of the data with one chromosome (i.e. chromosome 21) can be found on this [page](http://rfarouni.github.io/AltreDataRepo/). The  csv file is included and can be downloaded seperatly [here](https://raw.githubusercontent.com/rfarouni/AltreDataRepo/master/DNaseEncodeWindows.csv)


To download the entire data, please use a *file download manager* to download the files from the links listed below. After you download the files, modify the datapath column of the csv file so that all of the rows contain the file path pointing to the location of the data files on your local machine.

## BAM files:

### A549:
https://www.encodeproject.org/files/ENCFF001CLE/@@download/ENCFF001CLE.bam

https://www.encodeproject.org/files/ENCFF001CLJ/@@download/ENCFF001CLJ.bam
 
### SAEC:
https://www.encodeproject.org/files/ENCFF001EFI/@@download/ENCFF001EFI.bam

https://www.encodeproject.org/files/ENCFF001EFN/@@download/ENCFF001EFN.bam

## BED files:

### A549: 
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseA549HotspotsRep1.broadPeak.gz

http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseA549HotspotsRep2.broadPeak.gz


### SAEC:
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseSaecHotspotsRep1.broadPeak.gz

http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseSaecHotspotsRep2.broadPeak.gz

## Blacklisted genomic regions for functional genomics analysis

https://sites.google.com/site/anshulkundaje/projects/blacklists





##  To build from source (Optional)

To inculde the vignette in the installation package run

```{R}
install.packages("devtools")
biocLite("BiocCheck")
devtools::install(build_vignettes = TRUE)
BiocCheck("/home/rick/Documents/ALTRE")
```

Make sure that the *vignette.Rmd*file contains the following commands in the header:

```{R}
author: "...."
title: "ALTRE: vignette"
date: "`r Sys.Date()`"
output: 
  BiocStyle::html_document:
    toc: true
  BiocStyle::pdf_document:
    toc: true
vignette: >
  %\VignetteIndexEntry{ALTRE: vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
```

and the following command in the body:
```{R}
BiocStyle::markdown()
```
