# ALTRE R Package

[![Join the chat at https://gitter.im/ProjectALTRE/PublicLobby](https://badges.gitter.im/ProjectALTRE/PublicLobby.svg)](https://gitter.im/ProjectALTRE/PublicLobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Join the chat at https://gitter.im/ProjectALTRE/Lobby](https://badges.gitter.im/ProjectALTRE/Lobby.svg)](https://gitter.im/ProjectALTRE/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Build Status](https://travis-ci.org/Mathelab/ALTRE.svg?branch=master)](https://travis-ci.org/Mathelab/ALTRE)
[![Build status](https://ci.appveyor.com/api/projects/status/i7lbh9tl449hvnmj/branch/master?svg=true)](https://ci.appveyor.com/project/Mathelab/altre/branch/master)
[![codecov](https://codecov.io/gh/Mathelab/ALTRE/branch/master/graph/badge.svg)](https://codecov.io/gh/Mathelab/ALTRE)


## Installation From Github


To install ALTRE, run the following code in the R terminal

```{R}
# First, install the Bioconductor packages (dependencies) with these two lines           
source("http://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(c('org.Hs.eg.db', 'EnsDb.Hsapiens.v75', 'GO.db'))
# Second, install the devtools package for installing ALTRE from GitHub
install.packages("devtools") 
# Third, install the ALTRE package 
devtools::install_github("mathelab/ALTRE")
```

If you encounter an error when runing these 4 lines, please do the following:

### On Linux and Mac OS


When installing on Linux, installation might fail if the XML package cannot be installed. Installation failure  can be fixed by installing the libxml2, an XML C parser for Linux. Also to install the devtools R library, you also need to install several system dependencies. On Ubuntu this can be done by running the following line in the terminal:

```{R}
sudo apt-get install libxml2-dev libssl-dev libcurl4-openssl-dev gfortran
```
On Mac OS, the same dependecies can installed using the *brew* command. 


### On Windows


if you get an installation error then first run the following lines of code in the R console:

```{R}
install.packages(c("htmltools","httpuv","evaluate","markdown"))
```

### Installation Walk-through Screencast

![](inst/img/ALTREinstall.gif)


## Running and Launching the Shiny App

To launch the Shiny app inside R, run

```{R}
library(ALTRE)
runShinyApp()
```

### Shiny App How to Run Screencast

![](inst/img/ALTRErun.gif)


### Shiny App Preview


![](inst/img/ALTREprev.gif)

##Vignette 

A vignette (which provides an overview of the package via step-by-step guide through an example dataset) is here:
https://mathelab.github.io/ALTRE/vignette.html

## Data

A restricted subset of the data with one chromosome (i.e. chromosome 21) can be found on this [page](http://mathelab.github.io/ALTREsampledata/). The corresponding CSV file for input into ALTRE can be downloaded [here](https://raw.githubusercontent.com/mathelab/ALTREsampledata/master/DNaseEncodeExample.csv). Be sure that the CSV file and the data files are in the same folder when running analysis with ALTRE.

To download the entire data, please use a file download manager to download the files from the links listed below. 

### Alignment (in BAM format) files:

#### *A549*:
https://www.encodeproject.org/files/ENCFF001CLE/@@download/ENCFF001CLE.bam

https://www.encodeproject.org/files/ENCFF001CLJ/@@download/ENCFF001CLJ.bam
 
#### *SAEC*:
https://www.encodeproject.org/files/ENCFF001EFI/@@download/ENCFF001EFI.bam

https://www.encodeproject.org/files/ENCFF001EFN/@@download/ENCFF001EFN.bam

### Peak/hotspot (in BED format) files:

#### *A549*: 
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseA549HotspotsRep1.broadPeak.gz

http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseA549HotspotsRep2.broadPeak.gz


#### *SAEC*:
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseSaecHotspotsRep1.broadPeak.gz

http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseSaecHotspotsRep2.broadPeak.gz

### Blacklisted genomic regions for functional genomics analysis

https://sites.google.com/site/anshulkundaje/projects/blacklists



