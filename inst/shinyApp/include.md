---
title: "vignette"
author: "Elizabeth Baskin", "Ewy Mathe", "Rick Farouni"
date: "July 1, 2016"
---
<center> <h1>Welcome to ALTRE</h1> </center>
  

ALTRE(ALTered Regulatory Elements) is a software package that streamlines and simplifies analysis chromatin accessibility data, such as those derived from FAIRE-seq, ATAC-seq, and DHS-seq.  The analysis workflow consists of seven steps as displayed in the sidebar.

Chromatin accessibility data maps the location of regulatory elements (REs), including enhancers and promoters.  REs are involved in regulating gene transcription – they control genes and pathways that can be investigated as putative therapeutic targets, or they may serve as targets themselves. Chromatin accessibility assays such as FAIRE-seq, ATAC-seq, and DHS-seq can identify the location of regulatory regions genome-wide by identifying “open” chromatin (e.g. euchromatin).  Identifying regulatory regions that differ between cell types, such as cancerous and noncancerous cell lines and tissues, holds promise for identifying new mechanisms involved in cellular development and disease progression. While assays for defining regulatory regions are well established, there are currently few workflows that guide newcomers from aligned reads and peak calls to meaningful results (putative pathways and genes for further investigation). 

#### Questions
For any questions or issues that arise from using ALTRE, please contact XXX

#### How to Get Started and Load Data

To get started, load in a CSV (comma-separated-values) in step 1 on the left hand side.
The CSV file should contain the following five columns:
1. full path for the location of the alignment (BAM format) and peak (BED format) files
2. name of bamfiles
3. name of peakfiles
4. sample name
5. replicate number (a minimum of 2 replicates per sample is required to run the workflow)
NOTE: the column names for the CSV file should be "datapath","bamfiles","peakfiles","sample",and "replicate".  An example CSV file can be accssed [here](https://raw.githubusercontent.com/rfarouni/AltreDataRepo/master/DNaseEncodeWindows.csv).

For convenience and to test out ALTRE, we have created a restricted subset of the data with one chromosome (i.e. chromosome 21) can be found on this [page](http://rfarouni.github.io/AltreDataRepo/). 

To download the corresponding data in its entirety, please use a *file download manager* to download the files from the links listed below. After you download the files, modify the datapath column of the csv file so that all of the rows contain the file paths pointing to the location of the data files on your local machine.

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

## Blacklisted genomic regions for functional genomics analysis (optional, not used)

https://sites.google.com/site/anshulkundaje/projects/blacklists


