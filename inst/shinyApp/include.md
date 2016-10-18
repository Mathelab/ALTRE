---
title: "vignette"
author: "Elizabeth Baskin", "Ewy Mathe", "Rick Farouni"
date: "July 6, 2016"
---
<center> <h1>Welcome to ALTRE</h1> </center>
  
ALTRE(ALTered Regulatory Elements) is an R software package that streamlines and simplifies the analysis of data generated from genome-wide chromatin accessibility assays such as DNase-seq (a.k.a. DHS-seq), FAIRE-seq, ATAC-seq, and THS-seq. The analysis workflow consists of eight steps as displayed in the sidebar.

Chromatin accessibility data maps the location of regulatory elements (REs), including enhancers and promoters.  REs are involved in regulating gene transcription – they control genes and pathways that can be investigated as putative therapeutic targets, or they may serve as targets themselves. Chromatin accessibility assays allows us to identify the location of regulatory regions genome-wide by identifying “open” chromatin (i.e. euchromatin). Identifying regulatory regions that differ between cell types, such as cancerous and noncancerous cell lines and tissues, holds promise for discovering new mechanisms involved in cellular development and disease progression. While assays for defining regulatory regions are well established, currently there are few workflows that guide researchers though the process of analysing data, starting with aligned reads and peak calls, all the way to obtaining meaningful results such as determining a set of putative pathways and genes for further investigation. 

### How to Get Started and Load Data

__*It is important and all data to be analyzed (alignment files, peak files, meta information in CSV format) be in the same folder.*__

To get started, please load in a CSV (comma separated values) file in Step 1.
The CSV file should contain the following four columns:

1. name of bamfiles
2. name of peakfiles
3. sample name
4. replicate number (a minimum of 2 replicates per sample is required)

__NOTE:__ the column names for the CSV file should be "bamfiles", "peakfiles", "sample", and "replicate", in that order. An example CSV file can be accessed <a href="https://raw.githubusercontent.com/mathelab/AltreDataRepo/master/DNaseEncodeExample.csv" target="_blank">here</a> (right click save as).

To test the package, we are providing a restricted subset of the data with one chromosome (i.e. chromosome 21). The data can be found on this <a href="http://mathelab.github.io/ALTREsampledata/" target="_blank">page</a>. To analyze the entire data containing all chromosomes, please use a *file download manager* to download the files from the links listed on the same page. Remember to include all files (alignments, peaks, CSV) in the same folder for analysis with ALTRE.

### Notes on Plots and Processed Data

All plots can be downloaded at a high resolution in PDF, JPEG, SVG, or PNG.  <strong>This functionality works best in Chrome </strong>(in Firefox, the workflow could be halted when downloading). Within the Shiny app, users can control the names of the x-axis label, y-axis label, and main title.  Users can also choose different color schemes.  If users wish to control the size of the x-axis/y-axis labels and titles, they would need to use the command line.

Throughout the workflow, users can download processed data at the major steps, including after annotating peaks (step 3), after defining altered regions (step 6), and after pathway analysis (step 8).

### Questions
For any issues or questions that might arise, please file a new issue on <a href="https://github.com/Mathelab/ALTRE/issues" target="_blank">Github</a>.



