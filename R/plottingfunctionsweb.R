

#' plot highcharts
#' @export
#'
plotBarplot <- function() {

  samplepeaks <- data.frame(PeakType = c("ConsensusPeaks", "rep1", "rep2") ,
                    CellType = c(rep("A549",3), rep("SAEC",3)) ,
                    Count = c(143315, 276677, 350683,
                            202756, 278071, 328759))

  dat <- samplepeaks %>% split(levels(samplepeaks$CellType))
  p <- highchart() %>%
    hc_title(text = "Peak Counts by Cell Type",
             style = list(color = '#2E1717',
                          fontWeight = 'bold')) %>%
    hc_subtitle(text = "For bioreplicates and their merged consensus track") %>%
    hc_add_series(
      data = dat[[1]]$Count,
      name = names(dat[1]),
      type = "column",
      dataLabels = list(
        enabled = TRUE,
        rotation = 270,
        color = '#FFFFFF',
        y = 40
      )
    ) %>%
    hc_add_series(
      data = dat[[2]]$Count,
      name = names(dat[2]),
      dataLabels = list(
        enabled = TRUE,
        rotation = 270,
        color = '#FFFFFF',
        y = 40
      ),
      type = "column"
    ) %>%
    hc_yAxis(title = list(text = "Peak Counts"),
             labels = list(format = "{value}")) %>%
    hc_xAxis(categories = dat[[1]]$PeakType) %>%
    hc_legend(
      enabled = TRUE,
      layout = "vertical",
      align = "right",
      verticalAlign = "top",
      floating = TRUE,
      x = -20,
      y = 60
    ) %>%
    hc_tooltip(
      headerFormat = "<b>{series.name}_{point.key}</b><br>",
      pointFormat = "{point.y}",
      valueSuffix = ' peaks'
    ) %>%
    hc_exporting(enabled = TRUE)
  return(p)
}

#' plot heatmap
#' @export

plotHeatmap <- function() {
  dat <- data.frame(matrix(
    c(
      0.094,
      0.668,
      0.4153,
      0.4613
      ,
      0.1138,
      -0.3847,
      0.2671,
      0.1529
      ,
      0.1893,
      0.3303,
      0.5821,
      0.2632
      ,
      -0.0102,
      -0.4259,
      -0.5967,
      0.18
      ,
      0.1587,
      0.2948,
      0.153,
      -0.2208
      ,
      -0.4558,
      0.2244,
      0.6619,
      0.0457
      ,
      -0.6241,
      -0.3119,
      0.3642,
      0.2003
      ,
      -0.227,
      0.499,
      0.3067,
      0.3289
      ,
      0.7365,
      -0.0872,
      -0.069,
      -0.4252
      ,
      0.9761,
      0.4355,
      0.8663,
      0.8107
    ),
    ncol = 4
  ))
  x <- stats::dist(dat)

  p <- hchart(x) %>%
    hc_exporting(enabled = TRUE)

  return(p)
}

#' plot boxplot
#' @export
#'
plotBoxplot <- function() {
  categ <- c('A1', 'A2', 'B1', 'B2', 'B3')

  X <- list(
    c(760, 801, 848, 895, 965),
    c(733, 853, 939, 980, 1080),
    c(714, 762, 817, 870, 918),
    c(724, 802, 806, 871, 950),
    c(834, 836, 864, 882, 910)
  )


  p <- highchart() %>%
    hc_title(text = "Highcharts Box Plot Example",
             style = list(color = '#2E1717',
                          fontWeight = 'bold')) %>%
    hc_subtitle(text = "Subtitle") %>%
    hc_plotOptions(
      boxplot = list(
        fillColor = '#F0F0E0',
        lineWidth = 2,
        medianColor = '#0C5DA5',
        medianWidth = 3,
        stemColor = '#A63400',
        stemDashStyle = 'dot',
        stemWidth = 1,
        whiskerColor = '#3D9200',
        whiskerLength = '20%',
        whiskerWidth = 3
      )
    ) %>%
    hc_add_series(data = X,
                  name = 'Observations',
                  type = "boxplot") %>%
    hc_yAxis(title = list(text = "Observations"),
             labels = list(format = "{value}")) %>%
    hc_xAxis(categories = categ, title = "Experiment No.") %>%
    hc_tooltip(headerFormat = "<b>{point.key}</b><br>",
               pointFormat = "{point.y}") %>%
    hc_exporting(enabled = TRUE)
  return(p)
}

#' plot density
#' @export
plotDensity <- function() {
  diamonds <- data.frame(price = stats::rpois(1000, 100))

  p <- hchart(diamonds$price, color = "#147DA3", name = "Price") %>%
    hc_title(text = "Histogram") %>%
    hc_subtitle(text = "You can zoom me") %>%
    hc_exporting(enabled = TRUE)
  return(p)
}

#' plot scatter
#' @export
plotScatter <- function() {

dataraw <-
'DOK6 0.51 1.861e-08 0.0003053
 TBX5 -2.129 5.655e-08 0.0004191
 SLC32A1 0.9003 7.664e-08 0.0004191
 IFITM1 -1.687 3.735e-06 0.006809
 NUP93 0.3659 3.373e-06 0.006809
 EMILIN2 1.534 2.976e-06 0.006809
 TPX2 -0.9974 2.097e-06 0.006809
 LAMA2 -1.425 2.39e-06 0.006809
 CAV2 -1.052 3.213e-06 0.006809
 TNN -1.658 8.973e-06 0.01472
 POU3F4 1.181 1.062e-05 0.01584
 COL13A1 -1.647 1.394e-05 0.01592
 IFITM3 -1.61 1.202e-05 0.01592
 SHISA3 -1.477 1.31e-05 0.01592
 LOC285954 1.05 1.456e-05 0.01592
 VEPH1 1.137 2.211e-05 0.02267
 ARHGAP29 -1.526 3.675e-05 0.03547
 KIAA1755 -1.562 3.972e-05 0.0362
 LAMC3 -1.563 4.29e-05 0.03704
 ITM2A -1.398 4.972e-05 0.04078
 DTHD1 1.54 5.594e-05 0.04371
 RBMS1 -0.9139 6.688e-05 0.04988
 CEBPD -1.202 7.859e-05 0.05606
 DMBX1 -1.425 9.529e-05 0.06486
 PAPLN -1.253 9.883e-05 0.06486
 ADM -1.357 0.0001089 0.06872
 COL2A1 -1.187 0.0001424 0.07794
 HS3ST3A1 -1.004 0.0001388 0.07794
 DYSF -1.03 0.0001425 0.07794
 PI16 1.495 0.0001297 0.07794
 CDC42EP5 -1.355 0.0001581 0.08146
 SLC12A8 -0.9425 0.0001589 0.08146
 ZNF391 -1.024 0.0001913 0.09512
 GALNTL2 1.075 0.0002298 0.1109
 C4orf45 1.288 0.0002472 0.1159
 KIF18B -0.8849 0.0002551 0.1162
 KIF20A -0.9505 0.0002972 0.1318
 PDE1B 1.053 0.0003356 0.1449
 BCAN 1.117 0.0003698 0.1477
 APLNR -1.365 0.000378 0.1477
 CILP -1.11 0.0003582 0.1477
 TEC -1.373 0.0003701 0.1477
 KLF5 -0.8177 0.0004159 0.1578
 ACSS2 -0.5578 0.0004232 0.1578
 RAPGEF2 0.3371 0.0004513 0.1645
 C1orf51 -0.6451 0.0005237 0.1665
 IGF2-AS1 -1.235 0.0005835 0.1665
 RPLP0P2 -1.096 0.0005689 0.1665
 COTL1 -0.7376 0.0005886 0.1665
 MYO1D -0.8454 0.0005529 0.1665
 CIAO1 0.2695 0.000522 0.1665
 POU3F3 0.6857 0.000588 0.1665
 CFLAR -0.9694 0.0005598 0.1665
 BHLHE40 -1.127 0.0004785 0.1665
 PLSCR4 -1.317 0.0004978 0.1665
 HECW1 0.5135 0.0005373 0.1665
 KCNQ3 1.147 0.000483 0.1665
 TIMP1 -1.15 0.0005267 0.1665
 CAV1 -1.115 0.0006722 0.1869
 LTBP4 -0.8186 0.0006991 0.1912
 HDAC1 -0.4542 0.0007203 0.1937
 HSPG2 -1.218 0.0007963 0.1947
 CYR61 -1.102 0.0008071 0.1947
 SAP30 -0.02968 0.904 0.9994
 FBXO8 -0.07336 0.6573 0.9994
 CEP44 -0.1132 0.6782 0.9994
 HPGD -0.3122 0.4238 0.9994
 GLRA3 0.3489 0.3221 0.9994
 ADAM29 0.4294 0.2735 0.9994
 WDR17 0.1681 0.4038 0.9994
 SPATA4 0.4892 0.2052 0.9994
 SPCS3 -0.06255 0.7011 0.9994
 NEIL3 -0.4648 0.1156 0.9994
 AGA -0.4363 0.1089 0.9994
 MGC45800 -0.2157 0.4652 0.9994
 ODZ3 -0.3824 0.1282 0.9994
 FAM92A3 -0.4185 0.282 0.9994
 C4orf38 -0.05564 0.8868 0.9994
 CDKN2AIP 0.07392 0.7821 0.9994
 LOC389247 -0.06532 0.8677 0.9994
 ING2 -0.1931 0.1921 0.9994
 RWDD4 0.08427 0.6963 0.9994
 PDLIM3 0.195 0.5488 0.9994
 GLP1R -0.001237 0.9971 0.9996
 TCP1 -0.00206 0.9931 0.9996
 CCM2 -0.000938 0.9955 0.9996
 WBSCR16 0.001843 0.9928 0.9996
 LOC155060 -0.003811 0.9916 0.9996
 C8orf42 0.001161 0.9975 0.9996
 CSGALNACT1 -0.001458 0.9967 0.9996
 PPAPDC1B 0.0008906 0.9974 0.9996
 LINC00535 0.003013 0.9929 0.9996
 VPS28 0.0005436 0.9965 0.9996
 DNAJC25 -0.000678 0.9974 0.9996
 TSC1 0.001054 0.9962 0.9996
 RLIM -0.001962 0.9932 0.9996
 TAF9B -0.002379 0.9915 0.9996
 DDX3Y -0.001194 0.996 0.9996
 ZNRF2 0.0002994 0.9988 0.9998
 ZNF559 0.0002839 0.9994 0.9998'

datmat <- t(matrix(unlist(strsplit(unlist(strsplit(dataraw,"\n " ))," ")), 4, 100))
dat <- data.frame(log2FoldChange = as.numeric(datmat[ , 2]),
                     pvalue = as.numeric(datmat[ , 3]),
                     padj = as.numeric(datmat[ , 4]),
                     categ = rep(c("A", "B"), c(50, 50)))
rownames(dat) <- datmat[ , 1]


  p <- highchart() %>%
    hc_add_series_scatter(x = -log10(dat$pvalue),
                          y = dat$log2FoldChange,
                          color = dat$categ,
                          label = rownames(dat))  %>%
    hc_tooltip(headerFormat = "",
               pointFormat  = "<b>log2FC</b> = {point.x}<br> <b>-log10pvalue</b> = {point.y}<br>") %>%
    hc_exporting(enabled = TRUE)
  return(p)
}
