
shinyServer(function(input, output, session) {

  ############################################################################
   # Button to Load Data

  rootVolumes <- c(Home = normalizePath("~"), getVolumes()(), WD = '.')


  shinyFileChoose(input,'file',
                  roots = rootVolumes,
                  session = session)

  # Load data

  loadCSVObj <- reactive({
       loadCSVFile(
         req(
           as.character(
             parseFilePaths(
               rootVolumes,
               input$file)$datapath)
           ), input$csvsample1, input$csvsample2
         )

  })

  loadBedObj <- reactive({
      loadBedFiles(req(loadCSVObj()))
  })

  ############################################################################
  # function calls

  getConsensusObj <- eventReactive(input$buttonmerge, {
    withProgress(message = 'In Progress:',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = .1, detail = "Loading Peak Files")
                   loadBedOut <-  req(loadBedObj())
                   setProgress(value = 0.5, detail = "Merging Replicates")
                   consensusPeaksOut <- getConsensusPeaks(loadBedOut,
                                                       req(input$numOverlap))
                   setProgress(value = 1, detail = "Done!")
                   Sys.sleep(0.5)
                 })
    return(consensusPeaksOut)
  })

  combineAnnotateObj <- eventReactive(input$buttonannot, {
    withProgress(message = 'In Progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = .1, detail = "Retrieving TSS File")
                   TSSannot <- getTSS()
                   setProgress(value = 0.2, detail = "Annotating Peaks")
                   combineAnnotateOut <- combineAnnotatePeaks(
                     conspeaks = req(getConsensusObj()),
                     TSS = TSSannot,
                     merge = input$mergeradio,
                     distancefromTSSdist = input$distTSSdist,
                     distancefromTSSprox = input$distTSSprox,
                     mergedist = input$dist,
                     regionspecific = input$regionradio,
                     distancefromTSS = input$distTSS)
                   setProgress(value = 1, detail = "Done!")
                   Sys.sleep(0.5)
                 })
    return(combineAnnotateOut)
  })

  getCountsObj <- eventReactive(input$buttoncounts, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0.1, detail = "Retrieving counts")

                   if (input$operatingradio) {
                   getCountsOut <- getCountsFast(
                     annotpeaks = req(combineAnnotateObj()),
                     sampleinfo = req(loadCSVObj()),
                     reference = input$reference,
                     chrom = input$chr)
                   } else{
                     getCountsOut <- getCounts(
                       annotpeaks = req(combineAnnotateObj()),
                       sampleinfo = req(loadCSVObj()),
                       reference = input$reference,
                       chrom = input$chr)
                   }

                   setProgress(value = 1, detail = "Done!")


                   Sys.sleep(0.5)
                 })
    return(getCountsOut)
  })

  getAlteredObj <- eventReactive(input$buttondefine, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0.2, detail = "annotating peaks")
                   altred <-
                     countanalysis(
                       counts = req(getCountsObj()),
                       pval = input$alpha,
                       lfcvalue = input$lfcThreshold
                     )
                   setProgress(value = 1, detail = "Done!")
                   Sys.sleep(0.5)
                 })
    return(altred)
  })

  observeEvent(input$buttondefine, {
    getAlteredObj()
  })

  categAltreObj <- eventReactive(input$buttoncat, {
                   categAltreOut <-
                     categAltrePeaks(
                       analysisresults = req(getAlteredObj()),
                       lfctypespecific = input$lfcSpecific,
                       lfcshared = input$lfcShared,
                       pvaltypespecific = input$pvalueSpecific,
                       pvalshared = input$pvalueShared
                     )
    return(categAltreOut)
  })


  comparePeaksObj <- eventReactive(input$buttoncompare, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0.2, detail = "Comparing Methods")
                   comparePeaksOut <- comparePeaksAltre(
                     req(categAltreObj()))
                   setProgress(value = 1, detail = "Done!")
                   Sys.sleep(0.5)
                 })
    return(comparePeaksOut)
  })

  pathGREATObj <- eventReactive(input$buttongreat, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0.2, detail = "Running GREAT")
                   runGREATcalls <- runGREAT(req(categAltreObj()))
                   runGREATout <- processPathways(req(runGREATcalls))
                   setProgress(value = 1, detail = "Done!")
                   Sys.sleep(0.5)
                 })
    return(runGREATout)
  })
  ############################################################################
   #  get input

  output$chooseref <- renderUI({
    reflist <- unique(loadCSVObj()$sample)
    selectInput("reference",
                "Select which cell-type to act as reference",
                reflist ,
                selected = reflist[1])
  })


  output$chooseChrom <- renderUI({
    peaks <- req(combineAnnotateObj())
    chroChoices <- unique(as.character(GenomeInfoDb::seqnames(peaks[[1]])))
    selectInput("chr", "Choose Chromosome", chroChoices, selected = "chr21")

  })

  ### palettes
  colrs <-   rownames(RColorBrewer::brewer.pal.info)

  output$choosePalette1 <- renderUI({
    selectInput("palette1",
                "Select color palette",
                colrs ,
                selected = colrs[15])
  })

  output$choosePalette2 <- renderUI({
    selectInput("palette2",
                "Select color palette",
                colrs ,
                selected = colrs[15])
  })

  output$choosePalette3 <- renderUI({
    selectInput("palette3",
                "Select color palette",
                colrs ,
                selected = colrs[15])
  })

  output$choosePalette4 <- renderUI({
    selectInput("palette4",
                "Select color palette",
                colrs ,
                selected = colrs[15])
  })

  output$choosePalette5 <- renderUI({
    selectInput("palette5",
                "Select color palette",
                colrs ,
                selected = colrs[15])
  })

  ############################################################################
  # download buttons
  output$downloadAnnotate <- downloadHandler(
    filename = "annotatedRegions.csv",
    content = function(con) {
      writeConsensusRE(req(combineAnnotateObj()), con)
    }
  )
  output$downloadBED <- downloadHandler(
    filename = "AnnotatedTrack.bed",
    content = function(con) {
      writeBedFile(req(categAltreObj()), con)
    }
  )

  output$downloadCompareDT <- downloadHandler(
    filename = "dataTableRE.csv",
    content = function(con) {
      writeCompareRE(req(comparePeaksObj()), con)
    }
  )

  output$downloadGREAT <- downloadHandler(
    filename = "GREATpathways.zip",
    content = function(con) {
      writeGREATpath(req(pathGREATObj()), con)
    },
    contentType = "application/zip"
  )

  ####
  observeEvent(input$buttonstop, {
    stopApp(returnValue = invisible())
  })


  ############################################################################
  #tables
  output$table1 <- renderDataTable({
    if (!is.null(input$file)) {
     csvoutput <- loadCSVObj()[ , !(names(loadCSVObj()) %in% "datapath")]
     if (is.null(csvoutput)) {
	message("Check the format of the CSV file")
        data.frame(ERROR = "Check the format of the CSV file")
     }
     else {
	csvoutput[ , !(names(csvoutput) %in% "datapath")]
    }
   }
  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table2 <- renderDataTable({
    req(getConsensusObj())$consPeaksStats
  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table3 <- renderDataTable({
    req(combineAnnotateObj())$mergestats
  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table4 <- renderDataTable({
    req(categAltreObj())$stats
  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table5 <- renderDataTable({
    req(comparePeaksObj())$compareresults
  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table6 <- renderDataTable({
	mydf <- req(pathGREATObj())
   data.frame(Pathway = mydf$ExperimentSpecificByIntensity$stats[ , 1],
              Experiment_Specific = mydf$ExperimentSpecificByIntensity$stats[ , 2],
              Shared = mydf$Shared$stats[ , 2],
              Reference_Specific = mydf$ReferenceSpecificByIntensity$stats[ , 2])
 }, options = list(searching = FALSE,
                    paging = FALSE))

  ############################################################################
  # plots

  output$barplot <- highcharter::renderHighchart({
       plotConsensusPeaks(req(getConsensusObj()),
                       palette = input$palette1,
                       viewer = FALSE,
                       maintitle = input$consPlotTitle,
                       ylabel = input$consPlotylabel)
  })

  output$annotatebarplot <- renderUI({
    plotCombineAnnotatePeaks(req(combineAnnotateObj()),
                             viewer = FALSE,
                             palette = input$palette2,
                             leftmaintitle = input$combLeftPlotTitle,
                             rightmaintitle = input$combRightPlotTitle,
                             leftylabel = input$combLeftylabel,
                             rightylabel = input$combRightylabel)
  })

  output$densityplot <- highcharter::renderHighchart({
    plotGetCounts(req(getCountsObj()),
                  viewer = FALSE,
                  palette = input$palette3,
                  xlabel = input$countsxlabel,
                  ylabel = input$countsylabel,
                  maintitle = input$countsPlotTitle)
  })

  output$volcano <- renderUI({
                   plotCountAnalysis(req(categAltreObj()),
                                    viewer = FALSE,
                                   palette = input$palette4)
  })

  output$boxplotCounts <- highcharter::renderHighchart({
    plotDistCountAnalysis(req(categAltreObj()),
                          req(getCountsObj()),
                          viewer = FALSE,
                          palette = input$palette4)
  })

  output$pieplot <- renderUI({
    plotCompareMethodsAll(req(comparePeaksObj()),
                          viewer = FALSE,
                          palette = input$palette4,
                          title11 = input$title11,
                          title12 = input$title12,
                          title13 = input$title13,
                          title21 = input$title21,
                          title22 = input$title22,
                          title23 = input$title23)
  })

  output$heatplotGREAT <- highcharter::renderHighchart({
    plotGREATenrich(req(pathGREATObj()),
                    maintitle = "GO Molecular Function",
                    pathwaycateg = "GO_Molecular_Function")
  })


  ############################################################################
  # info boxes

  output$statusbox1 <- renderInfoBox({
    if (is.null(input$file)) {
      infoBox(
        "Status",
        "File Not Loaded Yet!",
        icon = icon("import", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE
      )}
    else if (!is.null(input$file)) {
    infoBox(
      "Status",
      HTML(paste("File Loading Complete.",
                 "You Can Proceed to Step 2.",
                 sep = "<br/>")),
      icon = icon("thumbs-up", lib = "glyphicon"),
      color = "green", fill = TRUE)
    }
  })

  output$statusbox2 <- renderInfoBox ({
   if (input$buttonmerge == 0) {
      infoBox(
        "Status",
        "Merge Button Not Clicked Yet!",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE)
     }
    else if (input$buttonmerge > 0 && is.null(input$file)) {
      infoBox(
        "Status",
        HTML(paste("Step 2 is Not Complete Yet!",
                   "Please Run Step 1 Before Proceeding!.",
                   sep = "<br/>")),
        icon = icon("warning-sign", lib = "glyphicon"),
        color = "red",
        fill = TRUE)
    }
    else if (input$buttonmerge > 0 && !is.null(getConsensusObj())) {
      infoBox(
        "Status",
        HTML(paste("Replicates Have Been Merged.",
                   "You Can Proceed to Step 3.",
                   sep = "<br/>")),
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "green",
        fill = TRUE)
    }
  })

  output$statusbox3 <- renderInfoBox({
    if (input$buttonannot == 0) {
      infoBox(
        "Status",
        "Annotate Button Not Clicked Yet!",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE)
    }
    else if (input$buttonannot > 0 && (input$buttonmerge == 0 || is.null(input$file))) {
      infoBox(
        "Status",
        HTML(paste("Step 2 is Not Complete Yet!",
                   "Please Run Previous Steps Before Proceeding!",
                   sep = "<br/>")),
        icon = icon("warning-sign", lib = "glyphicon"),
        color = "red",
        fill = TRUE)
    }
    else if (input$buttonannot > 0 && !is.null(combineAnnotateObj())) {
      infoBox(
        "Status",
        HTML(paste("Peaks Have Been Annotated.",
                   "You Can Proceed to Step 4.",
                   "(If You Change the Parameters,",
                    "Please Press Button Again.)",
                   sep = "<br/>")),
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "green",
        fill = TRUE)
    }
  })

  output$statusbox4 <- renderInfoBox({
    if (input$buttoncounts == 0) {
      infoBox(
        "Status",
        "Retrieve Counts Button Not Clicked Yet!",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE)
    }
    else if (input$buttoncounts > 0 && (input$buttonannot == 0 ||
                                        input$buttonmerge == 0 ||
                                         is.null(input$file))) {
      infoBox(
        "Status",
        HTML(paste("Step 3 is Not Complete Yet!",
                   "Please Run Previous Steps Before Proceeding!",
                   sep = "<br/>")),
        icon = icon("warning-sign", lib = "glyphicon"),
        color = "red",
        fill = TRUE)
    }
    else if (input$buttoncounts > 0 && !is.null(getCountsObj())) {
      infoBox(
        "Status",
        HTML(paste("Counts Have Been Retrieved.",
                   "You Can Proceed to Step 5.",
                   sep = "<br/>")),
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "green",
        fill = TRUE)
    }
  })

  output$statusbox5 <- renderInfoBox({
    if (input$buttondefine == 0) {
      infoBox(
        "Status",
        "Identify Altered Regions Button Not Clicked Yet!",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE
      )}
    else if (input$buttondefine > 0 &&  (input$buttoncounts == 0 ||
                                         input$buttonannot == 0 ||
                                         input$buttonmerge == 0 ||
                                         is.null(input$file))) {
      infoBox(
        "Status",
        HTML(paste("Step 4 is Not Complete Yet!",
                   "Please Run Previous Steps Before Proceeding!",
                   sep = "<br/>")),
        icon = icon("warning-sign", lib = "glyphicon"),
        color = "red",
        fill = TRUE)
    }
    else if (input$buttondefine > 0 && !is.null(getAlteredObj())) {
      infoBox(
        "Status",
        HTML(paste("Altered Regions Have Been Identified.",
                   "You Can Proceed to Step 6.",
                   sep = "<br/>")),
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "green",
        fill = TRUE
      )
    }
  })

  output$statusbox6 <- renderInfoBox({
    if (input$buttoncat == 0) {
      infoBox(
        "Status",
        "Categorize Regions Button Not Clicked Yet!",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE
      )}
    else if (input$buttoncat > 0 && (input$buttondefine == 0 ||
                                     input$buttoncounts == 0 ||
                                     input$buttonannot  == 0 ||
                                     input$buttonmerge  == 0 ||
                                     is.null(input$file))) {
      infoBox(
        "Status",
        HTML(paste("Step 5 is Not Complete Yet!",
                   "Please Run Previous Steps Before Proceeding!",
                   sep = "<br/>")),
        icon = icon("warning-sign", lib = "glyphicon"),
        color = "red",
        fill = TRUE)
    }
    else if (input$buttoncat > 0 && !is.null(categAltreObj())) {
      infoBox(
        "Status",
        HTML(paste("Altered Regions Have Been Categorized.",
                   "You Can Proceed to Step 7.",
                   sep = "<br/>")),
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "green",
        fill = TRUE
      )
    }
  })

  output$statusbox7 <- renderInfoBox({
    if (input$buttoncompare == 0) {
      infoBox(
        "Status",
        "Compare Methods Button Not Clicked Yet!",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE)
    }
    else if (input$buttoncompare > 0 && (input$buttoncat    == 0 ||
                                           input$buttondefine == 0 ||
                                           input$buttoncounts == 0 ||
                                           input$buttonannot  == 0 ||
                                           input$buttonmerge  == 0 ||
                                           is.null(input$file))) {
      infoBox(
        "Status",
        HTML(paste("Step 6 is Not Complete Yet!",
                   "Please Run Previous Steps Before Proceeding!",
                   sep = "<br/>")),
        icon = icon("warning-sign", lib = "glyphicon"),
        color = "red",
        fill = TRUE)
    }
    else if (input$buttoncompare > 0 && !is.null(comparePeaksObj())) {
      infoBox(
        "Status",
        "Method Comparison Completed.",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "green",
        fill = TRUE)
    }
  })

  output$statusbox9 <- renderInfoBox({
    if (input$buttongreat == 0) {
      infoBox(
        "Status",
        "Run GREAT Button Not Clicked Yet!",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE)
    }
    else if (input$buttongreat > 0 && (input$buttoncat    == 0 ||
                                        input$buttondefine == 0 ||
                                        input$buttoncounts == 0 ||
                                        input$buttonannot  == 0 ||
                                        input$buttonmerge  == 0 ||
                                        is.null(input$file))) {
      infoBox(
        "Status",
        HTML(paste("Step 6 is Not Complete Yet!",
                   "Please Run Previous Steps Before Proceeding!",
                   sep = "<br/>")),
        icon = icon("warning-sign", lib = "glyphicon"),
        color = "red",
        fill = TRUE)
    }
    else if (input$buttongreat > 0 && !is.null(pathGREATObj())) {
      infoBox(
        "Status",
        "RUN GREAT Completed.",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "green",
        fill = TRUE)
    }
  })



  ##########################################################
  })

