
shinyServer(function(input, output, session) {

  ############################################################################
  # Load data

  csvFile <- reactive({
      loadCSVFile(req(input$file)$datapath)
  })

  peaks <- reactive({
      loadBedFiles(req(csvFile()))
  })

  ############################################################################
  # function calls

  mergedPeaks <- eventReactive(input$buttonmerge, {
    withProgress(message = 'In Progress:',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = .1, detail = "Loading Peak Files")
                   peaks <-  req(peaks())
                   setProgress(value = 0.5, detail = "Merging Replicates")
                   consensusPeaks <- getConsensusPeaks(peaks, req(input$numOverlap))
                   setProgress(value = 1, detail = "Done!")
                   Sys.sleep(0.5)
                 })
    return(consensusPeaks)
  })

  annotatePeaks <- eventReactive(input$buttonannot, {
    withProgress(message = 'In Progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = .1, detail = "Retrieving TSS File")
                   TSSannot <- getTSS()
                   setProgress(value = 0.2, detail = "Annotating Peaks")
                   annotatedPeaks <- combineAnnotatePeaks(
                     conspeaks = mergedPeaks(),
                     TSS = TSSannot,
                     merge = input$mergeradio,
                     regionspecific = input$regionradio,
                     mergedist = input$dist,
                     mergedistenh = input$distenh,
                     mergedistprom = input$distprom,
                     distancefromTSS = input$distTSS)
                   setProgress(value = 1, detail = "Done!")
                   Sys.sleep(0.5)
                 })
    return(annotatedPeaks)
  })

  countsPeaks <- eventReactive(input$buttoncounts, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0.1, detail = "Retrieving counts")
                   countsSummary <- getcounts(
                     annotpeaks = annotatePeaks(),
                     sampleinfo = csvFile(),
                     reference = req(input$reference))
                   setProgress(value = 1, detail = "done!")
                   Sys.sleep(0.5)
                 })
    return(countsSummary)
  })

  alteredPeaks <- eventReactive(input$buttondefine, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0.2, detail = "annotating peaks")
                   altred <-
                     countanalysis(
                       counts = countsPeaks(),
                       pval = req(input$alpha),
                       lfcvalue = req(input$lfcThreshold)
                     )
                   setProgress(value = 1, detail = "done!")
                   Sys.sleep(0.5)
                 })
    return(altred)
  })


  compareMethods <- eventReactive(input$buttoncompare, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0.2, detail = "Comparing Methods")
                   compareResults <- resultsComparison(alteredPeaks(),
                                                        reference = "A549")
                   setProgress(value = 1, detail = "Done!")
                   Sys.sleep(0.5)
                 })
    return(compareResults)
  })

  pathewayOutputMF <- eventReactive(input$buttonpathwayMF, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0.2, detail = "MF: GO Enrichment Analysis")
                   MFenrich <-
                     pathenrich(
                       analysisresults = alteredPeaks(),
                       ontoltype = "MF",
                       enrichpvalfilt = req(input$pathpvaluecutoffMF)
                     )
                   setProgress(value = 1, detail = "Done!")
                   Sys.sleep(0.5)
                 })
    return(MFenrich)
  })

  pathewayOutputBP <- eventReactive(input$buttonpathwayBP, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0.2, detail = "BP: GO Enrichment Analysis")
                   BPenrich <-
                     pathenrich(
                       analysisresults = alteredPeaks(),
                       ontoltype = "BP",
                       enrichpvalfilt = req(input$pathpvaluecutoffBP)
                     )
                   setProgress(value = 1, detail = "Done!")
                   Sys.sleep(0.5)
                 })
    return(BPenrich)
  })
  ############################################################################
   #  get input

  output$chooseref <- renderUI({
    reflist <- unique(csvFile()$sample)
    selectInput("reference",
                "Reference Cell Type",
                reflist ,
                selected = reflist[1])
  })
  ############################################################################
  #tables
  output$table1 <- renderDataTable({
    if (!is.null(input$file)) {
      csvFile()[, -1]
    }
  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table2 <- renderDataTable({
    mergedPeaks()$consPeaksStats
  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table3 <- renderDataTable({
    annotatePeaks()$mergestats
  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table4 <- renderDataTable({
    compareMethods()
  }, options = list(searching = FALSE,
                    paging = FALSE))

  ############################################################################
  # plots

  output$barplot <- renderPlot({
    plotConsensusPeaks(mergedPeaks())
  })

  output$annotatebarplot <- renderPlot({
    plotCombineAnnotatePeaks(annotatePeaks())
  })

  output$densityplot <- renderPlot({
    plotgetcounts(countsPeaks())
  })

  output$volcanoplot <- renderPlot({
    plotCountAnalysis(alteredPeaks())
  })

  output$heatplotMF <- renderPlot({
      enrichHeatmap(pathewayOutputMF(), title = "GO:MF, p<0.01")
  })

  output$heatplotBP <- renderPlot({
        enrichHeatmap(pathewayOutputBP(), title = "GO:BP, p<0.01")
  })

  output$vennplot <- renderPlot({
    plotallvenn(compareMethods())
  })

  ############################################################################
  # info boxes

  output$statusbox1 <- renderInfoBox({
    if (is.null(input$file)) {
      infoBox(
        "Status",
        "File Not Loaded",
        icon = icon("import", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE
      )}
    else if (!is.null(input$file)) {
    infoBox(
      "Status",
      "File Loaded",
      icon = icon("thumbs-up", lib = "glyphicon"),
      color = "yellow", fill = TRUE)
    }
  })

  output$statusbox2 <- renderInfoBox ({
    if (input$buttonmerge == 0) {
      infoBox(
        "Status",
        "Merge Button Not Clicked Yet",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE)
      }
    else if(input$buttonmerge > 0) {
      infoBox(
        "Status",
        "Replicates Have Been Merged ",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow",
        fill = TRUE)
    }
  })

  output$statusbox3 <- renderInfoBox({
    if (input$buttonannot == 0) {
      infoBox(
        "Status",
        "Annotate Button Not Clicked Yet",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE)
      }
    else if (input$buttonannot > 0) {
      infoBox(
        "Status",
        "Peaks Have Been Annotated",
        icon = icon("thumbs-up",lib = "glyphicon"),
        color = "yellow",
        fill = TRUE)
    }
  })

  output$statusbox4 <- renderInfoBox({
    if (input$buttoncounts == 0){
      infoBox(
        "Status",
        "Retrieve Counts Button Not Clicked Yet",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE)
      }
    else if (input$buttoncounts > 0) {
      infoBox(
        "Status",
        "Counts Have Been Retrieved",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow",
        fill = TRUE)
    }
  })

  output$statusbox5 <- renderInfoBox({
    if (input$buttondefine == 0) {
      infoBox(
        "Status", "Define Altered Regions Button Not Clicked Yet",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE
      )}
    else if (input$buttondefine > 0) {
      infoBox(
        "Status", "Altered Regions Have Been Defined",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow",
        fill = TRUE
      )
    }
  })

  output$statusbox6 <- renderInfoBox({
    if (input$buttonpathwayMF == 0){
      infoBox(
        "Status",
        "Enrichment Analysis Button Not Clicked Yet",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE)
      }
    else if (input$buttonpathwayMF > 0) {
      infoBox(
        "Status",
        "Enrichment Analysis Has Been Run",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow",
        fill = TRUE)
    }
  })

  output$statusbox7 <- renderInfoBox({
    if (input$buttonpathwayBP == 0) {
      infoBox(
        "Status",
        "Enrichment Analysis BP Button Not Clicked Yet",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE)
      }
    else if(input$buttonpathwayBP > 0) {
      infoBox(
        "Status",
        "BP Enrichment Analysis Completed",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow",
        fill = TRUE
      )
    }
  })

  output$statusbox8 <- renderInfoBox({
    if (input$buttoncompare == 0) {
      infoBox(
        "Status",
        "Compare Methods Button Not Clicked Yet",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE)
      }
    else if (input$buttoncompare > 0) {
      infoBox(
        "Status",
        "Method Comparison Completed",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow",
        fill = TRUE)
    }
  })
  ##########################################################
  })

