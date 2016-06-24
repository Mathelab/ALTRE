# utility functions
################################################################

shinyServer(function(input, output, session) {
  ## tabName"definerep"

  #  tabPanel "Load CSV"

  inputFilePath <- reactive({
    if (!is.null(input$file)) {
      # User has not uploaded a file
      return(input$file)
    }
  })

  csvFile <- reactive({
    if (!is.null(input$file)) {
      inputPathObj <- inputFilePath()
      outputCSV <- loadCSVFile(inputPathObj$datapath)
      return(outputCSV)
    }
  })

  output$table1 <- renderDataTable({
    if (!is.null(input$file)) {
      csvFile()[, -1]
    }
  }, options = list(searching = FALSE,
                    paging = FALSE))


  output$chooseref <- renderUI({
    reflist <- unique(csvFile()$sample)
    selectInput("reference",
                "Reference Cell Type",
                reflist ,
                selected = reflist[1])
  })

#####################################################

  #  tabPanel "merge"

  peaks <- reactive({
    if (!is.null(input$file)) {
      fileIn <- csvFile()
      outputBed <- loadBedFiles(fileIn)
      return(outputBed)
    }
  })

  mergedPeaks <- eventReactive(input$buttonmerge, {
    withProgress(message = 'In progress:',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = .1, detail = "loading peak files")
                   peaks <-  peaks()
                   setProgress(value = 0.5, detail = "merging replicates")
                   consenpeaks <-
                     getConsensusPeaks(peaks, input$numOverlap)
                   setProgress(value = 1, detail = "done!")
                   Sys.sleep(0.5)
                 })
    return(consenpeaks)
  })

  annotatePeaks <- eventReactive(input$buttonannot, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = .1, detail = "running")
                   consPeaks <- mergedPeaks()
                   setProgress(value = .3, detail = "retrieving TSS file")
                   TSSannot <- getTSS()
                   setProgress(value = 0.5, detail = "annotating peaks")

                   consPeaksAnnotated <-
                     combineAnnotatePeaks(
                       conspeaks = consPeaks,
                       TSS = TSSannot,
                       merge = input$mergeradio,
                       regionspecific = input$regionradio,
                       mergedist = input$dist,
                       mergedistenh = input$distenh,
                       mergedistprom = input$distprom,
                       distancefromTSS = input$distTSS
                     )
                   setProgress(value = 1, detail = "done!")
                   Sys.sleep(0.5)
                 })

    return(consPeaksAnnotated)
  })


  countsPeaks <- eventReactive(input$buttoncounts, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = .1, detail = "running")
                   consPeaksAnnotated <- annotatePeaks()
                   #csvfile <- csvFile()
                   csvfile <- inputFilePath()$datapath
                   setProgress(value = 0.5, detail = "Retrieving counts")
                   counts_consPeaks <-
                     getcounts(
                       annotpeaks = consPeaksAnnotated,
                       csvfile = csvfile,
                       reference = input$reference
                     )
                   setProgress(value = 1, detail = "done!")
                   Sys.sleep(0.5)
                 })
    return(counts_consPeaks)
  })


  alteredPeaks <- eventReactive(input$buttondefine, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = .1, detail = "running")
                   counts_consPeaks <- countsPeaks()
                   setProgress(value = 0.5, detail = "annotating peaks")
                   altred_peaks <-
                     countanalysis(
                       counts = counts_consPeaks,
                       pval = input$alpha,
                       lfcvalue = input$lfcThreshold
                     )
                   setProgress(value = 1, detail = "done!")
                   Sys.sleep(0.5)
                 })

    return(altred_peaks)
  })

  pathewayOutputMF <- eventReactive(input$buttonpathwayMF, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = .1, detail = "running")
                   altre_peaks <- alteredPeaks()
                   setProgress(value = 0.3, detail = "GO Enrichment Analysis MF")
                   MFenrich <-
                     pathenrich(
                       analysisresults = altre_peaks,
                       ontoltype = "MF",
                       enrichpvalfilt = input$pathpvaluecutoffMF
                     )
                   setProgress(value = 1, detail = "done!")
                   Sys.sleep(0.5)
                 })

    return(MFenrich)
  })



  pathewayOutputBP <- eventReactive(input$buttonpathwayBP, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = .1, detail = "running")
                   altre_peaks <- alteredPeaks()
                   setProgress(value = 0.3, detail = "GO Enrichment Analysis BP")
                   BPenrich <-
                     pathenrich(
                       analysisresults = altre_peaks,
                       ontoltype = "BP",
                       enrichpvalfilt = input$pathpvaluecutoffBP
                     )
                   setProgress(value = 1, detail = "done!")
                   Sys.sleep(0.5)
                 })

    return(BPenrich)
  })

  output$table2 <- renderDataTable({
    mergedPeaks()$consPeaksStats

  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table3 <- renderDataTable({
    annotatePeaks()$mergestats
  }, options = list(searching = FALSE,
                    paging = FALSE))


  ####################################
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



  output$statusbox1 <- renderInfoBox({
    if(is.null(input$file)){
      infoBox(
        "Status", "File Not Loaded", icon = icon("import", lib = "glyphicon"),
        color = "aqua", fill = TRUE
      )}
    else if(!is.null(input$file)) {
    infoBox(
      "Status", "File Loaded", icon = icon("thumbs-up", lib = "glyphicon"),
      color = "yellow", fill = TRUE
    )
    }
  })


  output$statusbox2 <- renderInfoBox({
    if(input$buttonmerge==0){
      infoBox(
        "Status", "Merge Button Not Clicked Yet", icon = icon("flag", lib = "glyphicon"),
        color = "aqua", fill = TRUE
      )}
    else if(input$buttonmerge > 0) {
      infoBox(
        "Status", "Replicates Have Been Merged ", icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow", fill = TRUE
      )
    }
  })

  output$statusbox3 <- renderInfoBox({
    if(input$buttonannot==0){
      infoBox(
        "Status", "Annotate Button Not Clicked Yet",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua", fill = TRUE
      )}
    else if(input$buttonannot > 0) {
      infoBox(
        "Status", "Peaks Have Been Annotated",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow", fill = TRUE
      )
    }
  })

  output$statusbox4 <- renderInfoBox({
    if(input$buttoncounts==0){
      infoBox(
        "Status", "Retrieve Counts Button Not Clicked Yet",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua", fill = TRUE
      )}
    else if(input$buttoncounts > 0) {
      infoBox(
        "Status", "Counts Have Been Retrieved",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow", fill = TRUE
      )
    }
  })

  output$statusbox5 <- renderInfoBox({
    if(input$buttondefine==0){
      infoBox(
        "Status", "Define Altered Regions Button Not Clicked Yet",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua", fill = TRUE
      )}
    else if(input$buttondefine> 0) {
      infoBox(
        "Status", "Altered Regions Have Been Defined",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow", fill = TRUE
      )
    }
  })

  output$statusbox6 <- renderInfoBox({
    if(input$buttonpathwayMF==0){
      infoBox(
        "Status", "Enrichment Analysis Button Not Clicked Yet",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua", fill = TRUE
      )}
    else if(input$buttonpathwayMF > 0) {
      infoBox(
        "Status", "Enrichment Analysis Has Been Run",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow", fill = TRUE
      )
    }
  })

  output$statusbox7 <- renderInfoBox({
    if(input$buttonpathwayBP==0){
      infoBox(
        "Status", "Enrichment Analysis BP Button Not Clicked Yet",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua", fill = TRUE
      )}
    else if(input$buttonpathwayBP > 0) {
      infoBox(
        "Status", "BP Enrichment Analysis Completed",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "yellow", fill = TRUE
      )
    }
  })

  })

