
shinyServer(function(input, output, session) {

  ############################################################################
  # Load data

  loadCSVObj <- reactive({
      loadCSVFile(req(input$file)$datapath)
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
                     regionspecific = input$regionradio,
                     mergedist = input$dist,
                     mergedistenh = input$distenh,
                     mergedistprom = input$distprom,
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
                   getCountsOut <- getCounts(
                     annotpeaks = req(combineAnnotateObj()),
                     sampleinfo = req(loadCSVObj()),
                     reference = input$reference,
                     chrom = input$chr)
                   setProgress(value = 1, detail = "done!")
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
                   setProgress(value = 1, detail = "done!")
                   Sys.sleep(0.5)
                 })
    return(altred)
  })

  observeEvent(input$buttondefine, {
    getAlteredObj()
  })



  categAltreObj <- eventReactive(input$buttoncat, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0.2, detail = "categorizing altered peaks")
                   categAltreOut <-
                     categAltrePeaks(
                       analysisresults = req(getAlteredObj()),
                       lfctypespecific = input$lfcSpecific,
                       lfcshared = input$lfcShared,
                       pvaltypespecific = input$pvalueSpecific,
                       pvalshared = input$pvalueShared
                     )
                   setProgress(value = 1, detail = "done!")
                   Sys.sleep(0.5)
                 })
    return(categAltreOut)
  })


  comparePeaksObj <- eventReactive(input$buttoncompare, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0.2, detail = "Comparing Methods")
                   comparePeaksOut <- comparePeaksAltre(
                     req(categAltreObj()),
                     reference = req(input$reference)
                     )
                   setProgress(value = 1, detail = "done!")
                   Sys.sleep(0.5)
                 })
    return(comparePeaksOut)
  })

  pathenrichMFObj <- eventReactive(input$buttonpathwayMF, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0.2, detail = "MF: GO Enrichment Analysis")
                   MFenrich <-
                       pathenrich(
                         analysisresults = req(categAltreObj()),
                         ontoltype = "MF",
                         enrichpvalfilt = input$pathpvaluecutoffMF
                       )
                   setProgress(value = 1, detail = "Done!")
                   Sys.sleep(0.5)
                 })
    return(MFenrich)
  })

  pathenrichBPObj <- eventReactive(input$buttonpathwayBP, {
    withProgress(message = 'In progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0.2, detail = "BP: GO Enrichment Analysis")
                   BPenrich <-
                     pathenrich(
                       analysisresults = req(categAltreObj()),
                       ontoltype = "BP",
                       enrichpvalfilt = input$pathpvaluecutoffBP
                     )
                   setProgress(value = 1, detail = "Done!")
                   Sys.sleep(0.5)
                 })
    return(BPenrich)
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

  output$downloadData <- downloadHandler(
    filename = "AnnotatedTrack.bed",
    content = function(con) {
      writeBedFile(req(categAltreObj()), con)
    }
  )

  ############################################################################
  #tables
  output$table1 <- renderDataTable({
    if (!is.null(input$file)) {
      loadCSVObj()[, -1]
    }
  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table2 <- renderDataTable({
    getConsensusObj()$consPeaksStats
  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table3 <- renderDataTable({
    combineAnnotateObj()$mergestats
  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table4 <- renderDataTable({
    req(categAltreObj()$stats)
  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table5 <- renderDataTable({
    comparePeaksObj()
  }, options = list(searching = FALSE,
                    paging = FALSE))

  ############################################################################
  # plots

  output$barplot <- renderPlot({
    plotConsensusPeaks(getConsensusObj())
  })

  output$annotatebarplot <- renderPlot({
    plotCombineAnnotatePeaks(combineAnnotateObj())
  })

  output$densityplot <- renderPlot({
    plotgetcounts(getCountsObj())
  })

  output$volcanoplot <- renderPlot({
    plotCountAnalysis(req(categAltreObj()))
  })

  output$boxplot <- renderPlot({
    plotDistCountAnalysis(req(categAltreObj()), req(getCountsObj()))
  })

  output$heatplotMF <- renderPlot({
      enrichHeatmap(req(pathenrichMFObj()), title = "GO:MF, p<0.01")
  })

  output$heatplotBP <- renderPlot({
        enrichHeatmap(req(pathenrichBPObj()), title = "GO:BP, p<0.01")
  })

  output$vennplot <- renderPlot({
    plotallvenn(req(comparePeaksObj()))
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
      "File Loading Complete. You Can Proceed to Step 2.",
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
    else if(input$buttonmerge > 0 && !is.null(getConsensusObj())) {
      infoBox(
        "Status",
        "Replicates Have Been Merged. You Can Proceed to Step 3.",
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
    else if (input$buttonannot > 0 && !is.null(combineAnnotateObj())) {
      infoBox(
        "Status",
        "Peaks Have Been Annotated (If You Change the Parameters,
        Please Press Button Again). You Can Proceed to Step 4.",
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
    else if (input$buttoncounts > 0 && !is.null(getCountsObj())) {
      infoBox(
        "Status",
        "Counts Have Been Retrieved. You Can Proceed to Step 5.",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "green",
        fill = TRUE)
    }
  })

  output$statusbox5 <- renderInfoBox({
    if (input$buttondefine == 0) {
      infoBox(
        "Status", "Define Altered Regions Button Not Clicked Yet!",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE
      )}
    else if (input$buttondefine > 0 && !is.null(getAlteredObj())) {
      infoBox(
        "Status", "Altered Regions Have Been Defined.
        You Can Proceed to Step 6.",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "green",
        fill = TRUE
      )
    }
  })

  output$statusbox5b <- renderInfoBox({
    if (input$buttoncat == 0) {
      infoBox(
        "Status", "Categorize Altered Regions Button Not Clicked Yet!",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE
      )}
    else if (input$buttoncat > 0 && !is.null(categAltreObj())) {
      infoBox(
        "Status", "Altered Regions Have Been Categorized.
        You Can Proceed to Step 7.",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "green",
        fill = TRUE
      )
    }
  })

  output$statusbox6 <- renderInfoBox({
    if (input$buttonpathwayMF == 0) {
      infoBox(
        "Status",
        "MF Enrichment Analysis Button Not Clicked Yet!",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE)
      }
    else if (input$buttonpathwayMF > 0 && !is.null(pathenrichMFObj())) {
      infoBox(
        "Status",
        "MF Enrichment Analysis Has Been Run.",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "green",
        fill = TRUE)
    }
  })

  output$statusbox7 <- renderInfoBox({
    if (input$buttonpathwayBP == 0) {
      infoBox(
        "Status",
        "BP Enrichment Analysis Button Not Clicked Yet!",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
        fill = TRUE)
      }
    else if(input$buttonpathwayBP > 0 && !is.null(pathenrichBPObj())) {
      infoBox(
        "Status",
        "BP Enrichment Analysis Completed.",
        icon = icon("thumbs-up", lib = "glyphicon"),
        color = "green",
        fill = TRUE
      )
    }
  })

  output$statusbox8 <- renderInfoBox({
    if (input$buttoncompare == 0) {
      infoBox(
        "Status",
        "Compare Methods Button Not Clicked Yet!",
        icon = icon("flag", lib = "glyphicon"),
        color = "aqua",
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
  ##########################################################
  })

