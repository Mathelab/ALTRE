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
      csvFile()[,-1]
    }
  }, options = list(searching = FALSE,
                    paging = FALSE))

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
                     combineAnnotatePeaks(conspeaks = consPeaks,
                                          TSS = TSSannot)
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
                   csvfile <- csvFile()
                   setProgress(value = 0.5, detail = "annotating peaks")
                   counts_consPeaks <-
                     getcounts(
                       annotpeaks = consPeaksAnnotated,
                       csvfile = csvfile,
                       reference = "SAEC",
                       chrom = "chr21"
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
                     countanalysis(counts = counts_consPeaks,
                                   pval = 0.01,
                                   lfcvalue = 1)
                   setProgress(value = 1, detail = "done!")
                   Sys.sleep(0.5)
                 })

    return(altred_peaks)
  })


  # peaksdf <- reactive({
  #
  #   mergedpeaks <-cbind(PeaksFile = c("Consensus","Rep I","Rep II"),
  #         as.data.frame(mergedPeaks()$consPeaksStats))
  #   return(mergedpeaks)
  #
  # })

  output$table2 <- renderDataTable({
    mergedPeaks()$consPeaksStats

  }, options = list(searching = FALSE,
                    paging = FALSE))

  output$table3 <- renderDataTable({
    as.data.frame(annotatePeaks()$consPeaksAnnotated)

  }, options = list(searching = FALSE))


  ####################################
  # plots

  output$barplot <- renderPlot({
    plotConsensusPeaks(mergedPeaks()$consPeaksStats)
  })

  output$densityplot <- renderPlot({
    plotgetcounts(countsPeaks())
  })

  output$volcanoplot <- renderPlot({
    plotCountAnalysis(alteredPeaks())
  })


})
