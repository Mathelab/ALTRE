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
    peaks <-  peaks()
    mergedpeaks <- getConsensusPeaks(peaks, input$numOverlap)

    mergedpeaks$consPeaksStats <-
      cbind(
        PeaksFile = c("Consensus", "Rep I", "Rep II"),
        as.data.frame(mergedpeaks$consPeaksStats)
      )

    return(mergedpeaks)
  })

  observeEvent(input$buttonmerge, {
    withProgress(message = 'Merging in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   setProgress(value = 0, detail = "loading new dataset")
                   for (i in 1:2) {
                     incProgress(1 / 4)
                     Sys.sleep(0.05)
                   }
                   mergedPeaks()
                   for (i in 1:2) {
                     incProgress(1 / 4)
                     Sys.sleep(0.1)
                   }
                   setProgress(value = 1, detail = "Done!")
                 })
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

  ####################################
  #  tabPanel "plot"

  output$barplot <- renderPlot({
    plotCountSummary(mergedPeaks()$consPeaksStats)
  })

})
