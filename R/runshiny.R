#' run shiny app
#' @export
runShinyApp <- function() {
    appDir <- system.file("shinyApp", package = "ALTRE")
    if (appDir == "") {
        stop(" The ShinyApp directory was not found.
             Try re-installing `ALTRE`.",
             call. = FALSE)
    }
    shiny::runApp(appDir, display.mode = "normal")
}
