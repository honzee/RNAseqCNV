#' launches the shinyAppDemo app
#'
#' @export launchApp
#'
#' @return shiny application object
#'


# wrapper for shiny::shinyApp()
launchApp <- function() {
  shinyApp(ui = shinyAppUi, server = shinyAppServer)
}
