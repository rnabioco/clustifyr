#' Launch Shiny app version of clustifyr, 
#' may need to run install_clustifyr_app() at first time to install packages
#' @return instance of shiny app
#' @examples
#' \dontrun{
#' run_clustifyr_app()
#' }
#' @export
run_clustifyr_app <- function() {
  appDir <- system.file("shinyapp", package = "clustifyr")
  shiny::runApp(appDir, display.mode = "normal")
}
