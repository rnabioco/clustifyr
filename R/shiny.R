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

#' Install all packages needed for shiny app
#' @importFrom utils install.packages installed.packages read.csv
#' @return installation complete or error messages
#' @examples
#' \dontrun{
#' install_clustifyr_app()
#' }
#' @export
install_clustifyr_app <- function() {
  options(repos=structure(c(CRAN="http://cloud.r-project.org/")))
  pack <- read.csv(system.file("shinyapp", "data", "dependencies.csv", package = "clustifyr"))
  inst <- installed.packages()[,"Package"]
  pack_need <- !unlist(vapply(pack$package, function(x) {x %in% inst}, logical(1)))
  pack2 <- pack[pack_need, ]
  if (!("BiocManager" %in% inst)) {
    install.packages("BiocManager")
  }
  if (!("remotes" %in% inst)) {
    install.packages("remotes")
  }
  while (nrow(pack2) != 0) {
    message("Installing ", nrow(pack2), " package(s)...")
    if (pack2$source[1] == "CRAN") {
      install.packages(pack2$package[1], dependencies = TRUE, ask = FALSE)
    } else if (pack2$source[1] == "Bioconductor") {
      BiocManager::install(pack2$package[1], ask = FALSE, update = FALSE)
    } else {
      message("Installing latest GitHub devel version of clustifyr")
      remotes::install_github("rnabioco/clustifyr", dependencies = TRUE, upgrade = FALSE)
    }
    inst <- installed.packages()[,"Package"]
    pack_need <- !unlist(vapply(pack2$package, function(x) {x %in% inst}, logical(1)))
    pack2 <- pack2[pack_need, ]
  }
  message("All dependent packages installed")
}
