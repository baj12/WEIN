#' @importFrom shiny addResourcePath

.onLoad <- function(libname, pkgname) {
  # Create link to logo
  shiny::addResourcePath("ideal", system.file("www", package="WEIN"))
  
  shiny::addResourcePath("sbs", system.file("www", package = "shinyBS"))
}

