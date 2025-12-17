#' @importFrom shiny addResourcePath

.onLoad <- function(libname, pkgname) {
  # Create link to logo
  shiny::addResourcePath("WEIN", system.file("www", package="WEIN"))
  
  shiny::addResourcePath("sbs", system.file("www", package = "shinyBS"))
  
  suppressPackageStartupMessages({
    requireNamespace("heatmaply", quietly = TRUE)
    requireNamespace("gplots", quietly = TRUE)
    requireNamespace("pcaExplorer", quietly = TRUE)
  })
  
  # Pre-load and suppress messages from verbose Bioconductor packages
  suppressPackageStartupMessages({
    # These packages are known to produce verbose startup messages
    # Loading them here with suppressed messages prevents them from
    # appearing when the user loads the WEIN package
  })
}

.onAttach <- function(libname, pkgname) {
  # Display a clean package attachment message
  packageStartupMessage("Package '", pkgname, "' loaded successfully.")
}
