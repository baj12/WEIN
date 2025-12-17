#' WEIN: Web-based Engine for Interactive Next-generation sequencing analysis
#'
#' WEIN makes bulk RNAseq data analysis interactive, easy and reproducible.
#' The analysis of RNA-seq datasets is guided by the Shiny app as main component of
#' the package, which also provides a wide set of functions to efficiently extract
#' information from the existing data. The app can be also deployed on a Shiny
#' server, to allow its usage without any installation on the user's side.
#'
#' @import DESeq2
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame List
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom shiny actionButton bookmarkButton conditionalPanel downloadButton fileInput
#' @importFrom shiny fluidRow h2 h3 h4 hr p selectInput tags textInput updateSelectInput
#' @importFrom shiny updateSelectizeInput updateTextInput verbatimTextOutput wellPanel withProgress
#' @importFrom DT datatable dataTableOutput renderDataTable
#' @import shinydashboard
#' @importFrom AnnotationDbi mapIds keytypes
#' @importFrom shinyAce aceAutocomplete aceEditor getAceModes getAceThemes
#' updateAceEditor
#' @import BiocParallel
#' @import knitr
#' @import rmarkdown
#' @importFrom dplyr inner_join tbl_df filter mutate arrange last select
#' @importMethodsFrom GOstats hyperGTest summary
#' @import GO.db
#' @importFrom UpSetR upset fromList
#' @importFrom goseq getgo goseq nullp
#' @import pcaExplorer
#' @importFrom gplots venn
#' @importFrom IHW ihw
#' @importFrom rentrez entrez_summary
#' @importFrom limma goana topGO
#' @import topGO
#' @importFrom heatmaply heatmaply
#' @importFrom stringr str_count
#' @importFrom rintrojs introjs introjsUI
#' @importFrom shinyBS bsTooltip bsCollapse bsCollapsePanel
#' @importFrom ggrepel geom_text_repel
#' @importFrom base64enc dataURI
#' @importFrom grDevices dev.off pdf axisTicks
#' @import methods
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plotly plotlyOutput ggplotly renderPlotly
#'
#' @author
#' Bernd Jagla \email{bernd.jagla@@pasteur.fr}, 2025
#'
#' Maintainer: Bernd Jagla \email{bernd.jagla@@pasteur.fr}
#' @name WEIN-pkg
#' @docType package
#' @keywords internal
"_PACKAGE"
