#' WEIN: Web-based Engine for Interactive Next-generation sequencing analysis
#'
#' WEIN makes bulk RNAseq data analysis interactive, easy and reproducible.
#' This function launches the main application included in the package.
#'
#' @param dds_obj A \code{\link{DESeqDataSet}} object. If not provided, then a
#' \code{countmatrix} and a \code{expdesign} need to be provided. If none of
#' the above is provided, it is possible to upload the data during the
#' execution of the Shiny App
#' @param res_obj  A \code{\link{DESeqResults}} object. If not provided, it can
#' be computed during the execution of the application
#' @param annotation_obj A \code{data.frame} object, with row.names as gene
#' identifiers (e.g. ENSEMBL ids) and a column, \code{gene_name}, containing
#' e.g. HGNC-based gene symbols. If not provided, it can be constructed during
#' the execution via the org.eg.XX.db packages - these need to be installed
#' @param countmatrix A count matrix, with genes as rows and samples as columns.
#' If not provided, it is possible to upload the data during the execution of
#' the Shiny App
#' @param expdesign A \code{data.frame} containing the info on the covariates
#' of each sample. If not provided, it is possible to upload the data during the
#' execution of the Shiny App
#' @param gene_signatures A list of vectors, one for each pathway/signature. This 
#' is for example the output of the \code{\link{read_gmt}} function. The provided
#' object can also be replaced during runtime in the dedicated upload widget.
#'
#' @return A Shiny App is launched for interactive data exploration and
#' differential expression analysis
#'
#' @export WEIN
#' 
#' @examples
#' # with simulated data...
#' library(DESeq2)
#' dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8)
#' cm <- counts(dds)
#' cd <- colData(dds)
#'
#' # with the well known airway package...
#' library(airway)
#' data(airway)
#' airway
#' dds_airway <- DESeq2::DESeqDataSetFromMatrix(assay(airway),
#'                                              colData = colData(airway),
#'                                              design=~cell+dex)
#' \dontrun{
#'
#' idealImmunoTP()
#' idealImmunoTP(dds)
#' idealImmunoTP(dds_airway)
#'
#' dds_airway <- DESeq2::DESeq(dds_airway)
#' res_airway <- DESeq2::results(dds_airway)
#' idealImmunoTP(dds_airway, res_airway)
#' }
#'
#' 

# Import necessary functions from other files
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar dashboardBody
#' @importFrom shiny actionButton bookmarkButton conditionalPanel downloadButton fileInput 
#' fluidRow h2 h3 h4 hr p selectInput tags textInput updateSelectInput updateSelectizeInput 
#' updateTextInput verbatimTextOutput wellPanel withProgress
#' @importFrom shinyAce aceEditor
#' @importFrom shinyBS bsCollapse bsCollapsePanel bsTooltip dropdownMenu menuItem 
#' notificationItem
#' @importFrom DT dataTableOutput renderDataTable
#' @importFrom plotly plotlyOutput renderPlotly
#' @importFrom shinyjqui jqui_resizable
#' @importFrom rintrojs introjs introjsUI
#' @importFrom base64enc dataURI
#' @importFrom utils read.delim readLines
#' @importFrom grDevices barplot colorRampPalette dev.off pdf
#' @importFrom graphics text
#' @importFrom stats as.formula dist formula median model.matrix na.omit quantile 
#' relevel reshape setNames terms var
#' @importFrom methods is new
#' @importFrom BiocParallel MulticoreParam
#' @importFrom S4Vectors mcols
#' @importFrom IRanges IRanges
#' @importFrom SummarizedExperiment assay assays colData rowData
#' @importFrom GenomicRanges GRanges
#' @importFrom DESeq2 DESeq DESeqDataSet DESeqDataSetFromMatrix estimateSizeFactors 
#' results sizeFactors vst counts design resultsNames normTransform
#' @importFrom pcaExplorer plotDispEsts
#' @importFrom ggplot2 aes_string element_blank element_text geom_bar geom_histogram 
#' geom_hline geom_jitter geom_point geom_raster geom_text geom_vline ggplot 
#' scale_colour_manual scale_fill_gradient scale_x_discrete scale_y_continuous 
#' scale_y_log10 theme theme_bw xlab ylab
#' @importFrom ggrepel geom_text_repel
#' @importFrom reshape2 melt
#' @importFrom cowplot plot_grid
#' @importFrom gplots venn
#' @importFrom UpSetR upset fromList
#' @importFrom limma goana topGO
#' @importFrom IHW ihw
#' @importFrom rentrez entrez_summary
#' @importFrom AnnotationDbi mapIds
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr arrange filter inner_join mutate select summarise
#' @importFrom tidyr separate_rows unite
#' @importFrom stringr str_replace
#' @importFrom knitr kable
#' @importFrom rmarkdown render
#' @importFrom BiocGenerics rowVars
#' @importFrom Biobase featureNames
#' @importFrom methods is
NULL

#' Multi-dimensional PCA plot
#'
#' Creates a PCA plot with customizable dimensions
#'
#' @param object A DESeqTransform object
#' @param intgroup Character vector of interesting groups
#' @param ntop Number of top genes to use
#' @param returnData Logical, whether to return data instead of plot
#' @param pc1 First principal component to plot
#' @param pc2 Second principal component to plot
#' @return Either a ggplot object or PCA data
#' @export
multiAxPCA = function (object, intgroup = "condition", ntop = 500, returnData = FALSE, pc1=1, pc2=2) {
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, pc1], PC2 = pca$x[, pc2], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    return(percentVar)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 3) + xlab(paste0("PC",pc1,": ", round(percentVar[pc1] * 
                                                              100), "% variance")) + ylab(paste0("PC",pc2,": ", round(percentVar[pc2] * 
                                                                                                                        100), "% variance")) + coord_fixed()
}

#' Main function to launch the idealImmunoTP Shiny application
#'
#' Launches the interactive differential expression analysis application
#'
#' @param dds_obj A DESeqDataSet object
#' @param res_obj A DESeqResults object
#' @param annotation_obj A data frame with gene annotations
#' @param countmatrix A count matrix
#' @param expdesign A data frame with experimental design
#' @param gene_signatures A list of gene signatures
#' @return A Shiny application
#' @export
WEIN<- function(dds_obj = NULL,
                         res_obj = NULL,
                         annotation_obj = NULL,
                         countmatrix = NULL,
                         expdesign = NULL,
                         gene_signatures = NULL,
                         cur_species = NULL,
                         cur_type = NULL){
  
  if ( !requireNamespace('shiny',quietly = TRUE) ) {
    stop("idealImmunoTP requires 'shiny'. Please install it using
         install.packages('shiny')")
  }
  if ( !requireNamespace('tidyr',quietly = TRUE) ) {
    stop("idealImmunoTP requires 'tidyr'. Please install it using
         install.packages('tidyr')")
  }
  require(tidyr)
  require(shiny)
  require(shinydashboard)
  # create environment for storing inputs and values
  ## i need the assignment like this to export it up one level - i.e. "globally"
  ideal_env <<- new.env(parent = emptyenv())
  
  ## upload max 300mb files - can be changed if necessary
  options(shiny.maxRequestSize=600*1024^2)
  options(shiny.launch.browser = TRUE)
  app_data <- list(
    dds_obj = dds_obj,
    res_obj = res_obj,
    annotation_obj = annotation_obj,
    countmatrix = countmatrix,
    expdesign = expdesign,
    gene_signatures = gene_signatures,
    cur_species = cur_species,
    cur_type = cur_type
  )
  # Launch the app!
  shinyApp(ui = idealImmunoTP_ui, 
           server = idealImmunoTP_server( app_data) 
             )
}