#' WEIN Utilities
#'
#' Pure functions that can be used independently of the Shiny app for testing purposes.
#'
#' @docType package
#' @name utilities
#' @keywords internal
"_PACKAGE"

#' Read count matrix from file
#'
#' Reads a count matrix from a file, guessing the separator automatically.
#'
#' @param filepath Path to the file to read
#' @param sepguesser_function Function to use for guessing the separator (defaults to sepguesser)
#'
#' @return A matrix with count data
#' @export
#'
read_countmatrix <- function(filepath, sepguesser_function = sepguesser) {
  guessed_sep <- sepguesser_function(filepath)
  cm <- utils::read.delim(filepath, header = TRUE,
                          as.is = TRUE, sep = guessed_sep, quote = "",
                          row.names = 1,
                          check.names = FALSE)
  cm[is.na(cm)] <- 0
  return(cm)
}

#' Read metadata from file
#'
#' Reads metadata from a file, guessing the separator automatically.
#'
#' @param filepath Path to the file to read
#' @param sepguesser_function Function to use for guessing the separator (defaults to sepguesser)
#'
#' @return A data frame with metadata
#' @export
#'
read_metadata <- function(filepath, sepguesser_function = sepguesser) {
  guessed_sep <- sepguesser_function(filepath)
  expdesign <- utils::read.delim(filepath, header = TRUE,
                                 sep = guessed_sep, quote = "",
                                 check.names = FALSE, stringsAsFactors = TRUE)
  if (colnames(expdesign)[1] == "") {
    rownames(expdesign) <- expdesign[,1]
    expdesign <- expdesign[,-1]
  }
  return(expdesign)
}

#' Create DESeqDataSet from count matrix and metadata
#'
#' Creates a DESeqDataSet object from count matrix and experimental design data.
#'
#' @param countmatrix Matrix with count data
#' @param expdesign Data frame with experimental design
#' @param design_formula Formula for the experimental design
#' @param gene_filter Regular expression to filter genes (default: "")
#' @param dds_intercept Intercept value for the design (default: NULL)
#' @param design_factor Name of the design factor column (default: NULL)
#'
#' @return A DESeqDataSet object
#' @export
#'
create_dds <- function(countmatrix, expdesign, design_formula, gene_filter = "", 
                       dds_intercept = NULL, design_factor = NULL) {
  # Get common samples
  comSamples <- intersect(colnames(countmatrix), rownames(expdesign))
  
  # Apply gene filter if provided
  if (nchar(gene_filter) > 0) {
    countmatrix <- countmatrix[grep(gene_filter, rownames(countmatrix), invert = TRUE), ]
  }
  
  # Subset data to common samples
  expdesign <- expdesign[comSamples, ]
  countmatrix <- countmatrix[, comSamples]
  
  # Relevel factor if intercept is provided
  if (!is.null(dds_intercept) && !is.null(design_factor)) {
    expdesign[, design_factor] <- relevel(expdesign[, design_factor], ref = dds_intercept)
  }
  
  # Add pseudocount if all rows contain zeros
  if (all(rowSums(countmatrix == 0) > 0)) {
    countmatrix <- countmatrix + 1
  }
  
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = countmatrix,
    colData = expdesign,
    design = design_formula
  )
  
  # Estimate size factors
  dds <- estimateSizeFactors(dds)
  
  return(dds)
}

#' Generate co-occurrence plots
#'
#' Creates co-occurrence plots for experimental design visualization.
#'
#' @param expdesign Data frame with experimental design
#' @param design_formula Formula for the experimental design
#'
#' @return A ggplot object with co-occurrence plots
#' @export
#'
generate_cooccurrence_plots <- function(expdesign, design_formula) {
  # Validate inputs
  if (is.null(expdesign) || is.null(design_formula)) {
    return(NULL)
  }
  
  # Create model matrix
  expanded_formula <- terms(design_formula, data = expdesign)
  model_matrix <- model.matrix(expanded_formula, data = expdesign)
  
  # Prepare data for plotting
  longData <- reshape2::melt(model_matrix)
  longData <- longData[longData$value != 0, ]
  
  # Create main heatmap plot
  p1 <- ggplot(longData, aes(x = Var1, y = Var2)) + 
    geom_raster(aes(fill = value)) + 
    scale_fill_gradient(low = "grey90", high = "darkgrey") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.y = element_text(size = 14),
          legend.position = "none") 
  
  # Create sample counts plot
  sample_counts <- as.data.frame(colSums(model_matrix != 0))
  colnames(sample_counts) <- "count"
  sample_counts$names <- rownames(sample_counts)
  sample_counts$names <- factor(sample_counts$names, levels = levels(longData$Var2))
  sample_counts <- sample_counts[levels(longData$Var2), ]
  
  p2 <- ggplot(sample_counts, aes(x = names, y = count)) + 
    geom_bar(stat = "identity") + coord_flip() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  # Combine plots
  cowplot::plot_grid(p1, p2, rel_widths = c(4, 1))
}

#' Generate DE gene lists
#'
#' Extract upregulated, downregulated, and all DE genes from results.
#'
#' @param res_obj DESeqResults object
#' @param FDR False discovery rate threshold
#' @param annotation_obj Annotation data frame (optional)
#'
#' @return A list with UP, DOWN, and UPDOWN gene lists
#' @export
#'
generate_gene_lists <- function(res_obj, FDR = 0.05, annotation_obj = NULL) {
  # Convert results to DE genes
  res_tbl <- deseqresult2DEgenes(res_obj, FDR = FDR)
  
  if (nrow(res_tbl) < 1) {
    return(list(UP = NULL, DOWN = NULL, UPDOWN = NULL))
  }
  
  # Separate up and down regulated genes
  res_tbl_UP <- res_tbl[res_tbl$log2FoldChange > 0 & !is.na(res_tbl$padj), ]
  res_tbl_DOWN <- res_tbl[res_tbl$log2FoldChange < 0 & !is.na(res_tbl$padj), ]
  
  # Handle annotation if provided
  if (!is.null(annotation_obj) && "symbol" %in% colnames(res_obj)) {
    # Add symbols for UP genes
    if (nrow(res_tbl_UP) > 0) {
      res_tbl_UP$symbol <- annotation_obj$gene_name[
        match(res_tbl_UP$id, rownames(annotation_obj))]
      listUP <- res_tbl_UP$symbol
    } else {
      listUP <- NULL
    }
    
    # Add symbols for DOWN genes
    if (nrow(res_tbl_DOWN) > 0) {
      res_tbl_DOWN$symbol <- annotation_obj$gene_name[
        match(res_tbl_DOWN$id, rownames(annotation_obj))]
      listDOWN <- res_tbl_DOWN$symbol
    } else {
      listDOWN <- NULL
    }
    
    # Add symbols for all DE genes
    if (nrow(res_tbl) > 0) {
      res_tbl$symbol <- annotation_obj$gene_name[
        match(res_tbl$id, rownames(annotation_obj))]
      listUPDOWN <- res_tbl$symbol
    } else {
      listUPDOWN <- NULL
    }
  } else {
    listUP <- res_tbl_UP$symbol
    listDOWN <- res_tbl_DOWN$symbol
    listUPDOWN <- res_tbl$symbol
  }
  
  return(list(UP = listUP, DOWN = listDOWN, UPDOWN = listUPDOWN))
}