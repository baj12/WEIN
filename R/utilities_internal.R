#' Internal Utility Functions
#'
#' These functions are used internally by the WEIN package and are not meant
#' to be called directly by users.
#'
#' @docType package
#' @name utilities_internal
#' @keywords internal
"_PACKAGE"

#' Row scaling function
#'
#' Scale matrix rows by subtracting the mean and dividing by the standard deviation
#'
#' @param x A matrix or data frame to be scaled
#'
#' @return A scaled matrix with the same dimensions as the input
#' @keywords internal
#'
mat_rowscale <- function(x) {
  m <- apply(x, 1, mean, na.rm = TRUE)
  s <- apply(x, 1, sd, na.rm = TRUE)
  return((x - m) / s)
}

#' Get gene list
#'
#' Retrieve a gene list from the provided values based on the list name
#'
#' @param list_name Name of the gene list to retrieve
#' @param values Reactive values containing the gene lists
#'
#' @return The requested gene list or NULL if not found
#' @keywords internal
#'
get_gene_list <- function(list_name, values) {
  switch(list_name,
         "UP" = values$genelistUP(),
         "DOWN" = values$genelistDOWN(),
         "UPDOWN" = values$genelistUPDOWN(),
         "LIST1" = as.character(values$genelist1$`Gene.Symbol`),
         "LIST2" = as.character(values$genelist2$`Gene.Symbol`),
         "LIST3" = as.character(values$genelist3$`Gene.Symbol`),
         "LIST4" = as.character(values$genelist4$`Gene.Symbol`)
  )
}

#' Create GO term heatmap
#'
#' Create a heatmap for selected GO terms and genes
#'
#' @param topgo_data TopGO data containing GO terms and genes
#' @param selected_row Selected row index from the data
#' @param values Reactive values containing dds_obj, cur_species, and cur_type
#' @param annoSpecies_df Data frame with species information
#'
#' @return A heatmaply plot or NULL if no row is selected
#' @keywords internal
#'
create_goterm_heatmap <- function(topgo_data, selected_row, values, annoSpecies_df) {
  if(length(selected_row) == 0) return(NULL)
  
  # Extract genes and term (take first row if multiple selected)
  mygenes <- topgo_data[selected_row, ]$genes[1]
  myterm <- paste0(
    topgo_data[selected_row, ]$`GO.ID`, " - ",
    topgo_data[selected_row, ]$Term
  )
  
  # Convert gene symbols to IDs
  genevec <- unlist(strsplit(mygenes, split = ","))
  annopkg <- annoSpecies_df$pkg[annoSpecies_df$species == values$cur_species]
  genevec_ids <- as.character(mapIds(
    eval(parse(text = annopkg)),
    genevec,
    values$cur_type,
    "SYMBOL",
    multiVals = "first"
  ))
  
  # Get log2 transformed values
  log2things <- assay(normTransform(values$dds_obj))
  selectedLogvalues <- log2things[genevec_ids, ]
  
  # Determine row labels
  rowlabs <- if(length(genevec_ids) == length(genevec)) {
    genevec
  } else {
    genevec_ids
  }
  
  heatmaply(selectedLogvalues, scale = "row", labels_row = rowlabs, main = myterm)
}

#' Create gene list handler
#'
#' Create a reactive and observer for handling gene list uploads
#'
#' @param list_num Number of the gene list (1-4)
#' @param input Shiny input object
#' @param values Reactive values
#' @param ns Namespace function
#' @param read1stCol Function to read the first column of a file
#'
#' @return A reactive function for the gene list
#' @keywords internal
#'
create_gene_list_handler <- function(list_num, input, values, ns, read1stCol) {
  gl_name <- paste0("gl", list_num)
  genelist_name <- paste0("genelist", list_num)
  
  # Create reactive
  gl_reactive <- reactive({
    shiny::validate(
      need(!is.null(input[[gl_name]]), "No gene list file has been uploaded.")
    )
    gl <- read1stCol(input[[gl_name]]$datapath, values$dds_obj)
    shiny::validate(
      need(!is.null(gl),
           "Failed to read gene list file. Please ensure the file is properly formatted with gene identifiers in the first column.")
    )
    return(gl)
  })
  
  # Create observer
  observeEvent(input[[gl_name]], {
    gl <- gl_reactive()
    shiny::validate(
      need(!is.null(gl) && nrow(gl) >= 1, "Gene list is empty or invalid.")
    )
    mydf <- as.data.frame(gl, stringsAsFactors = FALSE)
    names(mydf) <- "Gene Symbol"
    values[[genelist_name]] <- mydf
  })
  
  return(gl_reactive)
}

#' Create UI output
#'
#' Create a UI output for displaying data tables
#'
#' @param ns Namespace function
#' @param values Reactive values
#' @param value_name Name of the value to display
#' @param title Title for the output
#' @param output_name Name of the output
#'
#' @return A UI output
#' @keywords internal
#'
create_ui_output <- function(ns, values, value_name, title, output_name) {
  renderUI({
    shiny::validate(
      need(!is.null(values[[value_name]]),
           paste("Value", value_name, "is not available"))
    )
    tagList(
      h4(title),
      DT::dataTableOutput(ns(output_name))
    )
  })
}

#' Create DataTable output
#'
#' Create a DataTable output for displaying data
#'
#' @param ns Namespace function
#' @param values Reactive values
#' @param value_name Name of the value to display
#' @param link_column Column to create links for
#'
#' @return A DataTable output
#' @keywords internal
#'
create_dt_output <- function(ns, values, value_name, link_column = "rownames") {
  DT::renderDataTable({
    shiny::validate(
      need(!is.null(values[[value_name]]),
           paste("Value", value_name, "is not available"))
    )
    mytbl <- values[[value_name]]
    
    if (link_column == "rownames") {
      rownames(mytbl) <- createLinkGO(rownames(mytbl))
    } else {
      mytbl[[link_column]] <- createLinkGO(mytbl[[link_column]])
    }
    
    mytbl
  }, escape = FALSE, options = list(scrollX = TRUE))
}