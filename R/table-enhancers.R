# table-enhancers.R



createLinkGO <- function(val) {
  sprintf('<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank" class="btn btn-primary">%s</a>',val,val)
}

createLinkENS  <- function(val, species="Mus_musculus") {
  paste0('<a href="http://www.ensembl.org/',species,'/Gene/Summary?g=',val,'" target="_blank" class="btn btn-primary">',val,'</a>')
}

createLinkGeneSymbol <- function(val) {
  # possibilities:
  # ncbi
  # genecards
  paste0('<a href="http://www.ncbi.nlm.nih.gov/gene/?term=',val,'[sym]" target="_blank" class="btn btn-primary">',val,'</a>')
}



geneinfo <- function(gene_id) {
  # the gene id has to be entrez_id
  
  # Implementation using rentrez to fetch gene information
  shiny::validate(
    need(!is.null(gene_id) && gene_id != "", "Empty gene ID provided")
  )
  
  tryCatch({
    entrezinfo <- rentrez::entrez_summary("gene", gene_id)
    return(entrezinfo)
  }, error = function(e) {
    warning("Failed to retrieve gene information for ID: ", gene_id, ". Error: ", e$message)
    shiny::validate(
      need(FALSE, paste("Failed to retrieve gene information for ID:", gene_id, ". Error:", e$message))
    )
  })
}




