# helpers.R





footer <- function(){
  tags$div(
    class = "panel-footer",
    style = "text-align:center",
    tags$div(
      class = "foot-inner",
      list(
        # hr(),
        "WEIN (Web-based Engine for Interactive Next-generation sequencing analysis) is based on the ideal project by Federico Marini ",
        tags$a(href="https://github.com/federicomarini/ideal", "GitHub"),
        ". This is further developed by Bernd Jagla at the Institut Pasteur.",
        br(),
        "",
        "License: ",tags$a(href="https://opensource.org/licenses/MIT","MIT"), br(),

        "Development of the WEIN package is on ",
        tags$a(href="https://github.com/baj12/idealImmunoTP", "GitHub"), "."
        
      )
    )
  )
}







############################# helper funcs #################################


read1stCol <- function (fileName,dds_obj){
  guessed_sep <- sepguesser(fileName)
  
  cm <- tryCatch({utils::read.delim(fileName, header = TRUE,
                                    as.is = TRUE, sep = guessed_sep, 
                                    # row.names = 1, # https://github.com/federicomarini/pcaExplorer/issues/1
                                    ## TODO: tell the user to use tsv, or use heuristics
                                    ## to check what is most frequently occurring separation character? -> see sepGuesser.R
                                    check.names = FALSE)
  }, error=function(e){
    cat(file = stderr(), paste(e,"\n"))
    return(NULL)
  }
  )
  if(is.null(cm)) return(NULL)
  if (ncol(cm) >1) {
    cm = cm[,1, drop = F]
  }
  cm = cm[cm[,1] %in% rownames(dds_obj),, drop = F]
  if(nrow(cm)<1) return(NULL)
  return(cm)
}



#' Make an educated guess on the separator character
#'
#' This function tries to guess which separator was used in a text delimited file
#'
#' @param file The name of the file which the data are to be read from
#' @param sep_list A vector containing the candidates for being identified as
#' separators. Defaults to \code{c(",", "\t", ";"," ")}
#'
#' @return A character value, corresponding to the guessed separator. One of ","
#' (comma), "\\t" (tab), ";" (semicolon)," " (whitespace)
#' @export
#'
#' @examples
#' sepguesser(system.file("extdata/design_commas.txt",package = "ideal"))
#' sepguesser(system.file("extdata/design_semicolons.txt",package = "ideal"))
#' sepguesser(system.file("extdata/design_spaces.txt",package = "ideal"))
#' mysep <- sepguesser(system.file("extdata/design_tabs.txt",package = "ideal"))
#'
#' # to be used for reading in the same file, without having to specify the sep
#'
sepguesser <- function(file, sep_list = c(",", "\t", ";"," ")) {
  separators_list = sep_list
  rl = readLines(file, warn = FALSE)
  rl = rl[rl != ""] # allow last line to be empty
  sephits_min = sapply(separators_list, function(x) min(stringr::str_count(rl, x))) #minimal number of separators on all lines
  sep = separators_list[which.max(sephits_min)]
  sep
}

sepguesser2 <- function(file, sep_list = c(",", "\t", ";"," ")) {
  separators_list = sep_list
  rl = readLines(file, warn = FALSE)
  rl = rl[rl != ""] # allow last line to be empty
  sephits_min = sapply(separators_list, function(x) min(stringr::str_count(rl, x))) #minimal number of separators on all lines
  
  counts <- sapply(separators_list, function(x) min(count.fields(textConnection(rl), sep=x)))
  
  sep = separators_list[which.max(counts)]
  # sep = separators_list[which.max(sephits_min)]
  sep
}





# combineTogether <- function(normCounts,resuTable,anns) {
#   combinedCountsAndRes <- inner_join(resuTable,normCounts,by="id")
#   anns2 <- anns[match(combinedCountsAndRes$id, anns[, 1]), ]
#   combinedCountsAndRes$Description <- anns2$description
#   return(combinedCountsAndRes)
# }
#
#
# combine_resucounts <- function(normCounts,resuTable) {
#   combinedCountsAndRes <- inner_join(resuTable,normCounts,by="id")
#   # anns2 <- anns[match(combinedCountsAndRes$id, anns[, 1]), ]
#   # combinedCountsAndRes$Description <- anns2$description
#   return(combinedCountsAndRes)
# }
#

# getGeneInfos <- function(obj, annopkg, idtype) {
#   # obj is a dds object...
#   ids <- rownames(obj)
#
#   mydf <- mapIds(eval(parse(text=annopkg)),keys=ids,column = "GENENAME",keytype = idtype)
#   mydf_2 <- AnnotationDbi::select(eval(parse(text=annopkg)),keys=ids,column = "GENENAME",keytype = idtype)
#
#   return(mydf_2)
# }
#
#

