ui_setup_server <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    # server info boxes -----------------------------------------------------------
    output$box_ddsobj <- renderUI({
      if(!is.null(values$dds_obj))
        return(valueBox("dds object",
                        paste0(nrow(values$dds_obj), " genes - ",ncol(values$dds_obj)," samples"),
                        icon = icon("list"),
                        color = "green",width = NULL))
      else
        return(valueBox("dds object",
                        "yet to create",
                        icon = icon("list"),
                        color = "red",width = NULL))
      
    })
    
    output$box_annobj <- renderUI({
      if(!is.null(values$annotation_obj))
        return(valueBox("Annotation",
                        paste0(nrow(values$annotation_obj), " genes - ",ncol(values$annotation_obj)," ID types"),
                        icon = icon("book"),
                        color = "green",width = NULL))
      else
        return(valueBox("Annotation",
                        "yet to create",
                        icon = icon("book"),
                        color = "red",width = NULL))
    })
    
    output$box_resobj <- renderUI({
      if(!is.null(values$res_obj)){
        DEregu <- sum(values$res_obj$padj < input$FDR & values$res_obj$log2FoldChange != 0, na.rm = TRUE)
        return(valueBox("DE genes",
                        paste0(DEregu, " DE genes - out of ",nrow(values$res_obj),""),
                        icon = icon("rectangle-list"),
                        color = "green",width = NULL))
      } else
        return(valueBox("DE genes",
                        "yet to create",
                        icon = icon("rectangle-list"),
                        color = "red",width = NULL))
    })
    
    
    # server ok objects -----------------------------------------------------------
    output$ok_cm <- renderUI({
      if (is.null(values$countmatrix))
        return(NULL)
      # icon("check",class = "icon-done") # this does not allow to set the size? go manually with..
      tags$div(HTML('<i class="fa fa-check fa-3x icon-done"></i>'))
    })
    output$ok_ed <- renderUI({
      if (is.null(values$dds_obj))
        return(NULL)
      # icon("check",class = "icon-done") # this does not allow to set the size? go manually with..
      tags$div(HTML('<i class="fa fa-check fa-3x icon-done"></i>'))
    })
    output$ok_dds <- renderUI({
      if (is.null(values$dds_obj))
        return(NULL)
      # icon("check",class = "icon-done") # this does not allow to set the size? go manually with..
      tags$div(HTML('<i class="fa fa-check fa-3x icon-done"></i>'))
    })
    output$ok_anno <- renderUI({
      if (is.null(values$annotation_obj))
        return(NULL)
      # icon("check",class = "icon-done") # this does not allow to set the size? go manually with..
      tags$div(HTML('<i class="fa fa-check fa-3x icon-done"></i>'))
    })
    output$ok_ddsRun <- renderUI({
      if(is.null(values$dds_obj))
        return(NULL)
      if (!"results" %in% mcols(mcols(values$dds_obj))$type) #
        return(NULL)
      tags$div(HTML('<i class="fa fa-check fa-3x icon-done"></i>'))
    })
    output$ok_resu <- renderUI({
      if (is.null(values$res_obj))
        return(NULL)
      # icon("check",class = "icon-done") # this does not allow to set the size? go manually with..
      tags$div(HTML('<i class="fa fa-check fa-3x icon-done"></i>'))
    })
    output$nonZeroCountsPlot <- renderPlot({
      # browser()
      if(is.null(values$dds_obj)) return(NULL)
      counts <- assays(values$dds_obj)[["counts"]]
      # colSums sums over each column producing a vector of counts.
      countsNz <- colSums(counts)
      # adjust the maximum for being able to plot numbers on top of the bars
      ylim <- c(0, 1.3 * max(countsNz))
      # las=2 rotates the labels
      par(mar = c(4, 10, 4, 2) + 0.1)
      xx <- barplot(countsNz,
                    xlim = ylim, main = "number of non-zero countsNz",
                    las = 2, horiz = T, cex.names = 0.75
      )
      text(
        y = xx, x = countsNz, label = prettyNum(countsNz, big.mark = ","),
        pos = 4, cex = 0.6, col = "darkgreen"
      )
    })
    output$alignedSequencesPlot <- renderPlot({
      if(is.null(values$dds_obj)) return(NULL)
      counts <- assays(values$dds_obj)[["counts"]]
      ylim <- c(0, 1.4 * max(colSums(counts)))
      op <- par(mar = c(4, 10, 4, 2) + 0.1)
      
      # Here we use a different orientation of the bar plot.
      xx <- barplot(colSums(counts),
                    xlim = ylim, main = "Number of Aligned Sequences",
                    horiz = T, las = 2, cex.names = 0.75
      )
      text(
        y = xx, x = colSums(counts), label = prettyNum(colSums(counts),
                                                       big.mark = ","
        ), pos = 4,
        cex = 0.8, col = "red"
      )
    })
    
  })
}
