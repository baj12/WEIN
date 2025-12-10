summary_plots_server <- function(id, values, annoSpecies_df, exportPlots) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    # not used?
    output$explore_res <- renderPrint({
      expfac <- attributes(terms.formula(design(values$dds_obj)))$term.labels
      expfac # plus, support up to four factors that are either there or not according to the length
    })
    
    output$plotma <- renderPlot({
      p <- plot_ma(values$res_obj,annotation_obj = values$annotation_obj,FDR = values$FDR)
      exportPlots$plot_ma <- p
      p
    })
    
    # Throttling for brush events to prevent excessive re-rendering
    brush_throttle <- reactiveVal(0)
    
    observeEvent(input$ma_brush, {
      # Throttle by 200ms
      invalidateLater(200, session)
      brush_throttle(isolate(brush_throttle()) + 1)
    })
    
    output$mazoom <- renderPlot({
      # Use throttle to control frequency
      brush_throttle()
      
      if(is.null(input$ma_brush)) return(ggplot() + annotate("text",label="click and drag to zoom in",0,0) + theme_bw())
      
      if(!is.null(values$annotation_obj))
        p <- plot_ma(values$res_obj,annotation_obj = values$annotation_obj,FDR = values$FDR) +
          coord_cartesian(xlim = c(input$ma_brush$xmin,input$ma_brush$xmax),
                          ylim = c(input$ma_brush$ymin,input$ma_brush$ymax)) +
          geom_text(aes_string(label="genename"),size=input$size_genelabels,hjust=0.25, vjust=-0.75)
      else
        p <-  plot_ma(values$res_obj,annotation_obj = values$annotation_obj,FDR = values$FDR) +
          coord_cartesian(xlim = c(input$ma_brush$xmin,input$ma_brush$xmax),
                          ylim = c(input$ma_brush$ymin,input$ma_brush$ymax))
      exportPlots$plot_mazoom <- p
      p
    })
    
     
    
    curData <- reactive({
      mama <- data.frame(mean=values$res_obj$baseMean,lfc=values$res_obj$log2FoldChange,padj = values$res_obj$padj,isDE= ifelse(is.na(values$res_obj$padj), FALSE, values$res_obj$padj < 0.10),ID=rownames(values$res_obj))
      mama$genename <- values$annotation_obj$gene_name[match(mama$ID,rownames(values$annotation_obj))]
      # mama$yesorno <- ifelse(mama$isDE,"yes","no")
      mama$yesorno <- ifelse(mama$isDE,"red","black")
      mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
      res <- brushedPoints(mama, input$ma_brush,xvar="logmean",yvar="lfc")
      res
    })
    
    
    # Throttling for click events to prevent excessive re-rendering
    click_throttle <- reactiveVal(0)
    
    observeEvent(input$mazoom_click, {
      # Throttle by 300ms
      invalidateLater(300, session)
      click_throttle(isolate(click_throttle()) + 1)
    })
    
    curDataClick <- reactive({
      # Use throttle to control frequency
      click_throttle()
      
      mama <- data.frame(mean=values$res_obj$baseMean,lfc=values$res_obj$log2FoldChange,padj = values$res_obj$padj,isDE= ifelse(is.na(values$res_obj$padj), FALSE, values$res_obj$padj < 0.10),ID=rownames(values$res_obj))
      mama$genename <- values$annotation_obj$gene_name[match(mama$ID,rownames(values$annotation_obj))]
      # mama$yesorno <- ifelse(mama$isDE,"yes","no")
      mama$yesorno <- ifelse(mama$isDE,"red","black")
      mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
      res <- nearPoints(mama, input$mazoom_click,threshold = 20, maxpoints = 1,
                        addDist = TRUE)
      res
    })
    
    
    
    
    output$ma_brush_out <- DT::renderDataTable({
      if(nrow(curData())==0)
        return(NULL)
      # browser()
      # curData()
      brushedObject <- curData()
      selectedGenes <- as.character(brushedObject$ID)
      rownames(brushedObject) = selectedGenes
      brushedObject
    },options=list(pageLength=100))
    
    output$selectedGenesMAplot <- renderText({
      brushedObject <- curData()
      selectedGenes <- as.character(brushedObject$ID)
      paste(selectedGenes, collapse = " ")
    })
    output$heatbrush <- renderPlotly({
      if((is.null(input$ma_brush))|is.null(values$dds_obj)) {
        return(NULL)
      }
      #
      
      brushedObject <- curData()
      selectedGenes <- as.character(brushedObject$ID)
      toplot <- assay(values$dds_obj)[selectedGenes,]
      rownames(toplot) <- values$annotation_obj$gene_name[match(rownames(toplot),rownames(values$annotation_obj))]
      
      if(input$pseudocounts) toplot <- log2(1+toplot)
      if(input$rowscale) toplot <- mat_rowscale(toplot)
      
      # Add early stopping for very large gene sets to prevent performance issues
      if(nrow(toplot) > 1000) {
        showNotification("Large gene set detected. Displaying first 1000 genes for performance.",
                             type = "warning", duration = 10)
        toplot <- toplot[1:min(1000, nrow(toplot)), , drop = FALSE]
      }
      
      heatmaply::heatmaply(toplot,cluster_cols = as.logical(input$heatmap_colv))
    })
    
    
    output$hpi_brush <- renderPlotly({
      # if((is.null(input$ma_brush))|is.null(values$dds_obj)) {
      #   # plot(100:1)
      # }
      #return(NULL)
      # browser()
      brushedObject <- curData()
      selectedGenes <- as.character(brushedObject$ID)
      toplot <- assay(values$dds_obj)[selectedGenes,]
      rownames(toplot) <- values$annotation_obj$gene_name[match(rownames(toplot),rownames(values$annotation_obj))]
      mycolss <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027","#a50026") # to be consistent with red/blue usual coding
      if(input$pseudocounts) toplot <- log2(1+toplot)
      if(input$rowscale) toplot <- mat_rowscale(toplot)
      if(nrow(toplot) <1){
        cat(file = stderr(), "\ntoplot nrow <1\n")
        return(NULL)
      }
      # Add early stopping for very large gene sets to prevent performance issues
      if(nrow(toplot) > 1000) {
        showNotification("Large gene set detected. Displaying first 1000 genes for performance.",
                             type = "warning", duration = 10)
        toplot <- toplot[1:min(1000, nrow(toplot)), , drop = FALSE]
      }
      
      heatmaply::heatmaply(toplot,Colv = as.logical(input$heatmap_colv),colors = mycolss, cexCol = 1)
    })
    
    
    output$deb <- renderPrint({
      # curDataClick()
      selectedGene <- curDataClick()$ID
      #         selectedGeneSymbol <- cm2$fromgtf[match(selectedGene,rownames(cm2))]
      #         # plotCounts(dds_cleaner,)
      #         genedata <- plotCounts(dds_cleaner,gene=selectedGene,intgroup = "condition",returnData = T)
      #         genedata
      # str(as.character(selectedGene))
      selectedGene
    })
    
    output$volcanoplot <- renderPlotly({
      p <- plot_volcano(values$res_obj, FDR = values$FDR)
      exportPlots$plot_volcanoplot <- p
      plotly::ggplotly(p)
    })
    
    # server genefinder --------------------------------------------------------
    output$genefinder_plot <- renderPlot({
      
      shiny::validate(
        need(
          length(values$color_by)>0,
          "Select an experimental factor in the Group/color by element in the sidebar"
        )
      )
      shiny::validate(
        need(
          !is.null(values$annotation_obj),
          "Please set the annotation"
        )
      )
      
      shiny::validate(
        need(!is.null(input$ma_brush),
             "Please select a region on the MA plot")
      )
      shiny::validate(
        need(!is.null(input$mazoom_click),
             "Please click on a gene in the zoomed MA plot")
      )
      
      
      selectedGene <- as.character(curDataClick()$ID)
      selectedGeneSymbol <- values$annotation_obj$gene_name[match(selectedGene,values$annotation_obj$gene_id)]
      
      p <- ggplotCounts(values$dds_obj, selectedGene, intgroup = values$color_by, annotation_obj=values$annotation_obj)
      
      if(input$ylimZero_genes)
        p <- p + ylim(0.1, NA)
      
      exportPlots$plot_genefinder <- p
      p
    })
    
    output$rentrez_infobox <- renderUI({
      shiny::validate(
        need(
          (nrow(curDataClick()) > 0),
          "Select a gene first to display additional info (retrieved from the NCBI/ENTREZ db website)"
        )
      )
      shiny::validate(
        need(
          (!is.null(values$cur_species)),
          "Select a species first in the Data Setup panel"
        )
      )
      # browser()
      selectedGene <- as.character(curDataClick()$ID)
      selgene_entrez <- mapIds(get(annoSpecies_df[values$cur_species,]$pkg),
                               selectedGene, "ENTREZID", values$cur_type)
      fullinfo <- geneinfo(selgene_entrez)
      
      # Build up link manually to paste under the info
      link_pubmed <- paste0('<a href="http://www.ncbi.nlm.nih.gov/gene/?term=',
                            selgene_entrez,
                            '" target="_blank" >Click here to see more at NCBI</a>')
      
      if(fullinfo$summary == "")
        return(HTML(paste0("<b>",fullinfo$name, "</b><br/><br/>",
                           fullinfo$description,"<br/><br/>",
                           link_pubmed
        )))
      else
        return(HTML(paste0("<b>",fullinfo$name, "</b><br/><br/>",
                           fullinfo$description, "<br/><br/>",
                           fullinfo$summary, "<br/><br/>",
                           link_pubmed
        )))
    })
    
    
    cur_combires <- reactive({
      shiny::validate(
        need(!is.null(values$res_obj),
             "Results object is not available")
      )
      
      
      normCounts <- as.data.frame(counts(estimateSizeFactors(values$dds_obj),normalized=TRUE))
      normCounts$id <- rownames(normCounts)
      res_df <- deseqresult2tbl(values$res_obj)
      
      combi_obj <- dplyr::inner_join(res_df,normCounts,by="id")
      combi_obj$symbol <- values$annotation_obj$gene_name[match(combi_obj$id,values$annotation_obj$gene_id)]
      
      
      if("symbol" %in% names(values$res_obj)) {
        sel_genes <- values$avail_symbols
        sel_genes_ids <- values$annotation_obj$gene_id[match(sel_genes,values$annotation_obj$gene_name)]
      } else {
        sel_genes_ids <- values$avail_ids
      }
      
      if(length(sel_genes_ids) > 0) {
        combi_obj[match(sel_genes_ids,combi_obj$id),]
      } else {
        combi_obj
      }
    })
    
    output$table_combi <- DT::renderDataTable({
      datatable(cur_combires(),options = list(scrollX=TRUE))
    })
    
     
    output$download_plot_ma <- downloadHandler(filename = function() {
      input$filename_plot_ma
    }, content = function(file) {
      ggsave(file, exportPlots$plot_ma, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    output$download_plot_mazoom <- downloadHandler(filename = function() {
      input$filename_plot_mazoom
    }, content = function(file) {
      ggsave(file, exportPlots$plot_mazoom, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    
     
    output$download_plot_volcanoplot <- downloadHandler(filename = function() {
      input$filename_plot_volcanoplot
    }, content = function(file) {
      ggsave(file, exportPlots$plot_volcanoplot, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    output$download_plot_genefinder <- downloadHandler(filename = function() {
      input$filename_plot_genefinder
    }, content = function(file) {
      ggsave(file, exportPlots$plot_genefinder, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    # base graphics plots
    output$download_plot_heatbrush <- downloadHandler(filename = function() {
      input$filename_plot_heatbrush
    }, content = function(file) {
      pdf(file)
      brushedObject <- curData()
      
      selectedGenes <- as.character(brushedObject$ID)
      toplot <- assay(values$dds_obj)[selectedGenes,]
      rownames(toplot) <- values$annotation_obj$gene_name[match(rownames(toplot),rownames(values$annotation_obj))]
      
      if(input$pseudocounts) toplot <- log2(1+toplot)
      
      if(input$rowscale) toplot <- mat_rowscale(toplot)
      
      heatmaply(toplot,cluster_cols = as.logical(input$heatmap_colv))
      dev.off()
    })
    # tbls
    output$downloadTblMabrush <- downloadHandler(
      filename = function() {
        "table_mabrush.csv"
      },
      content = function(file) {
        write.csv(curData(), file)
      }
    )
    
  })
}