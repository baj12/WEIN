count_overview_server <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    # server counts overview --------------------------------------------------------
    current_countmat <- reactive({
      cat(file = stderr(), paste("current_countmat\n"))
      
      if(input$countstable_unit=="raw_counts")
        return(counts(values$dds_obj,normalized=FALSE))
      if(input$countstable_unit=="normalized_counts")
        return(counts(values$dds_obj,normalized=TRUE))
      if(input$countstable_unit=="rlog_counts")
        return(NULL) ## see if it is worth to keep in here or explore possibility with fast vst
      if(input$countstable_unit=="log10_counts")
        return(log10(1 + counts(values$dds_obj,normalized=TRUE)))
      if(input$countstable_unit=="tpm_counts")
        return(NULL) ## TODO!: assumes length of genes/exons as known, and is currently not required in the dds
      
    })
    
    output$showcountmat <- DT::renderDataTable({
      datatable(current_countmat())
    })
    output$downloadData <- downloadHandler(
      filename = function() {
        paste0(input$countstable_unit,"table.csv")
      },
      content = function(file) {
        write.csv(current_countmat(), file)
      }
    )
    
    output$corrplot <- renderPlot({
      cat(file = stderr(), paste("corrplot\n"))
      
      if(input$compute_pairwisecorr)
        withProgress(
          pair_corr(current_countmat(),
                    method=input$corr_method,
                    log = input$corr_uselogs,
                    use_subset = input$corr_usesubset),
          message = "Preparing the plot",
          detail = "this can take a while..."
        )
    })
    # HEATCORR ======
    output$heatcorr <- renderPlotly({
      if(input$compute_pairwisecorr){
        input$compute_pairwisecorr
        values$avail_symbols
        values$dds_obj
        values$color_by
        heatmaply::heatmaply_cor(cor(current_countmat()))
      }
    })
    
    output$pcaPlot <- renderPlotly({
      shiny::validate(
        need(!is.null(values$dds_obj),
             "Provide or construct a dds object")
      )
      shiny::validate(
        need(!is.null(values$color_by),
             "Provide group to color by")
      )
      rld <- vst(values$dds_obj, blind = FALSE,nsub=10)
      p = multiAxPCA(rld, intgroup = values$color_by, ntop = 1000, pc1=as.numeric(input$pcaDim1), pc2=as.numeric(input$pcaDim2))
      p2 = p +aes(text=colnames(rld))
      ggplotly(p2, tooltip = "text")
    })
    
    output$pcaEVPlot <- renderPlotly({
      shiny::validate(
        need(!is.null(values$dds_obj),
             "Provide or construct a dds object")
      )
      shiny::validate(
        need(!is.null(values$color_by),
             "Provide group to color by")
      )
      # browser()
      # rld <- rlog(values$dds_obj, blind = FALSE)
      rld <- vst(values$dds_obj, blind = FALSE,nsub=10)
      # browser()
      # multiAxPCA = function (object, intgroup = "condition", ntop = 500, returnData = FALSE, pc1=1, pc2=2) 
      pD = multiAxPCA(rld, intgroup = values$color_by, ntop = 1000, pc1=1, pc2=2, returnData = T)
      df = data.frame(PC=seq(pD),EV=pD)
      # browser()
      p<-ggplot(data=df, aes(x=PC, y=EV)) +
        geom_bar(stat="identity")
      p
    })
    output$pcaPlot34 <- renderPlotly({
      shiny::validate(
        need(!is.null(values$dds_obj),
             "Provide or construct a dds object")
      )
      shiny::validate(
        need(!is.null(values$color_by),
             "Provide group to color by")
      )
      # browser()
      # rld <- rlog(values$dds_obj, blind = FALSE)
      rld <- vst(values$dds_obj, blind = FALSE,nsub=10)
      # browser()
      # multiAxPCA = function (object, intgroup = "condition", ntop = 500, returnData = FALSE, pc1=1, pc2=2) 
      p = multiAxPCA(rld, intgroup = values$color_by, ntop = 1000, pc1=3, pc2=4)
      p2 = p +aes(text=colnames(rld))
      ggplotly(p2, tooltip = "text")
    })
    output$pcaPlot56 <- renderPlotly({
      shiny::validate(
        need(!is.null(values$dds_obj),
             "Provide or construct a dds object")
      )
      shiny::validate(
        need(!is.null(values$color_by),
             "Provide group to color by")
      )
      # browser()
      # rld <- rlog(values$dds_obj, blind = FALSE)
      rld <- vst(values$dds_obj, blind = FALSE,nsub=10)
      # browser()
      # multiAxPCA = function (object, intgroup = "condition", ntop = 500, returnData = FALSE, pc1=1, pc2=2) 
      p = multiAxPCA(rld, intgroup = values$color_by, ntop = 1000, pc1=5, pc2=6)
      p2 = p +aes(text=colnames(rld))
      ggplotly(p2, tooltip = "text")
    })
    
    output$sizeFactorsPlot <- renderPlot({
      
      col <- RColorBrewer::brewer.pal(min(8,unique(colData(values$dds_obj)[, values$color_by[1]]) %>% length()), "Dark2")
      if(is.null(values$color_by)) {
        #TODO message about color_by
        return(NULL)
      }
      # browser()
      barplot(sizeFactors(values$dds_obj),
              main = "Size factors ",
              col = col[as.integer(colData(values$dds_obj)[, values$color_by[1]])],
              cex.axis = 1.2, cex.names = 0.8, las = 3
      )
      
    })
    output$pairwise_plotUI <- renderUI({
      req(input$compute_pairwisecorr)
      
      plotOutput(ns("corrplot"), height = "1000px")
      # )
    })
    
    output$heatcorr_plotUI <- renderUI({
      input$compute_pairwisecorr
      # if(!input$compute_pairwisecorr) return()
      # TODO check that values$color_by is set and set message otherwise
      plotlyOutput(ns("heatcorr"))
    })
    output$pca_Dim1 <-renderUI({
      req(input$compute_pairwisecorr)
      selectizeInput(ns("pcaDim1"), label = "Which PCA on X",
                     choices = c(1:20), selected = 1, multiple = F)
    })
    output$pca_Dim2 <-renderUI({
      req(input$compute_pairwisecorr)
      selectizeInput(ns("pcaDim2"), label = "Which PCA on Y",
                     choices = c(1:20), selected = 2, multiple = F)
    })
    output$pca_plotUI <- renderUI({
      req(input$compute_pairwisecorr)
      plotlyOutput(ns("pcaPlot"))
    })
    output$pcaEV_plotUI <- renderUI({
      req(input$compute_pairwisecorr)
      plotlyOutput(ns("pcaEVPlot"))
    })
    output$pca_plotUI34 <- renderUI({
      req(input$compute_pairwisecorr)
      plotlyOutput(ns("pcaPlot34"))
    })
    output$pca_plotUI56 <- renderUI({
      req(input$compute_pairwisecorr)
      plotlyOutput(ns("pcaPlot56"))
    })
    output$sizeFactors_plotUI    <- renderUI({
      req(input$compute_pairwisecorr)
      plotOutput(ns("sizeFactorsPlot"))
    }) 
    output$geneHeatmap_plotUI <- renderUI({
      cat(file =stderr(), paste("geneHeatmap: ",input$compute_pairwisecorr, "\n"))
      if(input$compute_pairwisecorr<1) return()
      plotlyOutput(ns("geneHeatmap"))
    })
    output$geneHeatmap_genesUI <- renderUI({
      cat(file =stderr(), paste("geneHeatmap: ",input$compute_pairwisecorr, "\n"))
      if(input$compute_pairwisecorr<1) return()
      verbatimTextOutput(ns("geneHeatmapgenes"))
    })
    
    output$geneHeatmap <- renderPlotly({
      cat(file =stderr(), paste("geneHeatmap renderPlot: ", input$compute_pairwisecorr,"\n"))
      # We define how many genes we want to look at
      nGenesHeatmap <- 20
      
      # ordered by total count over all experiments per gene
      select <- order(rowMeans(counts(values$dds_obj, normalized = TRUE)),
                      decreasing = TRUE
      )[1:nGenesHeatmap]
      
      # we transform the data into log2 space
      nt <- normTransform(values$dds_obj) # defaults to log2(x+1)
      
      # and extract the count data from the DESeq2 object.
      heatmapcounts <- assay(nt)[select, ]
      # save(file = "~/scShinyHubDebug/WEIN.RData", list = c(ls(), ls(envir = globalenv())))
      # load(file="~/scShinyHubDebug/WEIN.RData")
      # if("symbol" %in% names(dds_obj)) {
      #   p <- plot_ma(values$res_obj,
      #                intgenes = values$genelist_ma$`Gene Symbol`,annotation_obj = values$annotation_obj,FDR = input$FDR)
      # }
      # cluster_cols builds the dendrogram on the top
      if (!is.null(values$annotation_obj)){
        rownames(heatmapcounts) <- values$annotation_obj$gene_name[
          match(rownames(heatmapcounts),
                rownames(values$annotation_obj))]
      }
      heatmaply(heatmapcounts, cluster_rows = T, show_rownames = TRUE, cluster_cols = T)
      # p = ggplotify::as.ggplot(p)
      # p$theme = list()
      # # browser()
      # cat(file =stderr(), paste("geneHeatmap renderPlot: ", input$compute_pairwisecorr,"done.\n"))
      # # ggsave(
      # #   filename = paste0("test.",input$compute_pairwisecorr, ".png"), plot = p)
      # p
    })
    
    output$geneHeatmapgenes <- renderPrint({
      cat(file =stderr(), paste("geneHeatmapgenes renderPlot: ", input$compute_pairwisecorr,"\n"))
      # We define how many genes we want to look at
      nGenesHeatmap <- 20
      
      # ordered by total count over all experiments per gene
      select <- order(rowMeans(counts(values$dds_obj, normalized = TRUE)),
                      decreasing = TRUE
      )[1:nGenesHeatmap]
      paste(rownames(values$dds_obj[select, ]), collapse = " ")
    })
    
    # overview on number of detected genes on different threshold types
    output$detected_genes <- renderPrint({
      t1 <- rowSums(counts(values$dds_obj))
      t2 <- rowMeans(counts(values$dds_obj,normalized=TRUE))
      
      thresh_rowsums <- input$threshold_rowsums
      thresh_rowmeans <- input$threshold_rowmeans
      abs_t1 <- sum(t1 > thresh_rowsums)
      rel_t1 <- 100 * mean(t1 > thresh_rowsums)
      abs_t2 <- sum(t2 > thresh_rowmeans)
      rel_t2 <- 100 * mean(t2 > thresh_rowmeans)
      
      cat("Number of detected genes:\n")
      cat(abs_t1,"genes have at least a sample with more than",thresh_rowsums,"counts\n")
      cat(paste0(round(rel_t1,3),"%"), "of the",nrow(values$dds_obj),
          "genes have at least a sample with more than",thresh_rowsums,"counts\n")
      cat(abs_t2,"genes have more than",thresh_rowmeans,"counts (normalized) on average\n")
      cat(paste0(round(rel_t2,3),"%"), "of the",nrow(values$dds_obj),
          "genes have more than",thresh_rowsums,"counts (normalized) on average\n")
      cat("Counts are ranging from", min(counts(values$dds_obj)),"to",max(counts(values$dds_obj)))
    })
    
    observeEvent(input$featfilt_dds,
                 {
                   t1 <- rowSums(counts(values$dds_obj))
                   t2 <- rowMeans(counts(values$dds_obj,normalized=TRUE))
                   
                   thresh_rowsums <- input$threshold_rowsums
                   thresh_rowmeans <- input$threshold_rowmeans
                   
                   if(input$filter_crit == "row sums") {
                     filt_dds <- values$dds_obj[t1 > thresh_rowsums, ]
                   } else {
                     filt_dds <- values$dds_obj[t2 > thresh_rowmeans, ]
                   }
                   
                   # TODO: see if re-estimation of size factors is required
                   filt_dds <- estimateSizeFactors(filt_dds)
                   
                   values$dds_obj <- filt_dds
                   
                 })
    
  })
}