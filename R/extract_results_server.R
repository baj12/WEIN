extract_results_server <- function(id, values, annoSpecies_df, exportPlots) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    design_factors <- reactive({
      # rev(attributes(terms.formula(design(values$dds_obj)))$term.labels)
      cat(file = stderr(), "design_factors triggered\n")
      resultsNames(values$dds_obj)
      
    })
    
    output$choose_fac <- renderUI({
      selectInput(ns("choose_expfac"),label = "Choose the experimental factor to build the contrast (numerator)",
                  choices = c("",design_factors()), selected = "", multiple = TRUE)
    })
    output$choose_fac2 <- renderUI({
      selectInput(ns("choose_expfac2"),label = "Choose the experimental factor to build the contrast (denominator)",
                  choices = c("",design_factors()), selected = "", multiple = TRUE)
    })
    
    
    output$runresults <- renderUI({
      shiny::validate(
        need(input$choose_expfac!="",
             "Select a factor for the contrast first")
      )
      
      actionButton(ns("button_runresults"),"Extract the results!", icon = icon("spinner"), class = "btn btn-success")
    })
    
    observeEvent(input$button_runresults, {
      choose_expfac2 <- input$choose_expfac2
      if (is.null(choose_expfac2)) 
        choose_expfac2 = character()
      resultsNames(values$dds_obj)
      choose_expfac <- input$choose_expfac
      # choose_expfac <- c("Intercept", "STIMULUS_LPS_vs_antiCD3CD28" )
      # choose_expfac2 = "STIMULUS_null_vs_antiCD3CD28"
      # browser()
      withProgress(message="Computing the results...",
                   detail = "DE table on its way!",
                   value = 0,{
                     # handling the experimental covariate correctly to extract the results...
                     # if(is.factor(colData(values$dds_obj)[,input$choose_expfac])) {
                     if(input$resu_ihw) {
                       values$res_obj <- results(values$dds_obj,
                                                 # contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2),
                                                 contrast = list(c(choose_expfac),c(choose_expfac2)),
                                                 independentFiltering = input$resu_indfil, 
                                                 alpha = values$FDR,
                                                 filterFun = ihw)
                       
                       # if(input$resu_lfcshrink) {
                       #   incProgress(amount = 0.15,detail = "Results extracted. Shrinking the logFC now...")
                       #   values$res_obj <- lfcShrink(values$dds_obj,
                       #                               # contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2),
                       #                               contrast = list(c(input$choose_expfac),c(choose_expfac2)),
                       #                               res = values$res_obj)
                       #   incProgress(amount = 0.8,detail = "logFC shrunken, adding annotation info...")
                       # } else {
                       incProgress(amount = 0.9,detail = "logFC left unshrunken, adding annotation info...")
                       # }
                     } else {
                       values$res_obj <- results(values$dds_obj,
                                                 # contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2),
                                                 contrast = list(c(choose_expfac),c(choose_expfac2)),
                                                 independentFiltering = input$resu_indfil, 
                                                 alpha = values$FDR)
                       # if(input$resu_lfcshrink) {
                       #   incProgress(amount = 0.15,detail = "Results extracted. Shrinking the logFC now...")
                       #   
                       #   coef = c("Intercept", "STIMULUS_LPS_vs_antiCD3CD28")
                       #   values$res_obj <- lfcShrink(values$dds_obj,
                       #                               # contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2),
                       #                               contrast = list(c(input$choose_expfac),c(choose_expfac2)),
                       #                               res = values$res_obj,
                       #                               type = "ashr")
                       #   incProgress(amount = 0.8,detail = "logFC shrunken, adding annotation info...")
                       # } else {
                       incProgress(amount = 0.9,detail = "logFC left unshrunken, adding annotation info...")
                       # }
                     }
                     # }
                     # should not happen as all are factors,
                     # BJ commenting out
                     # if(class(colData(values$dds_obj)[,input$choose_expfac]) %in% c("integer","numeric"))
                     #   values$res_obj <- results(values$dds_obj,name = input$choose_expfac,
                     #                             independentFiltering = input$resu_indfil, 
                     #                             alpha = values$FDR
                     #                             # , addMLE = input$resu_lfcshrink
                     #   )
                     
                     # adding info from the annotation
                     if(!is.null(values$annotation_obj))
                       values$res_obj$symbol <- values$annotation_obj$gene_name[
                         match(rownames(values$res_obj),
                               rownames(values$annotation_obj))]
                   })
    })
    
    output$diyres_summary <- renderPrint({
      c1  = input$choose_expfac
      c2 = input$choose_expfac2
      # browser()
      shiny::validate(
        # need(input$choose_expfac!="" & input$fac1_c1 != "" & input$fac1_c2 != "" & input$fac1_c1 != input$fac1_c2 ,
        need(input$choose_expfac != ""  ,
             "Please select a coefficient to build the contrast"
        )
      )
      shiny::validate(
        need(!is.null(values$res_obj), "Parameters selected, please compute the results first")
      )
      # summary(results(values$dds_obj,contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2)))
      summary(values$res_obj,alpha = values$FDR)
    })
    
    output$printdds <- renderPrint({
      shiny::validate(
        need(!is.null(values$dds_obj),
             "Please provide a count matrix/dds object"
        )
      )
      
      values$dds_obj
      design(values$dds_obj)
    })
    
    output$printres <- renderPrint({
      shiny::validate(
        need(!is.null(values$res_obj),
             "Please provide a DESeqResults object"
        )
      )
      print(sub(".*p-value: (.*)","\\1",mcols(values$res_obj, use.names=TRUE)["pvalue","description"]))
      summary(values$res_obj,alpha = values$FDR) # use fdr shiny widget
    })
    
    
    output$store_result <- renderUI({
      if(is.null(values$res_obj))
        return(NULL)
      actionButton(ns("button_store_result"), "Store current results",class = "btn btn-primary")
    })
    
    observeEvent(input$button_store_result,
                 {
                   values$stored_res <- values$res_obj
                   # this is in such a way to store & compare later if some parameters are edited
                 })
    
    output$table_res <- DT::renderDataTable(server=T,{
      if(is.null(values$res_obj))
        return(NULL)
      mydf <- as.data.frame(values$res_obj[order(values$res_obj$padj),])#[1:500,]
      rownames(mydf) <- createLinkENS(rownames(mydf),species = annoSpecies_df$ensembl_db[match(values$cur_species,annoSpecies_df$species)]) ## TODO: check what are the species from ensembl and
      ## TODO: add a check to see if wanted?
      mydf$symbol <- createLinkGeneSymbol(mydf$symbol)
      # browser()
      datatable(mydf, extensions = 'Buttons', 
                options = list(dom = 'lfrtipB', buttons = c('csv')),
                escape = FALSE, filter = list(position = 'top', clear = FALSE))%>%
        formatRound(columns=c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'), digits=3)
    })
    
    # server resu diagnostics --------------------------------------------------------
    output$pvals_hist <- renderPlot({
      shiny::validate(
        need(!is.null(values$res_obj),message = "")
      )
      
      res_df <- as.data.frame(values$res_obj)
      res_df <- dplyr::filter(res_df, !is.na(pvalue))
      p <- ggplot(res_df, aes_string("pvalue")) +
        geom_histogram(binwidth = 0.01, boundary = 0) + theme_bw()
      
      # for visual estimation of the false discovery proportion in the first bin
      alpha <- binw <- values$FDR
      pi0 <- 2*mean(res_df$pvalue > 0.5)
      p <- p + geom_hline(yintercept = pi0 * binw * nrow(res_df), col = "steelblue") + 
        geom_vline(xintercept = alpha, col = "red")
      
      p <- p + ggtitle(
        label = "p-value histogram",
        subtitle = paste0(
          "Expected nulls = ", pi0 * binw * nrow(res_df), 
          " - #elements in the selected bins = ", sum(res_df$pvalue < alpha)
        ))
      
      exportPlots$plot_pvals_hist <- p
      p
      
    })
    
    output$pvals_hist_strat <- renderPlot({
      shiny::validate(
        need(!is.null(values$res_obj),message = "")
      )
      
      res_df <- as.data.frame(values$res_obj)
      res_df <- dplyr::filter(res_df, !is.na(pvalue))
      
      res_df <- mutate(
        res_df, 
        stratum = cut(baseMean, include.lowest = TRUE, 
                      breaks = signif(quantile(baseMean, probs = seq(0,1, length.out = 10)),2)))
      
      p <- ggplot(res_df, aes_string("pvalue")) +
        geom_histogram(binwidth = 0.01, boundary = 0) + 
        facet_wrap(~stratum) + 
        theme_bw()
      
      p <- p + ggtitle(
        label = "p-value histogram",
        subtitle = "stratified on the different value classes of mean expression values")
      
      exportPlots$plot_pvals_hist_strat <- p
      p
    })
    
    output$pvals_ss <- renderPlot({
      shiny::validate(
        need(!is.null(values$res_obj),message = "")
      )
      
      res_df <- as.data.frame(values$res_obj)
      res_df <- dplyr::filter(res_df, !is.na(pvalue))
      
      phi <- values$FDR
      res_df <- mutate(res_df, rank = rank(pvalue))
      m <- nrow(res_df)
      
      p <- ggplot(filter(res_df, rank <= 6000), 
                  aes_string(x = "rank", y = "pvalue")) + 
        geom_line() + 
        geom_abline(slope = phi/m, col = "red") + 
        theme_bw()
      
      p <- p + ggtitle(
        label = "Schweder-Spjotvoll plot",
        subtitle = paste0(
          "Intersection point at rank ", with(arrange(res_df,rank), last(which(pvalue <= phi * rank / m))))
      )
      exportPlots$plot_pvals_ss <- p
      p
    })
    # DE genes lists ----------------------------------------------------------
    values$genelistUP <- reactive({
      listUP <- tryCatch({
        # browser()
        res_tbl <- deseqresult2DEgenes(values$res_obj, FDR = values$FDR)
        if(nrow(res_tbl)<1) return(NULL)
        res_tbl_UP <- res_tbl[res_tbl$log2FoldChange > 0 & !is.na(res_tbl$padj),]
        # res_tbl_DOWN <- res_tbl[res_tbl$log2FoldChange < 0 & !is.na(res_tbl$padj),]
        
        if("symbol" %in% colnames(values$res_obj)) { 
          if(!is.null(values$annotation_obj)) {
            res_tbl_UP$symbol <- values$annotation_obj$gene_name[
              match(res_tbl_UP$id,
                    rownames(values$annotation_obj))]
            listUP <- res_tbl_UP$symbol
          } else {
            listUP <- NULL
          }
        } else {
          listUP <- res_tbl_UP$symbol
        }
        listUP},
        error=function(e)cat(file = stderr(), paste("genelistUP error ; ", e))
      )
      return(listUP)
    })
    
    values$genelistDOWN <- reactive({
      # browser()
      res_tbl <- deseqresult2DEgenes(values$res_obj, FDR = values$FDR)
      if(nrow(res_tbl)<1) return(NULL)
      # res_tbl_UP <- res_tbl[res_tbl$log2FoldChange > 0 & !is.na(res_tbl$padj),]
      res_tbl_DOWN <- res_tbl[res_tbl$log2FoldChange < 0 & !is.na(res_tbl$padj),]
      
      if("symbol" %in% colnames(values$res_obj)) { 
        if(!is.null(values$annotation_obj)) {
          res_tbl_DOWN$symbol <- values$annotation_obj$gene_name[
            match(res_tbl_DOWN$id,
                  rownames(values$annotation_obj))]
          listDOWN <- res_tbl_DOWN$symbol
        } else {
          listDOWN <- NULL
        }
      } else {
        listDOWN <- res_tbl_DOWN$symbol
      }
      return(listDOWN)
    })
    
    values$genelistUPDOWN <- reactive({
      # browser()
      res_tbl <- deseqresult2DEgenes(values$res_obj, FDR = values$FDR)
      if(nrow(res_tbl)<1) return(NULL)
      if("symbol" %in% colnames(values$res_obj)) { 
        if(!is.null(values$annotation_obj)) {
          res_tbl$symbol <- values$annotation_obj$gene_name[
            match(res_tbl$id,
                  rownames(values$annotation_obj))]
          listUPDOWN <- res_tbl$symbol
        } else {
          listUPDOWN <- NULL
        }
      } else {
        listUPDOWN <- res_tbl$symbol
      }
      return(listUPDOWN)
    })
    
    
    output$logfc_hist <- renderPlot({
      shiny::validate(
        need(!is.null(values$res_obj),message = "")
      )
      res_df <- as.data.frame(values$res_obj)
      res_df <- dplyr::filter(res_df, !is.na(pvalue))
      
      p <- ggplot(res_df, aes_string("log2FoldChange")) +
        geom_histogram(binwidth = 0.1) + theme_bw()
      
      p <- p + ggtitle(
        "Histogram of the log2 fold changes"
      )
      
      exportPlots$plot_logfc_hist <- p
      p
    })
    output$downloadTblResu <- downloadHandler(
      filename = function() {
        "table_results.csv"
      },
      content = function(file) {
        mydf <- as.data.frame(values$res_obj[order(values$res_obj$padj),])
        write.csv(mydf, file)
      }
    )
    output$download_plot_pvals_hist <- downloadHandler(filename = function() {
      input$filename_plot_pvals_hist
    }, content = function(file) {
      ggsave(file, exportPlots$plot_pvals_hist, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    output$download_plot_logfc_hist <- downloadHandler(filename = function() {
      input$filename_plot_logfc_hist
    }, content = function(file) {
      ggsave(file, exportPlots$plot_logfc_hist, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    
    
  })
}
  
