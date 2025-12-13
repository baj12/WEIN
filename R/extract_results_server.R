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
                  choices = c("",design_factors()), selected = ifelse(is.null(values$choose_expfac), "", values$choose_expfac), multiple = TRUE)
    })
    output$choose_fac2 <- renderUI({
      selectInput(ns("choose_expfac2"),label = "Choose the experimental factor to build the contrast (denominator)",
                  choices = c("",design_factors()), selected = ifelse(is.null(values$choose_expfac2), "", values$choose_expfac2), multiple = TRUE)
    })
    
    
    output$runresults <- renderUI({
      shiny::validate(
        need(input$choose_expfac!="",
             "Please select a factor for the contrast first")
      )
      
      actionButton(ns("button_runresults"),"Extract the results!", icon = icon("spinner"), class = "btn btn-success")
    })
    
    # Store the selected experimental factors in values object
    observeEvent(input$choose_expfac, {
      values$choose_expfac <- input$choose_expfac
    })
    
    observeEvent(input$choose_expfac2, {
      values$choose_expfac2 <- input$choose_expfac2
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
                       
                       incProgress(amount = 0.9,detail = "logFC left unshrunken, adding annotation info...")
                       # }
                     } else {
                       values$res_obj <- results(values$dds_obj,
                                                 # contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2),
                                                 contrast = list(c(choose_expfac),c(choose_expfac2)),
                                                 independentFiltering = input$resu_indfil, 
                                                 alpha = values$FDR)
                       incProgress(amount = 0.9,detail = "logFC left unshrunken, adding annotation info...")
                       # }
                     }

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
             "Please select a coefficient to build the contrast first"
        )
      )
      shiny::validate(
        need(!is.null(values$res_obj), "Parameters selected, please compute the results first by clicking the 'Compute Results' button")
      )
      # summary(results(values$dds_obj,contrast = c(input$choose_expfac, input$fac1_c1, input$fac1_c2)))
      summary(values$res_obj,alpha = values$FDR)
    })
    
    output$printdds <- renderPrint({
      shiny::validate(
        need(!is.null(values$dds_obj),
             "Please provide a count matrix or DESeqDataSet object to view results"
        )
      )
      
      values$dds_obj
      design(values$dds_obj)
    })
    
    output$printres <- renderPrint({
      shiny::validate(
        need(!is.null(values$res_obj),
             "Please provide a DESeqResults object or compute results first"
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
      # Check if species is selected and available in Ensembl
      if (!is.null(values$cur_species) && values$cur_species %in% annoSpecies_df$species) {
        ensembl_species <- annoSpecies_df$ensembl_db[match(values$cur_species, annoSpecies_df$species)]
        if (!is.na(ensembl_species) && ensembl_species != "") {
          rownames(mydf) <- createLinkENS(rownames(mydf), species = ensembl_species)
        }
      }
      mydf$symbol <- createLinkGeneSymbol(mydf$symbol)
      # browser()
      datatable(mydf, extensions = 'Buttons',
                options = list(dom = 'lfrtipB', buttons = c('csv'),
                              pageLength = 10, deferRender = TRUE),
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
      p <- ggplot(res_df, aes(x=.data[["pvalue"]])) +
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
      
      p <- ggplot(res_df, aes(x=.data[["pvalue"]])) +
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
                  aes(x = .data[["rank"]], y = .data[["pvalue"]])) +
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
        # Use the new utility function
        gene_lists <- generate_gene_lists(values$res_obj, values$FDR, values$annotation_obj)
        gene_lists$UP},
        error=function(e)cat(file = stderr(), paste("genelistUP error ; ", e))
      )
      return(listUP)
    })
    
    values$genelistDOWN <- reactive({
      # browser()
      # Use the new utility function
      gene_lists <- generate_gene_lists(values$res_obj, values$FDR, values$annotation_obj)
      gene_lists$DOWN
    })
    
    values$genelistUPDOWN <- reactive({
      # browser()
      # Use the new utility function
      gene_lists <- generate_gene_lists(values$res_obj, values$FDR, values$annotation_obj)
      gene_lists$UPDOWN
    })
    
    
    output$logfc_hist <- renderPlot({
      shiny::validate(
        need(!is.null(values$res_obj),message = "")
      )
      res_df <- as.data.frame(values$res_obj)
      res_df <- dplyr::filter(res_df, !is.na(pvalue))
      
      p <- ggplot(res_df, aes(x=.data[["log2FoldChange"]])) +
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
  
