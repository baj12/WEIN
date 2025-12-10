signature_explorer_server <- function(id, values, annoSpecies_df, exportPlots) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    # server signature explorer ------------------------------------------------------
    output$sig_ui_gmtin <- renderUI({
      fileInput(ns("sig_gmtin"),"gmt input file")
    })
    
    loaded_gmt <- reactive({
      if (is.null(input$sig_gmtin))
        return(NULL)
      mysigs <- read_gmt(input$sig_gmtin$datapath)
      return(mysigs)
    })
    
    observeEvent(input$sig_gmtin,
                 {
                   values$gene_signatures <- loaded_gmt()
                 })
    
    output$sig_ui_nrsigs <- renderUI({
      if(!is.null(values$gene_signatures))
        return(valueBox("Gene signatures",
                        paste0(length(values$gene_signatures), " gene signatures"),
                        icon = icon("list"),
                        color = "green",width = NULL))
      else
        return(valueBox("Gene signatures",
                        "yet to be loaded",
                        icon = icon("list"),
                        color = "red",width = NULL))
    })
    
    observeEvent(input$sig_button_computevst,
                 {
                   withProgress(message="Computing the variance stabilized transformed data...",
                                detail = "This step can take a little while",
                                value = 0,{
                                  values$vst_obj <- vst(values$dds_obj)
                                })
                 })
    
    output$sig_ui_selectsig <- renderUI({
      if(!is.null(values$gene_signatures))
        return(selectizeInput(ns("sig_selectsig"), label = "Select the gene signature",
                              choices = NULL, selected = NULL, multiple = FALSE))
      else
        return(NULL)
    })
    
    observe({
      updateSelectizeInput(session = session, inputId = 'sig_selectsig', choices = c(Choose = '', names(values$gene_signatures)), server = TRUE)
    })
    
    output$sig_sigmembers <- renderPrint({
      values$gene_signatures[[input$sig_selectsig]]
    })
    
    output$sig_ui_annocoldata <- renderUI({
      if(!is.null(values$dds_obj))
        return(selectizeInput(ns("sig_annocoldata"), label = "Select the colData to decorate",
                              choices = names(colData(values$dds_obj)),
                              selected = NULL, multiple = TRUE))
      else
        return(NULL)
    })
    
    
    output$sig_ui_id_data <- renderUI({
      if (is.null(values$dds_obj)) #
        return(NULL)
      validate(
        need(!is.null(values$cur_species), message = "Please specify the species in the Data Setup panel")
      )
      
      std_choices <- c("SYMBOL", "ENSEMBL","ENTREZID","REFSEQ")
      if (values$cur_species!=""){
        annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==values$cur_species]
        pkg_choices <- keytypes(get(annopkg))
        std_choices <- union(std_choices, pkg_choices)
      }
      selectInput(ns("sig_id_data"), "select the id type in your dds data", choices=std_choices)
    })
    
    output$sig_ui_id_sigs <- renderUI({
      if (is.null(values$gene_signatures)) #
        return(NULL)
      validate(
        need(!is.null(values$cur_species), message = "Please specify the species in the Data Setup panel")
      )
      
      std_choices <- c("SYMBOL", "ENSEMBL","ENTREZID","REFSEQ")
      if (values$cur_species!=""){
        annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==values$cur_species]
        pkg_choices <- keytypes(get(annopkg))
        std_choices <- union(std_choices, pkg_choices)
      }
      selectInput(ns("sig_id_sigs"), "select the id type in your signatures", choices=std_choices)
    })
    
    available_orgdb <- rownames(installed.packages())[
      grep(pattern = "^org.*db$",rownames(installed.packages()))]
    
    output$sig_ui_orgdbpkg <- renderUI({
      
      suggested_orgdb <- tryCatch(
        annoSpecies_df$pkg[annoSpecies_df$species==values$cur_species],
        error = function(e){return("")})
      selectInput(ns("sig_orgdbpkg"), "Select the organism package for matching", 
                  choices=c("",available_orgdb),selected = suggested_orgdb)
    })
    
    observeEvent(values$cur_species,
                 {
                   suggested_orgdb <- annoSpecies_df$pkg[annoSpecies_df$species==values$cur_species]
                   if(suggested_orgdb %in% available_orgdb)
                     updateSelectInput(session, inputId = "sig_orgdbpkg", selected = suggested_orgdb)
                 })
    
    observeEvent(input$sig_convert_setup,
                 {
                   require(input$sig_orgdbpkg,character.only=TRUE)
                   withProgress(message="Matching the identifiers",
                                detail = "Locating package", value = 0,{
                                  
                                  incProgress(0.1,detail = "Matching identifiers")
                                  
                                  x <- get(input$sig_orgdbpkg)
                                  # browser()
                                  if(any(is.null(input$sig_id_sigs), is.null(input$sig_id_data))){
                                    cat(file = stderr(), paste("\n\ndid you specifvy the annotations?\n\n"))
                                    return(NULL)
                                  }
                                  if (input$sig_id_sigs == "SYMBOL" & input$sig_id_data == "SYMBOL") {
                                    anno_vec = rownames(values$dds_obj)
                                    names(anno_vec) = rownames(values$dds_obj)
                                    values$anno_vec = anno_vec
                                  } else {
                                    values$anno_vec <- mapIds(x, rownames(values$dds_obj),
                                                              column = input$sig_id_sigs,
                                                              keytype = input$sig_id_data)
                                  }
                                })
                 })
    
    output$sig_convcheck <- renderPrint({
      head(values$anno_vec)
    })
    
    output$sig_heat <- renderPlotly({
      validate(
        need(!is.null(values$gene_signatures), message = "Please provide some gene signatures in gmt format"),
        need(!is.null(values$vst_obj), message = "Compute the vst transformed data"),
        need(!is.null(values$anno_vec), message = "Setup the conversion between data ids and signature ids"),
        need((!is.null(values$res_obj) | !input$sig_useDEonly),
             message = "Please compute the results first if you want to subset to DE genes only"),
        need(input$sig_selectsig!="", message = "Select a signature")
      )
      
      print(
        sig_heatmap(
          values$vst_obj,
          my_signature = values$gene_signatures[[input$sig_selectsig]],
          res_data = values$res_obj,
          FDR = values$FDR,
          de_only = input$sig_useDEonly,
          annovec = values$anno_vec,
          # anno_colData = colData(values$vst_obj)[,input$sig_annocoldata, drop = FALSE],
          title = names(values$gene_signatures)[match(input$sig_selectsig,names(values$gene_signatures))],
          cluster_rows = input$sig_clusterrows,
          cluster_cols = input$sig_clustercols,
          center_mean = input$sig_centermean,
          scale_row = input$sig_scalerow
        ))
      
    })
    
    output$sig_heat_genes <- renderPrint({
      validate(
        need(!is.null(values$gene_signatures), message = "Please provide some gene signatures in gmt format"),
        need(!is.null(values$vst_obj), message = "Compute the vst transformed data"),
        need(!is.null(values$anno_vec), message = "Setup the conversion between data ids and signature ids"),
        need((!is.null(values$res_obj) | !input$sig_useDEonly),
             message = "Please compute the results first if you want to subset to DE genes only"),
        need(input$sig_selectsig!="", message = "Select a signature")
      )
      annovec = values$anno_vec
      mydata <- assay(values$vst_obj)
      
      my_signature = values$gene_signatures[[input$sig_selectsig]]
      # save(file = "~/SCHNAPPsDebug/ideal.sig_heatmap.RData", list = c(ls()))
      # load("~/SCHNAPPsDebug/ideal.sig_heatmap.RData")
      signature_original_ids <- names(annovec)[match(my_signature,annovec)]
      
      sig_to_keep <- (signature_original_ids %in% rownames(mydata))#
      my_signature <- my_signature[sig_to_keep]
      signature_original_ids <- signature_original_ids[sig_to_keep]
      
      mydata_sig <- mydata[signature_original_ids,]
      
      # to avoid problems later, remove the ones non-expressed and with variance = 0
      to_remove <- apply(mydata_sig, 1, var) == 0
      mydata_sig <- mydata_sig[!to_remove,]
      mydata_sig <- mydata[signature_original_ids,]
      
      paste(rownames(mydata_sig), collapse = " ")
    })
    
    
  })
}