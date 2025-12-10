functional_analysis_server <- function(id, values, annoSpecies_df, exportPlots) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    # Helper function to create gene list reactive and observer
    create_gene_list_handler <- function(list_num) {
      gl_name <- paste0("gl", list_num)
      genelist_name <- paste0("genelist", list_num)
      
      # Create reactive
      gl_reactive <- reactive({
        if (is.null(input[[gl_name]])) {
          return(data.frame())
        } else {
          gl <- read1stCol(input[[gl_name]]$datapath, values$dds_obj)
          if(is.null(gl)) return(NULL)
          return(gl)
        }
      })
      
      # Create observer
      observeEvent(input[[gl_name]], {
        gl <- gl_reactive()
        if(is.null(gl) || nrow(gl) < 1) {
          values[[genelist_name]] <- data.frame()
          return(NULL)
        }
        mydf <- as.data.frame(gl, stringsAsFactors = FALSE)
        names(mydf) <- "Gene Symbol"
        values[[genelist_name]] <- mydf
      })
      
      return(gl_reactive)
    }
    
    # Create all gene list handlers
    gl_reactives <- lapply(1:4, create_gene_list_handler)
    names(gl_reactives) <- paste0("gl", 1:4)
    
    # Helper function to get gene list
    get_gene_list <- function(list_name) {
      switch(list_name,
             "UP" = values$genelistUP(),
             "DOWN" = values$genelistDOWN(),
             "UPDOWN" = values$genelistUPDOWN(),
             "LIST1" = as.character(values$genelist1$`Gene Symbol`),
             "LIST2" = as.character(values$genelist2$`Gene Symbol`),
             "LIST3" = as.character(values$genelist3$`Gene Symbol`),
             "LIST4" = as.character(values$genelist4$`Gene Symbol`)
      )
    }
    
    # Generic enrichment function
    perform_enrichment <- function(list_name, method = "goana") {
      # Validation
      if (is.null(values$cur_species) || values$cur_species == "") {
        showNotification("Please specify the species in the Data Setup panel", type = "warning")
        return(NULL)
      }
      
      genelist <- get_gene_list(list_name)
      if (is.null(genelist) || length(genelist) < 1) {
        showNotification(paste("The", list_name, "gene list is empty"), type = "warning")
        return(NULL)
      }
      
      tryCatch({
        organism <- annoSpecies_df[values$cur_species,]$species_short
        backgroundgenes <- rownames(values$dds_obj)[rowSums(counts(values$dds_obj)) > 0]
        annopkg <- annoSpecies_df[values$cur_species,]$pkg
        
        if (!require(annopkg, character.only = TRUE)) {
          stop("The package ", annopkg, " is not installed/available.")
        }
        
        if (method == "goana") {
          return(perform_goana_enrichment(genelist, backgroundgenes, organism, annopkg, list_name))
        } else if (method == "goseq") {
          return(perform_goseq_enrichment(genelist, backgroundgenes, annopkg, list_name))
        } else if (method == "topgo") {
          return(perform_topgo_enrichment(genelist, backgroundgenes, annopkg, list_name))
        }
      }, error = function(e) {
        showNotification(paste("Error during enrichment:", e$message), type = "warning")
        return(NULL)
      })
    }
    
    # GOANA enrichment
    perform_goana_enrichment <- function(genelist, backgroundgenes, organism, annopkg, list_name) {
      listGenesEntrez <- AnnotationDbi::mapIds(
        eval(parse(text = annopkg)), 
        keys = genelist,
        column = "ENTREZID", 
        keytype = "SYMBOL"
      )
      listBackgroundEntrez <- AnnotationDbi::mapIds(
        eval(parse(text = annopkg)), 
        keys = backgroundgenes,
        column = "ENTREZID", 
        keytype = values$cur_type
      )
      incProgress(0.1, detail = "IDs mapped")
      
      result <- limma::topGO(
        limma::goana(listGenesEntrez, listBackgroundEntrez, species = organism),
        ontology = input$go_cats[1],
        number = 200
      )
      
      incProgress(0.7, detail = "adding gene names to GO terms")
      go_ids <- rownames(result)
      allegs_list <- lapply(go_ids, function(arg) 
        AnnotationDbi::get(arg, get(paste0("org.", organism, ".egGO2ALLEGS"))))
      genes_list <- lapply(allegs_list, function(arg) 
        unlist(AnnotationDbi::mget(arg, get(paste0("org.", organism, ".egSYMBOL")))))
      DEgenes_list <- lapply(genes_list, function(arg) intersect(arg, genelist))
      
      result$genes <- unlist(lapply(DEgenes_list, function(arg) paste(arg, collapse = ",")))
      return(result)
    }
    
    # GOSEQ enrichment
    perform_goseq_enrichment <- function(genelist, backgroundgenes, annopkg, list_name) {
      assayed.genes.ids <- rownames(values$dds_obj)
      assayed.genes <- mapIds(
        get(annopkg),
        keys = assayed.genes.ids,
        column = "ENSEMBL",
        keytype = values$cur_type,
        multiVals = "first"
      )
      
      de.genes.ids <- mapIds(
        get(annopkg),
        keys = genelist,
        column = "ENSEMBL",
        keytype = "SYMBOL",
        multiVals = "first"
      )
      
      incProgress(0.1, detail = "IDs mapped")
      
      genesin <- unique(unname(de.genes.ids[!is.na(de.genes.ids)]))
      assayedGenes <- unique(unname(assayed.genes[!is.na(assayed.genes)]))
      
      length_data <- get(paste0(annoSpecies_df[values$cur_species,]$goseq_short, ".ensGene.LENGTH"))
      genesin <- genesin[genesin %in% length_data$Gene]
      assayedGenes <- assayedGenes[assayedGenes %in% length_data$Gene]
      
      result <- goseqTable(
        de.genes = genesin,
        assayed.genes = assayedGenes,
        genome = annoSpecies_df[values$cur_species,]$goseq_short,
        id = "ensGene",
        testCats = paste0("GO:", input$go_cats),
        FDR_GO_cutoff = 1,
        nTop = 200,
        addGeneToTerms = TRUE,
        orgDbPkg = annopkg
      )
      
      incProgress(0.89)
      return(result)
    }
    
    # TOPGO enrichment
    perform_topgo_enrichment <- function(genelist, backgroundgenes, annopkg, list_name) {
      bg_ids <- backgroundgenes[rowSums(counts(values$dds_obj)[backgroundgenes,]) > 0]
      bg_symbols <- if (values$cur_type == "SYMBOL") {
        bg_ids
      } else {
        mapIds(get(annopkg), keys = bg_ids, column = "SYMBOL", 
               keytype = values$cur_type, multiVals = "first")
      }
      
      incProgress(0.1, detail = "IDs mapped")
      
      if (length(unique(as.integer(bg_symbols %in% genelist))) == 1) {
        showNotification("All genes are the same", type = "warning")
        return(NULL)
      }
      
      result <- pcaExplorer::topGOtable(
        genelist, 
        bg_symbols,
        ontology = input$go_cats[1],
        mapping = annopkg,
        geneID = "symbol",
        addGeneToTerms = TRUE
      )
      
      incProgress(0.89)
      return(result)
    }
    
    # Create observers for all enrichment buttons
    list_types <- c("UP", "DOWN", "UPDOWN", paste0("LIST", 1:4))
    methods <- c("goana" = "", "goseq" = "_goseq", "topgo" = "_topgo")
    
    for (list_type in list_types) {
      for (method_name in names(methods)) {
        method_suffix <- methods[method_name]
        button_id <- paste0("button_enr", list_type, method_suffix)
        value_name <- paste0(
          ifelse(method_name == "goana", "gse_", paste0(method_name, "_")),
          tolower(list_type),
          ifelse(method_name == "goana", "", "")
        )
        if(method_name == "topgo") value_name <- paste0("topgo_", tolower(list_type))
        if(method_name == "goseq") value_name <- paste0("gse_", tolower(list_type), "_goseq")
        if(method_name == "goana") value_name <- paste0("gse_", tolower(list_type))
        
        local({
          lt <- list_type
          mn <- method_name
          vn <- value_name
          btn <- button_id
          
          observeEvent(input[[btn]], {
            msg <- paste(toupper(mn), "- Performing Gene Set Enrichment on", lt, "genes...")
            withProgress(message = msg, value = 0, {
              values[[vn]] <- perform_enrichment(lt, mn)
            })
          })
        })
      }
    }
    
    # Generic function to create UI output
    create_ui_output <- function(ns, value_name, title, output_name) {
      renderUI({
        if(is.null(values[[value_name]])) return(NULL)
        tagList(
          h4(title),
          DT::dataTableOutput(ns(output_name))
        )
      })
    }
    
    # Generic function to create DataTable output
    create_dt_output <- function(ns, value_name, link_column = "rownames", column_type = "GO") {
      DT::renderDataTable({
        if(is.null(values[[value_name]])) return(NULL)
        mytbl <- values[[value_name]]
        
        if (link_column == "rownames") {
          rownames(mytbl) <- createLinkGO(rownames(mytbl))
        } else {
          mytbl[[link_column]] <- createLinkGO(mytbl[[link_column]])
        }
        
        mytbl
      }, escape = FALSE, options = list(scrollX = TRUE))
    }
    
    # Generate all UI and DataTable outputs
    for (list_type in list_types) {
      for (method_name in names(methods)) {
        method_suffix <- methods[method_name]
        
        if(method_name == "topgo") {
          value_name <- paste0("topgo_", tolower(list_type))
          ui_name <- paste0("ui_DT_gse_", tolower(list_type), "_topgo")
          dt_name <- paste0("DT_gse_", tolower(list_type), "_topgo")
          title <- paste("topGO table -", list_type)
          link_col <- "GO.ID"
        } else if(method_name == "goseq") {
          value_name <- paste0("gse_", tolower(list_type), "_goseq")
          ui_name <- paste0("ui_DT_gse_", tolower(list_type), "_goseq")
          dt_name <- paste0("DT_gse_", tolower(list_type), "_goseq")
          title <- paste("goseq table -", list_type)
          link_col <- "category"
        } else {
          value_name <- paste0("gse_", tolower(list_type))
          ui_name <- paste0("ui_DT_gse_", tolower(list_type))
          dt_name <- paste0("DT_gse_", tolower(list_type))
          title <- paste("goana table -", list_type)
          link_col <- "rownames"
        }
        
        local({
          vn <- value_name
          uin <- ui_name
          dtn <- dt_name
          ttl <- title
          lc <- link_col
          
          output[[uin]] <- create_ui_output(ns, vn, ttl, dtn)
          output[[dtn]] <- create_dt_output(ns, vn, lc)
        })
      }
    }
    
    
    
    output$debuglists <- renderText({
      # length(gll_nonempty)
      # length(gll())
      # lapply(gll(),length)
      gll <- gll()
      txt = ""
      # save(file = "~/SCHNAPPsDebug/WEIN.RData", list = ls())
      # load("~/SCHNAPPsDebug/WEIN.RData")
      for (li in 1:length(gll)) {
        tx = paste(gll[[li]], sep = " ",collapse = " ")
        txt = paste(c(txt, names(gll)[li],tx), collapse = "\n")
      }
      txt
    })
    
    ## list of gene lists
    gll <- reactive({
      # browser()
      mylist <- list(listUP = values$genelistUP(),
                     listDOWN = values$genelistDOWN(),
                     listUPDOWN = values$genelistUPDOWN(),
                     list1 = as.character(values$genelist1$`Gene Symbol`),
                     list2 = as.character(values$genelist2$`Gene Symbol`),
                     list3 = as.character(values$genelist3$`Gene Symbol`),
                     list4 = as.character(values$genelist4$`Gene Symbol`))
      
      
      
      gll_nonempty <- mylist[!sapply(mylist,is.null)]
      
      # plus, add toggles to selectively keep only some lists?
      
      lists_tokeep <- names(mylist)[which(c(input$toggle_up,
                                            input$toggle_down,
                                            input$toggle_updown,
                                            input$toggle_list1,
                                            input$toggle_list2,
                                            input$toggle_list3,
                                            input$toggle_list4))]
      gll_final <- gll_nonempty[match(lists_tokeep,names(gll_nonempty))]
      return(gll_final) 
    })
    # Helper function to create GO term heatmap
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
    
    # Generate all heatmap outputs in a loop
    heatmap_types <- list(
      up = "up",
      down = "down", 
      updown = "updown",
      list1 = "list1",
      list2 = "list2",
      list3 = "list3",
      list4 = "list4"
    )
    
    for (hm_name in names(heatmap_types)) {
      local({
        hm <- hm_name
        topgo_name <- heatmap_types[[hm]]
        dt_input_name <- paste0("DT_gse_", hm, "_topgo_rows_selected")
        topgo_value_name <- paste0("topgo_", topgo_name)
        output_name <- paste0("goterm_heatmap_", hm, "_topgo")
        
        output[[output_name]] <- renderUI({
          s <- input[[dt_input_name]]
          
          # Show message if no row selected
          if(length(s) == 0) {
            return(div(
              style = "text-align: center; padding: 50px; color: #999; font-size: 16px;",
              icon("hand-pointer"),
              p("Please select a row from the table above to display the heatmap")
            ))
          }
          
          topgo_data <- values[[topgo_value_name]]
          if(is.null(topgo_data)) {
            return(div(
              style = "text-align: center; padding: 50px; color: #999; font-size: 16px;",
              p("No enrichment data available")
            ))
          }
          
          # Wrap plotlyOutput with jqui_resizable INSIDE renderUI
          shinyjqui::jqui_resizable(
            plotlyOutput(ns(paste0(output_name, "_plot")), height = "600px")
          )
        })
        
        # Separate renderPlotly for the actual heatmap
        output[[paste0(output_name, "_plot")]] <- renderPlotly({
          s <- input[[dt_input_name]]
          if(length(s) == 0) return(NULL)
          
          topgo_data <- values[[topgo_value_name]]
          if(is.null(topgo_data)) return(NULL)
          
          create_goterm_heatmap(topgo_data, s, values, annoSpecies_df)
        })
      })
    }    
    output$debugTable <- DT::renderDataTable(server =  TRUE,{
      gll <- gll()
      txt = ""
      # save(file = "~/SCHNAPPsDebug/WEIN.RData", list = ls())
      # load("~/SCHNAPPsDebug/WEIN.RData")
      upGll = UpSetR::fromList(gll)
      ugll = unique(unlist(gll))
      rownames(upGll) = ugll
      datatable(upGll, 
                filter = list(position = 'top', clear = FALSE))
    })
    output$debugTableSelected <- renderText({
      if(is.null(input$debugTable_rows_all)){
        return("noting")
      }
      gll <- gll()
      upGll = UpSetR::fromList(gll)
      ugll = unique(unlist(gll))
      rownames(upGll) = ugll
      return(paste(rownames(upGll)[input$debugTable_rows_all], collapse = " "))
    })
    
    output$vennlists <- renderPlot({
      shiny::validate(
        need(all(sapply(gll(),function(arg) !is.null(arg))),
             message = "Some lists are empty - make sure you extracted the results using the annotation object")
        
      )
      
      gplots::venn(gll())
    })
    
    output$upsetLists <- renderPlot({
      shiny::validate(
        need(sum(sapply(gll(),function(arg) length(arg)>0)) > 1,
             message = "Make sure you provide at least two sets")
      )
      UpSetR::upset(fromList(gll()))
    })
    
    
    output$printUPgenes <- renderPrint({
      print(head(values$genelistUP()))
      print(str(values$genelistUP()))
      organism <- annoSpecies_df[values$cur_species,]$species_short
      backgroundgenes <- rownames(values$dds_obj)[rowSums(counts(values$dds_obj))>0]
      inputType <- "SYMBOL" # will be replaced by input$...
      # annopkg <- paste0("org.",organism,".eg.db")
      annopkg <- annoSpecies_df[values$cur_species,]$pkg
      listGenesEntrez <- as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = values$genelistUP(),
                                                            column="ENTREZID", keytype=inputType))
      listBackgroundEntrez <- as.character(AnnotationDbi::mapIds(eval(parse(text=annopkg)), keys = backgroundgenes,
                                                                 column="ENTREZID", keytype=values$cur_type))
      
      # print(genelistUP())
      print(str(listGenesEntrez))
      print(class(listGenesEntrez))
      print(str(listBackgroundEntrez))
      print(class(listBackgroundEntrez))
      print(head(listGenesEntrez))
      print(head(listBackgroundEntrez))
      
    })
    
    output$download_plot_vennlists <- downloadHandler(filename = function() {
      input$filename_plot_vennlists
    }, content = function(file) {
      pdf(file)
      gplots::venn(gll())
      dev.off()
    })
    
    output$download_plot_upsetlists <- downloadHandler(filename = function() {
      input$filename_plot_upsetlists
    }, content = function(file) {
      pdf(file)
      UpSetR::upset(fromList(gll()))
      dev.off()
    })
    
    ## GO tbls topGO
    output$downloadGOTbl_up <- downloadHandler(
      filename = function() {
        "table_GOresults_up.csv"
      },
      content = function(file) {
        write.csv(values$topgo_up, file)
      }
    )
    output$downloadGOTbl_down <- downloadHandler(
      filename = function() {
        "table_GOresults_down.csv"
      },
      content = function(file) {
        write.csv(values$topgo_down, file)
      }
    )
    output$downloadGOTbl_updown <- downloadHandler(
      filename = function() {
        "table_GOresults_updown.csv"
      },
      content = function(file) {
        write.csv(values$topgo_updown, file)
      }
    )
    output$downloadGOTbl_l1 <- downloadHandler(
      filename = function() {
        "table_GOresults_list1.csv"
      },
      content = function(file) {
        write.csv(values$topgo_list1, file)
      }
    )
    output$downloadGOTbl_l2 <- downloadHandler(
      filename = function() {
        "table_GOresults_list2.csv"
      },
      content = function(file) {
        write.csv(values$topgo_list2, file)
      }
    )
    output$downloadGOTbl_l3 <- downloadHandler(
      filename = function() {
        "table_GOresults_list3.csv"
      },
      content = function(file) {
        write.csv(values$topgo_list3, file)
      }
    )
    output$downloadGOTbl_l4 <- downloadHandler(
      filename = function() {
        "table_GOresults_list4.csv"
      },
      content = function(file) {
        write.csv(values$topgo_list4, file)
      }
    )
    
    
  })
}