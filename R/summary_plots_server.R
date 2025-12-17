summary_plots_server <- function(id, values, annoSpecies_df, exportPlots) {
  # Define DEBUG variable - can be set to TRUE/FALSE to enable/disable debug output
  # Set options(wein.debug = TRUE) before launching the app to enable debug output
  DEBUG <- getOption("wein.debug", FALSE)
  
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    # Reactive value to store selected genes
    selected_genes <- reactiveVal(NULL)
    
    # Subsampled data reactive - used by all plots
    subsampled_res <- reactive({
      shiny::validate(
        need(!is.null(values$res_obj), "Results object is not available. Please generate results first."),
        need(!is.null(values$FDR), "FDR is not set.")
      )
      
      if (DEBUG) cat(file = stderr(), "=== SUBSAMPLING DATA ===\n")
      
      res_df <- as.data.frame(values$res_obj)
      res_df$ID <- rownames(res_df)
      
      # Always keep significant genes
      sig_idx <- which(!is.na(res_df$padj) & res_df$padj < values$FDR)
      nonsig_idx <- which(is.na(res_df$padj) | res_df$padj >= values$FDR)
      
      if (DEBUG) cat(file = stderr(), "Total genes:", nrow(res_df),
                     "Significant:", length(sig_idx),
                     "Non-significant:", length(nonsig_idx), "\n")
      
      # If too many genes, subsample non-significant
      if (nrow(res_df) > input$max_points) {
        keep_nonsig <- min(input$max_points - length(sig_idx), length(nonsig_idx))
        if (keep_nonsig > 0 && length(nonsig_idx) > keep_nonsig) {
          sampled_nonsig <- sample(nonsig_idx, keep_nonsig)
          keep_idx <- c(sig_idx, sampled_nonsig)
        } else {
          keep_idx <- c(sig_idx, nonsig_idx)
        }
        
        res_df <- res_df[keep_idx, ]
        if (DEBUG) cat(file = stderr(), "Subsampled to:", nrow(res_df), "genes\n")
      }
      # browser()
      # Convert back to DESeqResults-like object
      res_sub <- values$res_obj[res_df$ID, ]
      res_sub
    }) %>% bindCache(values$res_obj, values$FDR, input$max_points)
    
    # output explore_res ----
    # not used?
    output$explore_res <- renderPrint({
      if (DEBUG) cat(file = stderr(), "explore_res\n")
      
      expfac <- attributes(terms.formula(design(values$dds_obj)))$term.labels
      expfac # plus, support up to four factors that are either there or not according to the length
    })
    
    # output plotma ----
    output$plotma <- shiny::renderPlot({
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== PLOTMA START ===\n")
      # Validate that required data is available
      if (DEBUG) cat(file = stderr(), paste0("res_obj exists: ", !is.null(values$res_obj), "\n"))
      if (DEBUG) cat(file = stderr(), paste0("annotation_obj exists: ", !is.null(values$annotation_obj), "\n"))
      if (DEBUG) cat(file = stderr(), paste0("FDR value: ", values$FDR, "\n"))
      
      shiny::validate(
        need(!is.null(values$res_obj), "Results object is not available. Please generate results first."),
        need(!is.null(values$annotation_obj), "Annotation object is not available. Please set annotation first."),
        need(!is.null(values$FDR), "FDR is not set.")  # ← Add this
      )
      
      if (DEBUG) cat(file = stderr(), "in plotma\n")
      # browser()
      p <-     plot_ma(res_obj = subsampled_res(),annotation_obj = values$annotation_obj,FDR = values$FDR)
      
      exportPlots$plot_ma <- p
      # ggsave("test.pdf", p)
      # cat(file = stderr(), "MA plot saved as test.pdf\n")
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== PLOTMA END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      p
    })
    
    # Throttling for brush events to prevent excessive re-rendering
    # brush_throttle <- reactiveVal(0)
    
    # observeEvent(input$ma_brush, {
    #   # Throttle by 200ms
    #   invalidateLater(200, session)
    #   brush_throttle(isolate(brush_throttle()) + 1)
    # })
    ma_brush_debounced <- debounce({
      if (DEBUG) cat(file = stderr(), "ma_brush_debounced\n")
      reactive(input$ma_brush)
    }, 500)  # 500ms delay
    
    # 2. Add observer to UPDATE text input when brushing (append genes)
    observeEvent(ma_brush_debounced(), {
      if (!is.null(ma_brush_debounced())) {
        if (DEBUG) cat(file = stderr(), "in brushed\n")
        brushedObject <- curData()
        if (nrow(brushedObject) > 0) {
          selectedGenes <- as.character(brushedObject$ID)
          
          # Get gene names
          geneNames <- values$annotation_obj$gene_name[match(selectedGenes, 
                                                             rownames(values$annotation_obj))]
          displayNames <- ifelse(!is.na(geneNames) & geneNames != "", geneNames, selectedGenes)
          
          # Append to existing genes in text input
          current_text <- input$gene_list_input
          if (is.null(current_text) || current_text == "") {
            new_text <- paste(displayNames, collapse = " ")
          } else {
            existing_genes <- unlist(strsplit(current_text, "\\s+"))
            all_genes <- unique(c(existing_genes, displayNames))
            new_text <- paste(all_genes, collapse = " ")
          }
          
          # Update text input
          updateTextAreaInput(session, "gene_list_input", value = new_text)
          
          # Also update selected_genes reactive
          all_genes_list <- unlist(strsplit(new_text, "\\s+"))
          selected_genes(all_genes_list[all_genes_list != ""])
        }
        if (DEBUG) cat(file = stderr(), "leaving brushed\n")
        
      }
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
    
    
    
    # output mazoom ----
    output$mazoom <- renderPlot({
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== MAZOOM START ===\n")
      shiny::validate(
        need(!is.null(values$res_obj), "Results object is not available. Please generate results first."),
        need(!is.null(values$annotation_obj), "Annotation object is not available. Please set annotation first.")
      )
      # brush_throttle()
      if (DEBUG) cat(file = stderr(), "mazoom\n")
      
      # Suppress par() warnings that commonly occur in Shiny reactive contexts
      # p <- suppressWarnings({
      if (DEBUG) cat(file = stderr(), "in zoomed MA plot\n")
      manual_genes <- selected_genes()
      
      if (!is.null(manual_genes) && length(manual_genes) > 0) {
        # Create base plot
        p <- plot_ma(subsampled_res(), annotation_obj = values$annotation_obj,
                     FDR = values$FDR)
        # browser()
        mama <- data.frame(mean=subsampled_res()$baseMean,
                           lfc=subsampled_res()$log2FoldChange,
                           padj = subsampled_res()$padj,
                           isDE= ifelse(is.na(subsampled_res()$padj), FALSE,
                                        subsampled_res()$padj < 0.10),
                           ID=rownames(subsampled_res()))
        mama$genename <- values$annotation_obj$gene_name[match(mama$ID,
                                                               rownames(values$annotation_obj))]
        mama$logmean <- log10(mama$mean)
        
        matched_by_name <- which(mama$genename %in% manual_genes)
        matched_by_id <- which(mama$ID %in% manual_genes)
        matched_indices <- unique(c(matched_by_name, matched_by_id))
        
        if (length(matched_indices) > 0) {
          selected_data <- mama[matched_indices, ]
          selected_data <- selected_data[!is.na(selected_data$logmean) &
                                           !is.na(selected_data$lfc), ]
          
          # Add small points
          p <- p +
            geom_point(data = selected_data,
                       aes(x = logmean, y = lfc),
                       color = "steelblue",
                       size = 1.5,
                       alpha = 0.8)
          
          # Add labels with ggrepel
          p <- p +
            ggrepel::geom_text_repel(
              data = selected_data,
              aes(x = logmean, y = lfc, label = genename),
              color = "steelblue",
              size = 4,
              fontface = "bold",
              max.overlaps = 50,
              box.padding = 0.5,
              point.padding = 0.3,
              min.segment.length = 0
            )
          
          x_range <- range(selected_data$logmean, na.rm = TRUE)
          y_range <- range(selected_data$lfc, na.rm = TRUE)
          x_padding <- max(0.5, 0.15 * diff(x_range))
          y_padding <- max(0.5, 0.15 * diff(y_range))
          
          p <- p + coord_cartesian(
            xlim = c(x_range[1] - x_padding, x_range[2] + x_padding),
            ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)
          )
        }
      } else {
        if(is.null(ma_brush_debounced())) {
          return(ggplot() +
                   annotate("text", label="click and drag to zoom in", 0, 0) +
                   theme_bw())
        }
        
        p <- plot_ma(subsampled_res(), annotation_obj = values$annotation_obj,
                     FDR = values$FDR) +
          coord_cartesian(xlim = c(ma_brush_debounced()$xmin, ma_brush_debounced()$xmax),
                          ylim = c(ma_brush_debounced()$ymin, ma_brush_debounced()$ymax))
      }
      exportPlots$plot_mazoom <- p
      # })
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== MAZOOM END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      p
    })
    
    
    
    
    # Function to validate gene names against the DESeq object
    validate_genes <- function(gene_list, dds_obj, annotation_obj) {
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== VALIDATE_GENES START ===\n")
      
      if (is.null(dds_obj) || is.null(annotation_obj)) {
        end_time <- Sys.time()
        if (DEBUG) cat(file = stderr(), paste0("=== VALIDATE_GENES END: ",
                                               round(as.numeric(end_time - start_time, units = "secs"), 2),
                                               " seconds ===\n"))
        return(list(valid = character(0), invalid = gene_list))
      }
      
      # Get all available gene names from annotation
      all_gene_names <- annotation_obj$gene_name
      all_gene_ids <- rownames(annotation_obj)
      
      # Check if genes are in the annotation object
      valid_genes <- gene_list[gene_list %in% c(all_gene_names, all_gene_ids)]
      invalid_genes <- gene_list[!gene_list %in% c(all_gene_names, all_gene_ids)]
      
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== VALIDATE_GENES END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      return(list(valid = valid_genes, invalid = invalid_genes))
    }
    
    # Observe add genes button
    observeEvent(input$add_genes_button, {
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== ADD_GENES_BUTTON START ===\n")
      # Parse gene list from input (split by whitespace)
      gene_input <- input$gene_list_input
      if (is.null(gene_input) || gene_input == "") {
        showNotification("Please enter gene names to add.", type = "warning")
        return()
      }
      
      # Split by whitespace
      new_genes <- unlist(strsplit(gene_input, "\\s+"))
      new_genes <- new_genes[new_genes != ""]
      
      if (length(new_genes) == 0) {
        showNotification("No valid gene names found in input.", type = "warning")
        return()
      }
      
      # Validate genes
      validation_result <- validate_genes(new_genes, values$dds_obj, values$annotation_obj)
      
      # Notify user about invalid genes
      if (length(validation_result$invalid) > 0) {
        showNotification(paste("Invalid gene names:", paste(validation_result$invalid, collapse = ", ")),
                         type = "error", duration = 10)
      }
      
      # Add valid genes to selection
      if (length(validation_result$valid) > 0) {
        current_genes <- selected_genes()
        if (is.null(current_genes)) {
          selected_genes(validation_result$valid)
          # Update the text area to show all selected genes
          updateTextAreaInput(session, ns("gene_list_input"), value = paste(validation_result$valid, collapse = " "))
        } else {
          # Combine with existing genes, removing duplicates
          combined_genes <- unique(c(current_genes, validation_result$valid))
          selected_genes(combined_genes)
          # Update the text area to show all selected genes
          updateTextAreaInput(session, ns("gene_list_input"), value = paste(combined_genes, collapse = " "))
        }
        showNotification(paste("Added", length(validation_result$valid), "genes to selection."), type = "message")
      } else if (length(validation_result$invalid) == 0) {
        showNotification("No valid genes found in input.", type = "warning")
      }
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== ADD_GENES_BUTTON END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
    })
    
    
    
    curData <- reactive({
      shiny::validate(
        need(!is.null(values$res_obj), "Results object is not available. Please generate results first."),
        need(!is.null(values$annotation_obj), "Annotation object is not available. Please set annotation first.")
      )
      
      if (DEBUG) cat(file = stderr(), "=== CURDATA START ===\n")
      start_time <- Sys.time()
      
      # Create base data frame
      mama <- data.frame(
        mean = subsampled_res()$baseMean,
        lfc = subsampled_res()$log2FoldChange,
        padj = subsampled_res()$padj,
        isDE = ifelse(is.na(subsampled_res()$padj), FALSE, subsampled_res()$padj < 0.10),
        ID = rownames(subsampled_res())
      )
      mama$genename <- values$annotation_obj$gene_name[match(mama$ID, rownames(values$annotation_obj))]
      mama$yesorno <- ifelse(mama$isDE, "red", "black")
      mama$logmean <- log10(mama$mean)
      
      # Get brushed genes
      brushed <- brushedPoints(mama, ma_brush_debounced(), xvar = "logmean", yvar = "lfc")
      brushed_ids <- brushed$ID
      
      # Get manually selected genes
      manual_genes <- selected_genes()
      manual_ids <- character(0)
      
      if (!is.null(manual_genes) && length(manual_genes) > 0) {
        matched_by_name <- which(mama$genename %in% manual_genes)
        matched_by_id <- which(mama$ID %in% manual_genes)
        matched_indices <- unique(c(matched_by_name, matched_by_id))
        manual_ids <- mama$ID[matched_indices]
      }
      
      # COMBINE both brushed and manual
      all_ids <- unique(c(brushed_ids, manual_ids))
      
      # LIMIT TO 100 with warning
      if (length(all_ids) > 100) {
        showNotification(
          paste0("Too many genes selected (", length(all_ids), 
                 "). Showing only the first 100 genes."),
          type = "warning",
          duration = 5
        )
        all_ids <- all_ids[1:100]
      }
      
      if (length(all_ids) > 0) {
        res <- mama[mama$ID %in% all_ids, ]
      } else {
        res <- mama[numeric(0), ]
      }
      
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== CURDATA END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      
      res
    })
    # Observe clear genes button
    observeEvent(input$clear_genes_button, {
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== CLEAR_GENES_BUTTON START ===\n")
      
      selected_genes(NULL)
      updateTextAreaInput(session, "gene_list_input", value = "")
      showNotification("Gene selection cleared.", type = "message")
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== CLEAR_GENES_BUTTON END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
    })
    
    
    # Throttling for click events to prevent excessive re-rendering
    # click_throttle <- reactiveVal(0)
    
    # observeEvent(input$mazoom_click, {
    #   # Throttle by 300ms
    #   invalidateLater(300, session)
    #   click_throttle(isolate(click_throttle()) + 1)
    # })
    
    observeEvent( input$mazoom_click, {
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== observeEvent click START ===\n")
      if (DEBUG) cat(file = stderr(), "Click coords:",
                     input$mazoom_click$x, input$mazoom_click$y, input$mazoom_click$ID, "\n")
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== observeEvent click END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      # browser()
    })
    
    clickOut <- reactive({
      if (DEBUG) cat(file = stderr(), "Click clickOut:\n")
      input$mazoom_click
    })
    # curDataClick reactive ----
    curDataClick <- reactive({
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== CURDATACLICK START ===\n")
      
      # Use throttle to control frequency
      shiny::validate(
        need(!is.null(values$res_obj), "Results object is not available. Please generate results first."),
        need(!is.null(values$annotation_obj), "Annotation object is not available. Please set annotation first."),
        need(!is.null(clickOut()), "Please click on a gene in the zoomed MA plot to select it"),
        need(!is.null(input$ma_brush), "Please brush an area in the MA plot to select genes")
      )
      
      mama <- data.frame(
        mean=subsampled_res()$baseMean,
        lfc=subsampled_res()$log2FoldChange,
        padj = subsampled_res()$padj,
        isDE= ifelse(is.na(subsampled_res()$padj), FALSE, subsampled_res()$padj < 0.10),
        ID=rownames(subsampled_res())
      )
      mama$genename <- values$annotation_obj$gene_name[match(mama$ID,rownames(values$annotation_obj))]
      # mama$yesorno <- ifelse(mama$isDE,"yes","no")
      mama$yesorno <- ifelse(mama$isDE,"red","black")
      mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
      
      res <- nearPoints(mama,
                        clickOut(),
                        threshold = 20,
                        maxpoints = 1,
                        addDist = T)
      
      if (DEBUG) cat(file = stderr(), "nearPoints returned", nrow(res), "rows\n")
      if (nrow(res) > 0) {
        if (DEBUG) cat(file = stderr(), "Found gene:", res$ID[1], "\n")
        if (DEBUG) cat(file = stderr(), "Gene name:", res$genename[1], "\n")
      } else {
        if (DEBUG) cat(file = stderr(), "No points found near click!\n")
      }
      
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== CURDATACLICK END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      res
    })
    
    
    
    # output ma_brush_out ----
    output$ma_brush_out <- DT::renderDataTable({
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== MA_BRUSH_OUT START ===\n")
      
      brushedObject <- curData()
      if(nrow(brushedObject)==0)
        return(NULL)
      # browser()
      # curData()
      if (DEBUG) cat(file = stderr(), "ma_brush_out\n")
      selectedGenes <- as.character(brushedObject$ID)
      rownames(brushedObject) = selectedGenes
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== MA_BRUSH_OUT END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      brushedObject
    },options=list(pageLength=100))
    
    # output selectedGenesMAplot ----
    output$selectedGenesMAplot <- renderText({
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== SELECTEDGENESMAPLOT START ===\n")
      
      brushedObject <- curData()
      if (DEBUG) cat(file = stderr(), "selectedGenesMAplot\n")
      selectedGenes <- as.character(brushedObject$ID)
      
      # If we have an annotation object, try to use gene names instead of IDs
      if (!is.null(values$annotation_obj) && nrow(brushedObject) > 0) {
        geneNames <- values$annotation_obj$gene_name[match(selectedGenes, rownames(values$annotation_obj))]
        # Use gene names where available, fallback to IDs
        displayNames <- ifelse(!is.na(geneNames) & geneNames != "", geneNames, selectedGenes)
        end_time <- Sys.time()
        if (DEBUG) cat(file = stderr(), paste0("=== SELECTEDGENESMAPLOT END: ",
                                               round(as.numeric(end_time - start_time, units = "secs"), 2),
                                               " seconds ===\n"))
        paste(displayNames, collapse = " ")
      } else {
        end_time <- Sys.time()
        if (DEBUG) cat(file = stderr(), paste0("=== SELECTEDGENESMAPLOT END: ",
                                               round(as.numeric(end_time - start_time, units = "secs"), 2),
                                               " seconds ===\n"))
        paste(selectedGenes, collapse = " ")
      }
    })
    # output heatbrush ----
    output$heatbrush <- renderPlotly({
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== HEATBRUSH START ===\n")
      
      # Check if we have manually selected genes
      manual_genes <- selected_genes()
      if (is.null(values$dds_obj)) {
        end_time <- Sys.time()
        if (DEBUG) cat(file = stderr(), paste0("=== HEATBRUSH END: ",
                                               round(as.numeric(end_time - start_time, units = "secs"), 2),
                                               " seconds ===\n"))
        return(NULL)
      }
      
      if (DEBUG) cat(file = stderr(), "heatbrush\n")
      brushedObject <- curData()
      selectedGenes <- as.character(brushedObject$ID)
      
      if (length(selectedGenes) == 0) {
        end_time <- Sys.time()
        if (DEBUG) cat(file = stderr(), paste0("=== HEATBRUSH END: ",
                                               round(as.numeric(end_time - start_time, units = "secs"), 2),
                                               " seconds ===\n"))
        return(NULL)
      }
      
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
      
      p <- heatmaply::heatmaply(toplot,cluster_cols = as.logical(input$heatmap_colv))
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== HEATBRUSH END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      p
    })
    
    
    
    # output hpi_brush ----
    output$hpi_brush <- renderPlotly({
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== HPI_BRUSH START ===\n")
      
      # Check if we have manually selected genes
      manual_genes <- selected_genes()
      if (is.null(values$dds_obj)) {
        end_time <- Sys.time()
        if (DEBUG) cat(file = stderr(), paste0("=== HPI_BRUSH END: ",
                                               round(as.numeric(end_time - start_time, units = "secs"), 2),
                                               " seconds ===\n"))
        return(NULL)
      }
      
      if (DEBUG) cat(file = stderr(), "hpi_brush\n")
      brushedObject <- curData()
      selectedGenes <- as.character(brushedObject$ID)
      
      if (length(selectedGenes) == 0) {
        end_time <- Sys.time()
        if (DEBUG) cat(file = stderr(), paste0("=== HPI_BRUSH END: ",
                                               round(as.numeric(end_time - start_time, units = "secs"), 2),
                                               " seconds ===\n"))
        return(NULL)
      }
      
      toplot <- assay(values$dds_obj)[selectedGenes,]
      rownames(toplot) <- values$annotation_obj$gene_name[match(rownames(toplot),rownames(values$annotation_obj))]
      mycolss <- c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#fee090","#fdae61","#f46d43","#d73027","#a50026") # to be consistent with red/blue usual coding
      if(input$pseudocounts) toplot <- log2(1+toplot)
      if(input$rowscale) toplot <- mat_rowscale(toplot)
      if(nrow(toplot) <1){
        if (DEBUG) cat(file = stderr(), "\ntoplot nrow <1\n")
        end_time <- Sys.time()
        if (DEBUG) cat(file = stderr(), paste0("=== HPI_BRUSH END: ",
                                               round(as.numeric(end_time - start_time, units = "secs"), 2),
                                               " seconds ===\n"))
        return(NULL)
      }
      # Add early stopping for very large gene sets to prevent performance issues
      if(nrow(toplot) > 1000) {
        showNotification("Large gene set detected. Displaying first 1000 genes for performance.",
                         type = "warning", duration = 10)
        toplot <- toplot[1:min(1000, nrow(toplot)), , drop = FALSE]
      }
      
      p <- heatmaply::heatmaply(toplot,Colv = as.logical(input$heatmap_colv),colors = mycolss, cexCol = 1)
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== HPI_BRUSH END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      p
    })
    
    # output deb ----
    output$deb <- renderPrint({
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== DEB START ===\n")
      
      curDataClick()
      if (DEBUG) cat(file = stderr(), "deb\n")
      selectedGene <- curDataClick()$ID
      #         selectedGeneSymbol <- cm2$fromgtf[match(selectedGene,rownames(cm2))]
      #         # plotCounts(dds_cleaner,)
      #         genedata <- plotCounts(dds_cleaner,gene=selectedGene,intgroup = "condition",returnData = T)
      #         genedata
      # str(as.character(selectedGene))
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== DEB END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      selectedGene
    })
    
    # output volcanoplot ----
    output$volcanoplot <- renderPlot({
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== VOLCANOPLOT START ===\n")
      
      manual_genes <- selected_genes()
      
      # Only plot a subset of non-significant points
      res_for_plot <- subsampled_res()
      p = NULL
      # Suppress par() warnings that commonly occur in Shiny reactive contexts
      p <- suppressWarnings({
        if (DEBUG) cat(file = stderr(), "volcanoplot\n")
        if (!is.null(manual_genes) && length(manual_genes) > 0) {
          # Create base plot
          p <- plot_volcano(res_for_plot, FDR = values$FDR,
                            annotation_obj = values$annotation_obj)
          
          # Get data for selected genes
          res_df <- as.data.frame(res_for_plot)
          res_df$ID <- rownames(res_df)
          res_df$genename <- values$annotation_obj$gene_name[match(res_df$ID,
                                                                   rownames(values$annotation_obj))]
          
          matched_by_name <- which(res_df$genename %in% manual_genes)
          matched_by_id <- which(res_df$ID %in% manual_genes)
          matched_indices <- unique(c(matched_by_name, matched_by_id))
          
          if (length(matched_indices) > 0) {
            selected_data <- res_df[matched_indices, ]
            selected_data <- selected_data[!is.na(selected_data$log2FoldChange) &
                                             !is.na(selected_data$padj), ]
            
            # Add small points for selected genes
            p <- p +
              geom_point(data = selected_data,
                         aes(x = log2FoldChange, y = -log10(padj)),
                         color = "steelblue",
                         size = 2,
                         alpha = 0.8)
            
            # Add labels if checkbox is checked
            if (input$show_gene_names) {
              p <- p +
                ggrepel::geom_text_repel(
                  data = selected_data,
                  aes(x = log2FoldChange, y = -log10(padj), label = genename),
                  color = "steelblue",
                  size = 4,
                  fontface = "bold",
                  max.overlaps = 50,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  min.segment.length = 0
                )
            }
          }
        } else {
          p <- plot_volcano(res_for_plot, FDR = values$FDR,
                            annotation_obj = values$annotation_obj)
        }
        
        exportPlots$plot_volcanoplot <- p
        end_time <- Sys.time()
        if (DEBUG) cat(file = stderr(), paste0("=== VOLCANOPLOT END: ",
                                               round(as.numeric(end_time - start_time, units = "secs"), 2),
                                               " seconds ===\n"))
        # ggsave("test.pdf", p)
        # cat(file = stderr(), "MA plot saved as test.pdf\n")
        p
      })
      p
    })
    
    # server genefinder --------------------------------------------------------
    # output genefinder_plot ----
    output$genefinder_plot <- renderPlot({
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== GENEFINDER_PLOT START ===\n")
      
      shiny::validate(
        need(
          !is.null(values$annotation_obj),
          "Please set the annotation"
        )
      )
      if (DEBUG) cat(file = stderr(), "genefinder_plot\n")
      
      shiny::validate(
        need(!is.null(input$mazoom_click),
             "Please click on a gene in the zoomed MA plot to select it")
      )
      
      # TRY to get clicked data
      clickData <- tryCatch({
        curDataClick()
      }, error = function(e) {
        return(NULL)
      })
      
      # Validate we got data
      shiny::validate(
        need(!is.null(clickData) && nrow(clickData) > 0,
             "No gene found at click location. Try clicking closer to a point.")
      )        
      selectedGene <- as.character(curDataClick()$ID)
      
      
      selectedGeneSymbol <- values$annotation_obj$gene_name[match(selectedGene,values$annotation_obj$gene_id)]
      
      # Suppress par() warnings that commonly occur in Shiny reactive contexts
      p <- suppressWarnings({
        ggplotCounts(values$dds_obj, selectedGene, intgroup = values$color_by, annotation_obj=values$annotation_obj)
      })
      
      if(input$ylimZero_genes)
        p <- p + coord_cartesian(ylim = c(0.1, NA))
      
      exportPlots$plot_genefinder <- p
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== GENEFINDER_PLOT END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      p
    })
    
    # output rentrez_infobox ----
    output$rentrez_infobox <- renderUI({
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== RENTREZ_INFOBOX START ===\n")
      
      # Check if we have manually selected genes
      manual_genes <- selected_genes()
      # Check if we have a click
      shiny::validate(
        need(!is.null(input$mazoom_click),
             "Select a gene first to display additional info")
      )
      
      # TRY to get clicked data
      clickData <- tryCatch({
        curDataClick()
      }, error = function(e) {
        return(NULL)
      })
      
      shiny::validate(
        need(!is.null(clickData) && nrow(clickData) > 0,
             "Select a gene first to display additional info (retrieved from the NCBI/ENTREZ db website)")
      )
      
      selectedGene <- as.character(curDataClick()$ID)
      
      
      if (DEBUG) cat(file = stderr(), "rentrez_infobox\n")
      shiny::validate(
        need(
          (!is.null(values$cur_species)),
          "Select a species first in the Data Setup panel"
        )
      )
      
      # Handle the case where cur_type might be SYMBOL or other invalid keytype
      annopkg <- get(annoSpecies_df[values$cur_species,]$pkg)
      selgene_entrez <- NULL
      
      # If cur_type is SYMBOL, we need to convert it to a valid keytype
      if (values$cur_type == "SYMBOL") {
        # When cur_type is SYMBOL, the selectedGene is already a gene symbol
        # We need to map from SYMBOL to ENTREZID using the annotation object
        # First, try to find the ENTREZID in the annotation object
        gene_idx <- match(selectedGene, values$annotation_obj$gene_name)
        if (!is.na(gene_idx)) {
          selgene_entrez <- rownames(values$annotation_obj)[gene_idx]
        }
        
        # If that doesn't work, try to use mapIds with "SYMBOL" as column and keytype as rownames
        if (is.null(selgene_entrez) || selgene_entrez == selectedGene) {
          tryCatch({
            selgene_entrez <- mapIds(annopkg, selectedGene, "ENTREZID", "SYMBOL", multiVals = "first")
          }, error = function(e) {
            # If that fails, try with the rownames as keytype
            tryCatch({
              selgene_entrez <- mapIds(annopkg, selectedGene, "ENTREZID", keytype = keytypes(annopkg)[1], multiVals = "first")
            }, error = function(e2) {
              # If all fails, use the selectedGene directly (might be an ENTREZID already)
              selgene_entrez <- selectedGene
            })
          })
        }
      } else {
        # For other keytypes, use mapIds normally
        tryCatch({
          selgene_entrez <- mapIds(annopkg, selectedGene, "ENTREZID", values$cur_type, multiVals = "first")
        }, error = function(e) {
          # If mapIds fails, try with the first available keytype
          tryCatch({
            selgene_entrez <- mapIds(annopkg, selectedGene, "ENTREZID", keytype = keytypes(annopkg)[1], multiVals = "first")
          }, error = function(e2) {
            # If all fails, use the selectedGene directly
            selgene_entrez <- selectedGene
          })
        })
      }
      
      # Get gene information from entrez
      fullinfo <- NULL
      if (!is.null(selgene_entrez) && selgene_entrez != "") {
        tryCatch({
          fullinfo <- geneinfo(selgene_entrez)
        }, error = function(e) {
          # If geneinfo fails, create a minimal info object
          fullinfo <- list(
            name = if (!is.null(selgene_entrez)) selgene_entrez else selectedGene,
            description = "Gene information not available",
            summary = ""
          )
        })
      } else {
        # Create minimal info if no entrez id
        fullinfo <- list(
          name = selectedGene,
          description = "Gene information not available",
          summary = ""
        )
      }
      
      # Build up link manually to paste under the info
      link_pubmed <- paste0('<a href="http://www.ncbi.nlm.nih.gov/gene/?term=',
                            if (!is.null(selgene_entrez)) selgene_entrez else selectedGene,
                            '" target="_blank" >Click here to see more at NCBI</a>')
      
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== RENTREZ_INFOBOX END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      if (is.null(fullinfo) || fullinfo$summary == "")
        return(HTML(paste0("<b>", if (!is.null(fullinfo$name)) fullinfo$name else selectedGene, "</b><br/><br/>",
                           if (!is.null(fullinfo$description)) fullinfo$description else "No description available", "<br/><br/>",
                           link_pubmed
        )))
      else
        return(HTML(paste0("<b>", if (!is.null(fullinfo$name)) fullinfo$name else selectedGene, "</b><br/><br/>",
                           if (!is.null(fullinfo$description)) fullinfo$description else "No description available", "<br/><br/>",
                           fullinfo$summary, "<br/><br/>",
                           link_pubmed
        )))
    })
    
    
    # output table_combi ----
    cur_combires <- reactive({
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== CUR_COMBIRES START ===\n")
      shiny::validate(
        need(!is.null(values$res_obj), "Results object is not available. Please generate results first."),
        need(!is.null(values$dds_obj), "DESeq object is not available. Please run DESeq analysis first."),
        need(!is.null(values$annotation_obj), "Annotation object is not available. Please set annotation first.")
      )
      if (DEBUG) cat(file = stderr(), "cur_combires\n")
      
      normCounts <- as.data.frame(counts(estimateSizeFactors(values$dds_obj), normalized=TRUE))
      normCounts$id <- rownames(normCounts)
      res_df <- deseqresult2tbl(subsampled_res())
      
      combi_obj <- dplyr::inner_join(res_df, normCounts, by="id")
      combi_obj$symbol <- values$annotation_obj$gene_name[match(combi_obj$id, values$annotation_obj$gene_id)]
      
      if("symbol" %in% names(subsampled_res())) {
        sel_genes <- values$avail_symbols
        sel_genes_ids <- values$annotation_obj$gene_id[match(sel_genes, values$annotation_obj$gene_name)]
      } else {
        sel_genes_ids <- values$avail_ids
      }
      
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== CUR_COMBIRES END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      if(length(sel_genes_ids) > 0) {
        combi_obj[match(sel_genes_ids, combi_obj$id), ]
      } else {
        combi_obj
      }
    })
    
    # Then define output only once:
    output$table_combi <- DT::renderDataTable({
      start_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), "=== TABLE_COMBI START ===\n")
      
      if (DEBUG) cat(file = stderr(), "values$res_obj changed\n")
      dt <- datatable(cur_combires(), options = list(scrollX=TRUE))
      end_time <- Sys.time()
      if (DEBUG) cat(file = stderr(), paste0("=== TABLE_COMBI END: ",
                                             round(as.numeric(end_time - start_time, units = "secs"), 2),
                                             " seconds ===\n"))
      dt
    })
    
    
    # output download_plot_ma ----
    output$download_plot_ma <- downloadHandler(filename = function() {
      input$filename_plot_ma
    }, content = function(file) {
      ggsave(file, exportPlots$plot_ma, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    # output download_plot_mazoom ----
    output$download_plot_mazoom <- downloadHandler(filename = function() {
      input$filename_plot_mazoom
    }, content = function(file) {
      ggsave(file, exportPlots$plot_mazoom, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    
    
    # output download_plot_volcanoplot ----
    output$download_plot_volcanoplot <- downloadHandler(filename = function() {
      input$filename_plot_volcanoplot
    }, content = function(file) {
      ggsave(file, exportPlots$plot_volcanoplot, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    # output download_plot_genefinder ----
    output$download_plot_genefinder <- downloadHandler(filename = function() {
      input$filename_plot_genefinder
    }, content = function(file) {
      ggsave(file, exportPlots$plot_genefinder, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    # base graphics plots
    # output download_plot_heatbrush ----
    output$download_plot_heatbrush <- downloadHandler(filename = function() {
      input$filename_plot_heatbrush
    }, content = function(file) {
      pdf(file)
      # Check if we have manually selected genes
      manual_genes <- selected_genes()
      if (!is.null(manual_genes) && length(manual_genes) > 0) {
        # Create data frame with all manually selected genes
        mama <- data.frame(mean=subsampled_res()$baseMean,lfc=subsampled_res()$log2FoldChange,padj = subsampled_res()$padj,isDE= ifelse(is.na(subsampled_res()$padj), FALSE, subsampled_res()$padj < 0.10),ID=rownames(subsampled_res()))
        mama$genename <- values$annotation_obj$gene_name[match(mama$ID,rownames(values$annotation_obj))]
        
        # Filter to only include manually selected genes
        # First try to match by gene name, then by ID
        matched_by_name <- which(mama$genename %in% manual_genes)
        matched_by_id <- which(mama$ID %in% manual_genes)
        matched_indices <- unique(c(matched_by_name, matched_by_id))
        
        if (length(matched_indices) > 0) {
          selectedGenes <- mama$ID[matched_indices]
        } else {
          selectedGenes <- character(0)
        }
      } else {
        brushedObject <- curData()
        selectedGenes <- as.character(brushedObject$ID)
      }
      
      if (length(selectedGenes) > 0) {
        toplot <- assay(values$dds_obj)[selectedGenes,]
        rownames(toplot) <- values$annotation_obj$gene_name[match(rownames(toplot),rownames(values$annotation_obj))]
        
        if(input$pseudocounts) toplot <- log2(1+toplot)
        
        if(input$rowscale) toplot <- mat_rowscale(toplot)
        
        heatmaply(toplot,cluster_cols = as.logical(input$heatmap_colv))
      }
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