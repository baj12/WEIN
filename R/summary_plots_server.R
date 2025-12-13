summary_plots_server <- function(id, values, annoSpecies_df, exportPlots) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    # Reactive value to store selected genes
    selected_genes <- reactiveVal(NULL)
    
    
    # not used?
    output$explore_res <- renderPrint({
      expfac <- attributes(terms.formula(design(values$dds_obj)))$term.labels
      expfac # plus, support up to four factors that are either there or not according to the length
    })
    
    output$plotma <- renderPlot({
          # Validate that required data is available
          shiny::validate(
            need(!is.null(values$res_obj), "Results object is not available. Please generate results first."),
            need(!is.null(values$annotation_obj), "Annotation object is not available. Please set annotation first.")
          )
          
      cat(file = stderr(), "in plotma\n")
      # browser()
          p <- suppressWarnings({
            plot_ma(values$res_obj,annotation_obj = values$annotation_obj,FDR = values$FDR)
          })
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
    
    # 2. Add observer to UPDATE text input when brushing (append genes)
    observeEvent(input$ma_brush, {
      if (!is.null(input$ma_brush)) {
      cat(file = stderr(), "in brushed\n")
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
      }
    }, ignoreInit = TRUE)
    
    
    
    output$mazoom <- renderPlot({
          brush_throttle()
          
          # Suppress par() warnings that commonly occur in Shiny reactive contexts
          p <- suppressWarnings({
            manual_genes <- selected_genes()
            
            if (!is.null(manual_genes) && length(manual_genes) > 0) {
              # Create base plot
              p <- plot_ma(values$res_obj, annotation_obj = values$annotation_obj,
                           FDR = values$FDR)
              
              mama <- data.frame(mean=values$res_obj$baseMean,
                                 lfc=values$res_obj$log2FoldChange,
                                 padj = values$res_obj$padj,
                                 isDE= ifelse(is.na(values$res_obj$padj), FALSE,
                                              values$res_obj$padj < 0.10),
                                 ID=rownames(values$res_obj))
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
              if(is.null(input$ma_brush)) {
                return(ggplot() +
                         annotate("text", label="click and drag to zoom in", 0, 0) +
                         theme_bw())
              }
              
              p <- plot_ma(values$res_obj, annotation_obj = values$annotation_obj,
                           FDR = values$FDR) +
                coord_cartesian(xlim = c(input$ma_brush$xmin, input$ma_brush$xmax),
                                ylim = c(input$ma_brush$ymin, input$ma_brush$ymax))
            }
            exportPlots$plot_mazoom <- p
            p
          })
        })
    
    
    
   
    # Function to validate gene names against the DESeq object
    validate_genes <- function(gene_list, dds_obj, annotation_obj) {
      if (is.null(dds_obj) || is.null(annotation_obj)) {
        return(list(valid = character(0), invalid = gene_list))
      }
      
      # Get all available gene names from annotation
      all_gene_names <- annotation_obj$gene_name
      all_gene_ids <- rownames(annotation_obj)
      
      # Check if genes are in the annotation object
      valid_genes <- gene_list[gene_list %in% c(all_gene_names, all_gene_ids)]
      invalid_genes <- gene_list[!gene_list %in% c(all_gene_names, all_gene_ids)]
      
      return(list(valid = valid_genes, invalid = invalid_genes))
    }
    
    # Observe add genes button
    observeEvent(input$add_genes_button, {
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
    })
    
    curData <- reactive({
      # If we have manually selected genes, use those
      manual_genes <- selected_genes()
      if (!is.null(manual_genes) && length(manual_genes) > 0) {
        # Create data frame with all manually selected genes
        mama <- data.frame(mean=values$res_obj$baseMean,lfc=values$res_obj$log2FoldChange,padj = values$res_obj$padj,isDE= ifelse(is.na(values$res_obj$padj), FALSE, values$res_obj$padj < 0.10),ID=rownames(values$res_obj))
        mama$genename <- values$annotation_obj$gene_name[match(mama$ID,rownames(values$annotation_obj))]
        
        # Filter to only include manually selected genes
        # First try to match by gene name, then by ID
        matched_by_name <- which(mama$genename %in% manual_genes)
        matched_by_id <- which(mama$ID %in% manual_genes)
        matched_indices <- unique(c(matched_by_name, matched_by_id))
        
        if (length(matched_indices) > 0) {
          res <- mama[matched_indices, ]
        } else {
          res <- mama[numeric(0), ]  # Empty data frame with same structure
        }
        return(res)
      }
      
      # Otherwise, use brushed data as before
      mama <- data.frame(mean=values$res_obj$baseMean,lfc=values$res_obj$log2FoldChange,padj = values$res_obj$padj,isDE= ifelse(is.na(values$res_obj$padj), FALSE, values$res_obj$padj < 0.10),ID=rownames(values$res_obj))
      mama$genename <- values$annotation_obj$gene_name[match(mama$ID,rownames(values$annotation_obj))]
      # mama$yesorno <- ifelse(mama$isDE,"yes","no")
      mama$yesorno <- ifelse(mama$isDE,"red","black")
      mama$logmean <- log10(mama$mean) # TO ALLOW FOR BRUSHING!!
      res <- brushedPoints(mama, input$ma_brush,xvar="logmean",yvar="lfc")
      res
    })
    
    
    # Observe clear genes button
    observeEvent(input$clear_genes_button, {
      selected_genes(NULL)
      updateTextAreaInput(session, "gene_list_input", value = "")
      showNotification("Gene selection cleared.", type = "message")
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
        brushedObject <- curData()
        if(nrow(brushedObject)==0)
          return(NULL)
        # browser()
        # curData()
        selectedGenes <- as.character(brushedObject$ID)
        rownames(brushedObject) = selectedGenes
        brushedObject
      },options=list(pageLength=100))
      
      output$selectedGenesMAplot <- renderText({
        brushedObject <- curData()
        selectedGenes <- as.character(brushedObject$ID)
        
        # If we have an annotation object, try to use gene names instead of IDs
        if (!is.null(values$annotation_obj) && nrow(brushedObject) > 0) {
          geneNames <- values$annotation_obj$gene_name[match(selectedGenes, rownames(values$annotation_obj))]
          # Use gene names where available, fallback to IDs
          displayNames <- ifelse(!is.na(geneNames) & geneNames != "", geneNames, selectedGenes)
          paste(displayNames, collapse = " ")
        } else {
          paste(selectedGenes, collapse = " ")
        }
      })
      output$heatbrush <- renderPlotly({
        # Check if we have manually selected genes
        manual_genes <- selected_genes()
        if (is.null(values$dds_obj)) {
          return(NULL)
        }
        
        brushedObject <- curData()
        selectedGenes <- as.character(brushedObject$ID)
        
        if (length(selectedGenes) == 0) {
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
        
        heatmaply::heatmaply(toplot,cluster_cols = as.logical(input$heatmap_colv))
      })
      
      
      
      output$hpi_brush <- renderPlotly({
        # Check if we have manually selected genes
        manual_genes <- selected_genes()
        if (is.null(values$dds_obj)) {
          return(NULL)
        }
        
        brushedObject <- curData()
        selectedGenes <- as.character(brushedObject$ID)
        
        if (length(selectedGenes) == 0) {
          return(NULL)
        }
        
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
      
      output$volcanoplot <- renderPlot({  # Changed from renderPlotly
              manual_genes <- selected_genes()
              
              # Suppress par() warnings that commonly occur in Shiny reactive contexts
              p <- suppressWarnings({
                if (!is.null(manual_genes) && length(manual_genes) > 0) {
                  # Create base plot
                  p <- plot_volcano(values$res_obj, FDR = values$FDR,
                                    annotation_obj = values$annotation_obj)
                  
                  # Get data for selected genes
                  res_df <- as.data.frame(values$res_obj)
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
                  p <- plot_volcano(values$res_obj, FDR = values$FDR,
                                    annotation_obj = values$annotation_obj)
                }
                
                exportPlots$plot_volcanoplot <- p
                p
              })
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
              
              # Allow gene finder to work with manually selected genes or brushed/clicked genes
              manual_genes <- selected_genes()
              if (!is.null(manual_genes) && length(manual_genes) > 0) {
                # Use the first manually selected gene if available
                selectedGene <- manual_genes[1]
                # Try to match by gene name first, then by ID
                matched_by_name <- which(values$annotation_obj$gene_name %in% selectedGene)
                if (length(matched_by_name) > 0) {
                  selectedGene <- rownames(values$annotation_obj)[matched_by_name[1]]
                }
              } else {
                shiny::validate(
                  need(!is.null(input$ma_brush),
                       "Please select a region on the MA plot by clicking and dragging")
                )
                shiny::validate(
                  need(!is.null(input$mazoom_click),
                       "Please click on a gene in the zoomed MA plot to select it")
                )
                
                selectedGene <- as.character(curDataClick()$ID)
              }
              
              selectedGeneSymbol <- values$annotation_obj$gene_name[match(selectedGene,values$annotation_obj$gene_id)]
              
              # Suppress par() warnings that commonly occur in Shiny reactive contexts
              p <- suppressWarnings({
                ggplotCounts(values$dds_obj, selectedGene, intgroup = values$color_by, annotation_obj=values$annotation_obj)
              })
              
              if(input$ylimZero_genes)
                p <- p + coord_cartesian(ylim = c(0.1, NA))
              
              exportPlots$plot_genefinder <- p
              p
            })
      
      output$rentrez_infobox <- renderUI({
        # Check if we have manually selected genes
        manual_genes <- selected_genes()
        if (!is.null(manual_genes) && length(manual_genes) > 0) {
          # Use the first manually selected gene
          selectedGene <- manual_genes[1]
          # Try to match by gene name first, then by ID
          matched_by_name <- which(values$annotation_obj$gene_name %in% selectedGene)
          if (length(matched_by_name) > 0) {
            selectedGene <- rownames(values$annotation_obj)[matched_by_name[1]]
          }
        } else {
          shiny::validate(
            need(
              (nrow(curDataClick()) > 0),
              "Select a gene first to display additional info (retrieved from the NCBI/ENTREZ db website)"
            )
          )
          selectedGene <- as.character(curDataClick()$ID)
        }
        
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
        # Check if we have manually selected genes
        manual_genes <- selected_genes()
        if (!is.null(manual_genes) && length(manual_genes) > 0) {
          # Create data frame with all manually selected genes
          mama <- data.frame(mean=values$res_obj$baseMean,lfc=values$res_obj$log2FoldChange,padj = values$res_obj$padj,isDE= ifelse(is.na(values$res_obj$padj), FALSE, values$res_obj$padj < 0.10),ID=rownames(values$res_obj))
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