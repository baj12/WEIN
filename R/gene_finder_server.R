gene_finder_server <- function(id, values, annoSpecies_df, exportPlots) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    # Helper function (define at top of server)
    has_valid_genes <- function(x) {
      if (is.null(x)) return(FALSE)
      if (length(x) == 0) return(FALSE)
      if (all(is.na(x))) return(FALSE)
      if (all(x == "")) return(FALSE)
      return(TRUE)
    }
    
    # Validates gene availability at index
    validate_gene_at_index <- function(values, index) {
      ordinal <- c("first", "second", "third", "fourth")[index]
      
      shiny::validate(
        need(
          has_valid_genes(values$avail_symbols) || has_valid_genes(values$avail_ids),
          paste("Select at least", ifelse(index == 1, "a", paste("a", ordinal)), "gene to plot")
        ),
        need(
          length(values$avail_symbols) >= index || length(values$avail_ids) >= index,
          paste("Select at least", ifelse(index == 1, "a", paste("a", ordinal)), "gene to plot")
        )
      )
    }
    validate_color_by <- function(values) {
      shiny::validate(
        need(
          length(values$color_by) > 0,
          "Select an experimental factor in the Group/color by element in the sidebar"
        )
      )
    }
    
    # Gets gene ID at specified index, handling symbol/id conversion
    get_gene_id_at_index <- function(values, index) {
      if (length(values$avail_symbols) >= index) {
        # Got symbol, look for ID
        mysym <- values$avail_symbols[index]
        myid <- values$annotation_obj$gene_id[match(mysym, values$annotation_obj$gene_name)]
      } else {
        # Got ID, optionally look for symbol
        myid <- values$avail_ids[index]
        if (!is.null(values$annotation_obj)) {
          mysym <- values$annotation_obj$gene_name[match(myid, values$annotation_obj$gene_id)]
        } else {
          mysym <- ""
        }
      }
      
      return(myid)
    }
    create_gene_plot <- function(values, input, index) {
          validate_color_by(values)
          validate_gene_at_index(values, index)
          
          myid <- get_gene_id_at_index(values, index)
          
          # Suppress par() warnings that commonly occur in Shiny reactive contexts
          p <- suppressWarnings({
            ggplotCounts(
              values$dds_obj,
              myid,
              intgroup = values$color_by,
              annotation_obj = values$annotation_obj
            )
          })
          
          if (input$ylimZero_genefinder) {
            p <- p + ylim(0.1, NA)
          }
          
          return(p)
        }
    
    # Create all 4 plot outputs using loop
    lapply(1:4, function(i) {
      output_name <- paste0("bp", i)
      export_name <- paste0("plotbp", i)
      
      output[[output_name]] <- renderPlot({
        p <- create_gene_plot(values, input, i)
        exportPlots[[export_name]] <- p
        p
      })
    })
    
    
    cur_combires_list <- reactive({
      if(is.null(values$res_obj))
        return(NULL)
      
      normCounts <- as.data.frame(counts(estimateSizeFactors(values$dds_obj),normalized=TRUE))
      normCounts$id <- rownames(normCounts)
      res_df <- deseqresult2tbl(values$res_obj)
      
      combi_obj <- dplyr::inner_join(res_df,normCounts,by="id")
      combi_obj$symbol <- values$annotation_obj$gene_name[match(combi_obj$id,values$annotation_obj$gene_id)]
      
      
      if("symbol" %in% names(values$res_obj)) {
        sel_genes <- values$genelist_ma$`Gene Symbol`
        sel_genes_ids <- values$annotation_obj$gene_id[match(sel_genes,values$annotation_obj$gene_name)]
      } else {
        # sel_genes_ids <- values$genelist_ma$`Gene Symbol`
      }
      
      if(length(sel_genes_ids) > 0) {
        combi_obj[match(sel_genes_ids,combi_obj$id),]
      } else {
        combi_obj
      }
    })
    

    output$plotCoefficients<- renderPlot({
          shiny::validate(
            need(
              length(values$color_by)>0,
              "Select an experimental factor in the Group/color by element in the sidebar"
            )
          )
          shiny::validate(
            need(
              (length(values$avail_symbols)>0 | length(values$avail_ids)>0),
              "Select at least a gene to plot"
            )
          )
          # browser()
          mysym <- values$avail_symbols[1]
          myid <- values$annotation_obj$gene_id[match(mysym, values$annotation_obj$gene_name)]
          # Suppress par() warnings that commonly occur in Shiny reactive contexts
          suppressWarnings({
            plotCoefficients(values$dds_obj, myid, legend = T)
          })
          
          
        })
    
    # server report editor --------------------------------------------------------
    ### yaml generation
    rmd_yaml <- reactive({
      paste0("---",
             "\ntitle: '", input$report_title,
             "'\nauthor: '", input$report_author,
             "'\ndate: '", Sys.Date(),
             "'\noutput:\n  html_document:\n    toc: ", input$report_toc, "\n    number_sections: ", input$report_ns, "\n    theme: ", input$report_theme, "\n---\n\n",collapse = "\n")
    })
    
    output$ma_highlight <- renderPlot({
          shiny::validate(
            need(!is.null(values$res_obj),message = "Please generate the results object in the Extract Results panel to display the plot and show the combined tables")
          )
          
          # Suppress par() warnings that commonly occur in Shiny reactive contexts
          suppressWarnings({
            if("symbol" %in% names(values$res_obj)) {
              p <- plot_ma(values$res_obj,
                           intgenes = values$avail_symbols,annotation_obj = values$annotation_obj,FDR = values$FDR)
            } else {
              p <- plot_ma(values$res_obj,
                           intgenes = values$avail_ids,annotation_obj = values$annotation_obj,FDR = values$FDR)
            }
          })
          
          exportPlots$plot_mahighlight <- p
          p
        })
    
    output$ma_hl_list <- renderPlot({
          shiny::validate(
            need(!is.null(values$genelist_ma),message = "Please select genes in the MA plot to generate a list to plot here")
          )
          shiny::validate(
            need("symbol" %in% names(values$res_obj),
                 message = "Please ensure your results object has gene symbol annotation. This requires setting up annotation in the Data Setup panel.")
          )
          if(is.null(values$genelist_ma))
            return(NULL)
          
          # Suppress par() warnings that commonly occur in Shiny reactive contexts
          suppressWarnings({
            if("symbol" %in% names(values$res_obj)) {
              p <- plot_ma(values$res_obj,
                           intgenes = values$genelist_ma$`Gene Symbol`,annotation_obj = values$annotation_obj,FDR = values$FDR)
            } else {
              # plot_ma(values$res_obj,
              # intgenes = values$genelist_ma,annotation_obj = values$annotation_obj)
              return(NULL)
            }
          })
          exportPlots$plot_mahllist <- p
          p
        })
    
    observeEvent(input$gl_ma,
                     {
                       gl = gl_ma()
                       if(is.null(gl)) {values$genelist_ma = data.frame(); return(NULL)}
                       if(nrow(gl)<1) {values$genelist_ma = data.frame(); return(NULL)}
                       # If gl is already a data frame with gene_id column, convert to Gene Symbol
                       if("gene_id" %in% names(gl)) {
                         mydf <- data.frame("Gene Symbol" = gl$gene_id, stringsAsFactors=FALSE)
                       } else {
                         mydf <- as.data.frame(gl,stringsAsFactors=FALSE)
                         names(mydf) <- "Gene Symbol"
                       }
                       values$genelist_ma <- mydf
                     })
    
    
    
    gl_ma <- reactive({
      if (is.null(input$gl_ma)) {
        # User has not uploaded a file yet
        return(data.frame())
      } else {
        gl_ma <- read1stCol(input$gl_ma$datapath, values$dds_obj)
        # browser()
        return(gl_ma)
      }
    })
    
    output$download_plotbp1 <- downloadHandler(filename = function() {
      input$filename_plotbp1
    }, content = function(file) {
      ggsave(file, exportPlots$plotbp1, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    
    output$download_plotbp2 <- downloadHandler(filename = function() {
      input$filename_plotbp2
    }, content = function(file) {
      ggsave(file, exportPlots$plotbp2, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    
    output$download_plotbp3 <- downloadHandler(filename = function() {
      input$filename_plotbp3
    }, content = function(file) {
      ggsave(file, exportPlots$plotbp3, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    
    output$download_plotbp4 <- downloadHandler(filename = function() {
      input$filename_plotbp4
    }, content = function(file) {
      ggsave(file, exportPlots$plotbp4, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    
    output$downloadTblCombi <- downloadHandler(
      filename = function() {
        "table_combi.csv"
      },
      content = function(file) {
        write.csv(cur_combires(), file)
      }
    )
    output$downloadTblCombiList <- downloadHandler(
      filename = function() {
        "table_combilist.csv"
      },
      content = function(file) {
        write.csv(cur_combires_list(), file)
      }
    )
    output$download_plot_mahighlight <- downloadHandler(filename = function() {
      input$filename_plot_mahighlight
    }, content = function(file) {
      ggsave(file, exportPlots$plot_mahighlight, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    
    
    output$download_plot_mahllist <- downloadHandler(filename = function() {
      input$filename_plot_mahllist
    }, content = function(file) {
      ggsave(file, exportPlots$plot_mahllist, width = values$export_width,
             height = values$export_height, units = "cm")
    })
    output$table_combi_list <- DT::renderDataTable({
      if(is.null(values$genelist_ma))
        return(NULL)
      datatable(cur_combires_list(),options = list(scrollX=TRUE))
    })
    
  })
}
