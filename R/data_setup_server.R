data_setup_server <- function(id, values, annoSpecies_df) {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    # output dt_cm ----
    output$dt_cm <- DT::renderDataTable({
      shiny::validate(
        need(!is.null(values$countmatrix), "No count matrix data uploaded yet. Please upload a count matrix file to proceed.")
      )
      datatable(values$countmatrix, options = list(scrollX = TRUE, scrollY = "400px"))
    })
    
    # output dt_ed ----
    output$dt_ed <- DT::renderDataTable({
      shiny::validate(
        need(!is.null(values$expdesign), "No experimental design data uploaded yet. Please upload a metadata file to proceed.")
      )
      datatable(values$expdesign, options = list(scrollX = TRUE))
    })
    
    
    
    # output ui_step2 ----
    output$ui_step2 <- renderUI({
      cat(file = stderr(), paste("ui_step2\n"))
      shiny::validate(
        need(!is.null(values$expdesign) && !is.null(values$countmatrix),
             "Please provide both experimental design and count matrix files to proceed")
      )
      
      box(width = 12, title = "Step 2s", status = "warning", solidHeader = TRUE,
          tagList(
            # as in https://groups.google.com/forum/#!topic/shiny-discuss/qQ8yICfvDu0
            h2("Select the DE design and create the DESeqDataSet object"),
            fluidRow(
              column(
                width = 6,
                
                uiOutput(ns("ddsdesign")),
                uiOutput(ns("ddsintercept")),
                # checkboxInput(inputId = "multiplyDesign", label = "multiply first two arguments",value = FALSE),
                textInput(ns("geneFilter"), "Reg. expr. to fileter genes", value = values$geneFilter),
                uiOutput(ns("ui_diydds")),
                hr(),
                # uiOutput("ok_dds"),
                verbatimTextOutput(ns("debugdiy"))
              ),
              column(
                width = 6,
                plotOutput(ns("cooccurrence_matrix_plot")) %>% shinyjqui::jqui_resizable(),
                
              )
            )
          ))
    })
    
    
    # Plot cooccurrence matrix ----------------------------------------------
    # output cooccurrence_matrix_plot ----
    output$cooccurrence_matrix_plot <- shiny::renderPlot({
          shiny::validate(
            shiny::need(
              input$dds_design != "" ,
              "Please select a design formula where all terms appear in the sample data"
            )
          )
          # Suppress par() warnings that commonly occur in Shiny reactive contexts
          suppressWarnings({
            if (is.null(cooccurrenceplots())) {
              NULL
            } else {
              # cowplot::plot_grid(plotlist = cooccurrenceplots(),
              #                    ncol = 1)
              cooccurrenceplots()
            }
          })
        })
    
    # The following comes from ExploreModelMatrix package
    ## ----------------------------------------------------------------------- ##
    ## Create co-occurrence plot
    ## ----------------------------------------------------------------------- ##
    ## 
    cooccurrenceplots <- reactive({
      cat(file = stderr(), paste("cooccurrenceplots\n"))
      # browser()
      shiny::validate(
        need(!is.null(values$expdesign), "Experimental design data is not available"),
        need(!is.null(input$dds_design), "Please select a design formula"),
        need(input$dds_design[1] != "*", "Please select a valid design formula")
      )
      if(input$dds_design[1] == "*"){
        updateSelectInput(session, inputId = "dds_design", selected = "")
        shiny::validate(
          need(FALSE, "Please select a valid design formula")
        )
      }
      # if(!all(input$dds_design %in% colnames(values$expdesign)))
      #   return(NULL)
      dsgString = paste(input$dds_design, collapse = " + ")
      dsgString = str_replace(dsgString, fixed(" + * + "), " * ")
      shiny::validate(
        need(!is.null(tryCatch( as.formula(paste0("~", dsgString)),
                             error = function(e) return(NULL) )), "Invalid design formula. Please check your design selection.")
      )
      dsgn <- as.formula(paste0("~", dsgString))
      
      # Use the new utility function
      generate_cooccurrence_plots(values$expdesign, dsgn)
    })
    
    observeEvent(input$geneFilter, {
      cat(file = stderr(), paste("input$geneFilter\n"))
      # browser()
      if (is.null(values$geneFilter)) {
        updateTextInput(session, inputId = "geneFilter", value = "^MT-|^RP")
        return()
      }
      shiny::validate(
        need(input$geneFilter != values$geneFilter, "Gene filter has not changed")
      )
      values$geneFilter = input$geneFilter
      # browser()
      values$dds_obj = NULL
    })
    
    # output ddsdesign ----
    output$ddsdesign <- renderUI({
      cat(file = stderr(), paste("ddsdesign\n"))
      shiny::validate(
        need(!is.null(values$expdesign),
             "Please upload experimental design data to proceed with design selection")
      )
      poss_covars <- colnames(values$expdesign)
      cat(file = stderr(), "selectINPUT dds_design\n")
      # browser()
      choices = unique(c(values$dds_design, poss_covars, "*"))
      selectInput(ns('dds_design'), label = 'Select the design for your experiment: ',
                  choices = choices, selected = values$dds_design, multiple = TRUE, selectize = T)
    })
    
    
    observeEvent(input$help_format, {
      showModal(modalDialog(
        title = "Format specifications for WEIN",
        includeMarkdown(system.file("extdata", "datainput.md",package = "WEIN")),
        h4("Example:"),
        tags$img(
          src = base64enc::dataURI(file = system.file("www", "help_dataformats.png",package = "pcaExplorer"), mime = "image/png"),
          width = 750
        ),
        easyClose = TRUE,
        footer = NULL,
        size = "l"
      ))
    })
    
    observeEvent(input$btn_loaddemo,withProgress(
      message = "Loading demo data",
      detail = "Loading airway count and metadata information", value = 0,
      {
        aw <- requireNamespace("airway",quietly = TRUE)
        incProgress(0.2,detail = "`airway` package loaded")
        if(aw) {
          data(airway,package="airway",envir = environment())
          
          cm_airway <- assay(airway)
          incProgress(0.7, detail = "Count matrix loaded")
          ed_airway <- as.data.frame(colData(airway))
          
          values$countmatrix <- cm_airway
          values$expdesign <- ed_airway
          incProgress(0.3, detail = "Experimental metadata loaded")
          # just to be sure, erase the annotation and the rest
          # browser()
          values$dds_obj <- NULL
          values$annotation_obj <- NULL
          values$res_obj <- NULL
          showNotification("All components for generating the DESeqDataset object have been loaded, proceed to Step 2!",
                           type = "message")
        } else {
          showNotification("The 'airway' package is currently not installed. Please do so by executing BiocManager::install('airway') before launching WEIN()",type = "warning")
        }
      })
    )
    
    # output ddsintercept ----
    output$ddsintercept <- renderUI({
      shiny::validate(
        need(!is.null(values$expdesign) && !is.null(input$dds_design),
             "Please provide both experimental design data and select a design formula")
      )
      if(input$dds_design[1] == "*") {
        return(NULL)
      }
      poss_covars <- levels(values$expdesign[,input$dds_design[1]])
      selectInput(ns('dds_intercept'), label = 'Select the intercept for your experiment: ',
                  choices = c(NULL, poss_covars), selected = values$dds_intercept, multiple = FALSE)
    })
    
    observeEvent(input$dds_design, {
      if(length(input$dds_design) > 0 && input$dds_design[1] == "*") {
        updateSelectInput(session, ns("dds_design"), selected = "")
      }
    }, ignoreInit = TRUE)
    
    # output ui_stepanno ----
    output$ui_stepanno <- renderUI({
      cat(file = stderr(), paste("ui_stepanno\n"))
      shiny::validate(
        need(!is.null(values$dds_obj),
             "Please create a DESeqDataSet object first before removing outliers")
      )
      
      
      box(width = 12, title = "Optional Step", status = "info", solidHeader = TRUE,
          tagList(
            h2("Create the annotation data frame for your dataset"),
            
            fluidRow(
              column(
                width = 8,
                uiOutput(ns("ui_selectspecies")),
                verbatimTextOutput(ns("speciespkg")),
                uiOutput(ns("ui_idtype")),
                verbatimTextOutput(ns("printDIYanno"))
                
              )
            )
            ,
            uiOutput(ns("ui_getanno"))
          )
      )
    })
    # output printDIYanno ----
    output$printDIYanno <- renderPrint({
      print(head(values$annotation_obj))
    })
    
    
    output$speciespkg <- renderText({
      shiny::validate(
        need(!is.null(values$dds_obj), "DESeqDataSet object is not available. Please create it first."),
        need(values$cur_species!="",
             "Please select a species - this requires the corresponding annotation package to be installed"
        )
      )
      annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==values$cur_species]
      shiny::validate(
        need(require(annopkg,character.only=TRUE),
             paste0("The package ",annopkg, " is not installed/available. Try installing it with BiocManager::install('",annopkg,"')"))
      )
      retmsg <- paste0(annopkg," - package available and loaded")
      # if (!require(annopkg,character.only=TRUE)) {
      # stop("The package",annopkg, "is not installed/available. Try installing it with BiocManager::install() ?")
      # }
      retmsg <- paste0(retmsg," - ",gsub(".eg.db","",gsub("org.","",annopkg)))
      retmsg
    })
    
    
    
    # output ui_getanno ----
    output$ui_getanno <- renderUI({
      shiny::validate(
        need(!is.null(values$dds_obj), "DESeqDataSet object is not available. Please create it first."),
        need(values$cur_species != "",
             "Please select a species first in the panel above")
      )
      actionButton(ns("button_getanno"),"Retrieve the gene symbol annotation for the uploaded data", class = "btn btn-primary")
    })
    
    
    observeEvent(input$button_getanno,
                 {
                   withProgress(message="Retrieving the annotation...",
                                detail = "Locating package", value = 0,{
                                  # browser()
                                  annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
                                  incProgress(0.1,detail = "Matching identifiers")
                                  values$cur_species <- input$speciesSelect
                                  values$cur_type <- input$idtype
                                  
                                  if (input$idtype == "SYMBOL") {
                                    values$annotation_obj = data.frame(gene_id = rownames(values$dds_obj), gene_name = rownames(values$dds_obj), 
                                                                       stringsAsFactors = FALSE, row.names = rownames(values$dds_obj))
                                  } else {
                                    tryCatch({
                                      annotation_obj <- get_annotation_orgdb(values$dds_obj,orgdb_species = annopkg, idtype = input$idtype)
                                      values$annotation_obj <- annotation_obj
                                      # and also, set the species in the reactiveValues
                                    },
                                    error=function(e) {
                                      showNotification(
                                        paste("Warning! The annotation object was not generated,",
                                              "because of an error in the underlying `mapIds` function:",
                                              "-----", e), type = "warning")
                                    })
                                  }
                                })
                   
                 })
    
    observeEvent(input$speciesSelect,{
      values$cur_species <- input$speciesSelect
      curr_idtype <- values$cur_type
      updateSelectInput(session, inputId = ns("idtype"), selected = curr_idtype)
    }
    )
    
    observeEvent(input$button_rundeseq,
                 {
                   withProgress(message="Running DESeq on your data...",
                                detail = "This step might take a while", value = 0,{
                                  # trick to keep species info while still changing the dds_obj
                                  curr_species <- input$speciesSelect
                                  incProgress(0.1)
                                  
                                  # if(input$nrcores == 1){
                                  #   pa = FALSE
                                  #   bp = bpparam()
                                  # } else {
                                  pa = FALSE
                                  # Dynamically determine number of workers based on system
                                  workers <- parallel::detectCores() - 1
                                  # Ensure at least 1 worker
                                  workers <- max(1, workers)
                                  bp = MulticoreParam(workers = workers)
                                  # }
                                  # leave open option for computing in parallel?
                                  values$dds_obj <- tryCatch({
                                    DESeq(values$dds_obj,
                                          parallel = pa,
                                          BPPARAM = bp)},
                                    # BPPARAM = MulticoreParam(workers = input$nrcores))},
                                    error = function(e) {
                                      showNotification(
                                        paste(
                                          "Error during creation of DDS object",
                                          "-----", e
                                        ),
                                        type = "error"
                                      )
                                      return(values$dds_obj)
                                    }
                                  )
                                  incProgress(0.89)
                                  updateSelectInput(session, inputId = ns("speciesSelect"), selected = curr_species)
                                })
                 })
    
    # server run deseq --------------------------------------------------------
    # output rundeseq ----
    output$rundeseq <- renderUI({
      shiny::validate(
        need(!is.null(values$dds_obj), "DESeqDataSet object is not available. Please create it first.")
      )
      actionButton(ns("button_rundeseq"),"Run DESeq!", icon = icon("spinner"), class = "btn btn-success")
    })
    
    output$ui_step3 <- renderUI({
      shiny::validate(
        need(!is.null(values$dds_obj), "DESeqDataSet object is not available. Please create it first.")
      )
      box(width = 12, title = "Step 3", status = "success", solidHeader = TRUE,
          tagList(
            h2("Run DESeq!"),
            
            # fluidRow(
            #   column(
            #     width = 4,
            #     uiOutput("ui_nrcores")
            #   )
            # ),
            
            uiOutput(ns("rundeseq")),
            verbatimTextOutput(ns("printDIYresults")),
            uiOutput(ns("ui_stepend"))
          )
      )
    })
    
    # output ui_stepend ----
    output$ui_stepend <- renderUI({
      shiny::validate(
        need(!is.null(values$dds_obj), "DESeqDataSet object is not available. Please create it first."),
        need("results" %in% mcols(mcols(values$dds_obj))$type, "DESeq analysis results are not available. Please run DESeq analysis first.")
      )
      
      tagList(
        h2("Good to go!")
        # ,
        # box(width = 6, title = "Diagnostic plot", status = "info", solidHeader = TRUE,
        #     collapsible = TRUE, collapsed = TRUE,
        #     plotOutput("diagno_dispests"))
      )
    })
    
    
    # output printDIYresults ----
    output$printDIYresults <- renderPrint({
      shiny::validate(
        need(!is.null(values$dds_obj),
             "Please provide or construct a DESeqDataSet object to view results")
      )
      shiny::validate(
        need("results" %in% mcols(mcols(values$dds_obj))$type ,
             "DESeq results not found. Please run DESeq() first using the 'Run DESeq!' button above"
        )
      )
      summary(DESeq2::results(values$dds_obj), alpha = values$FDR)
    })
    
    # output ui_selectspecies ----
    output$ui_selectspecies <- renderUI({
      shiny::validate(
        need(!is.null(values$dds_obj), "DESeqDataSet object is not available. Please create it first.")
      )
      selectInput(ns("speciesSelect"),label = "Select the species of your samples - it will also be used for enhancing result tables",
                  choices = annoSpecies_df$species,selected="Human")
    })
    # output ui_idtype ----
    output$ui_idtype <- renderUI({
      shiny::validate(
        need(!is.null(values$dds_obj), "DESeqDataSet object is not available. Please create it first.")
      )
      
      std_choices <- c("ENSEMBL","ENTREZID","REFSEQ","SYMBOL")
      if(!"speciesSelect" %in% names(input)) {
        return (NULL)
      }
      if (input$speciesSelect!=""){
        annopkg <- annoSpecies_df$pkg[annoSpecies_df$species==input$speciesSelect]
        require(annopkg,character.only=TRUE)
        pkg_choices <- keytypes(get(annopkg))
        std_choices <- union(std_choices, pkg_choices)
      }
      selectInput(ns("idtype"), "select the id type in your data", choices=std_choices, selected = "SYMBOL")
    })
    
    # output ui_diydds ----
    output$ui_diydds <- renderUI({
      shiny::validate(
        need(!is.null(values$expdesign) && !is.null(values$countmatrix) && !is.null(input$dds_design),
             "Please provide experimental design, count matrix, and select a design formula to generate the dds object")
      )
      actionButton(ns("button_diydds"),"Generate the dds object", class = "btn btn-success")
    })
    
    observeEvent(input$button_diydds,{
      if(!is.null(values$countmatrix) & !is.null(values$expdesign))
        values$dds_obj <- diyDDS()
    })
    
    # output ui_stepoutlier ----
    output$ui_stepoutlier <- renderUI({
      shiny::validate(
        need(!is.null(values$dds_obj),
             "Please create a DESeqDataSet object first by uploading data and selecting a design formula")
      )
      # browser()
      box(
        width = 12, title = "Optional Step", status = "info", solidHeader = TRUE,
        tagList(
          h2("Remove sample(s) from the current dataset - suspected outliers!"),
          
          fluidRow(
            column(
              width = 8,
              uiOutput(ns("ui_selectoutliers")),
              uiOutput(ns("outliersout")),
              verbatimTextOutput(ns("printremoved"))
            )
          )
        )
      )
    })
    
    output$printremoved <- renderPrint({
      print(values$removedsamples)
    })
    
    # output outliersout ----
    output$outliersout <- renderUI({
      shiny::validate(
        need(!is.null(values$dds_obj), "DESeqDataSet object is not available. Please create it first.")
      )
      actionButton(ns("button_outliersout"),"Recompute the dds without some samples",class = "btn btn-primary")
    })
    
    observeEvent(input$button_outliersout,{
      cat(file = stderr(), paste("input$button_outliersout\n"))
      
      withProgress({
        allsamples <- colnames(values$countmatrix)
        outliersamples <- input$selectoutliers
        
        keptsamples <- setdiff(allsamples,outliersamples)
        colData <- values$expdesign[keptsamples,]
        colData[, input$dds_design[1]] = relevel(colData[, input$dds_design[1]], ref = values$dds_intercept)
        dds <- tryCatch({DESeqDataSetFromMatrix(countData = values$countmatrix[,keptsamples],
                                                colData = colData,
                                                design  = design(values$dds_obj)
                                                # design=as.formula(paste0("~",paste(input$dds_design, collapse=" + ")))
        )}, error = function(e) {
          showNotification(
            paste(
              "Error during creation of DDS object",
              "-----", e
            ),
            type = "error"
          )
          return(NULL)
        })
        if (is.null(dds)) {
          return(NULL)
        }
        dds <- estimateSizeFactors(dds)
        
        # return(dds)
        # re-create the dds and keep track of which samples were removed
        values$removedsamples <- input$selectoutliers
        
        curr_species <- input$speciesSelect
        values$dds_obj <- dds
        updateSelectInput(session, inputId = ns("speciesSelect"), selected = curr_species)
        
        # accordingly, reset the results
        values$res_obj <- NULL},
        message = "Removing selected samples from the current dataset")
    })
    
    
    # output upload_count_matrix ----
    output$upload_count_matrix <- renderUI({
      if (values$input_provided) {
        return(fluidRow(column(
          width = 12,
          tags$li("You already provided a count matrix or a DESeqDataSet object as input. You can check your input data in the collapsible box here below."), offset = 2)))
      } else {
        return(fileInput(inputId = ns("uploadcmfile"),
                         label = "Upload one or more count matrix file(s)",
                         accept = c("text/csv", "text/comma-separated-values",
                                    "text/tab-separated-values", "text/plain",
                                    ".csv", ".tsv", ".xls"), multiple = TRUE))
      }
    })
    
    output$upload_metadata_ui <- renderUI({
      if (values$input_provided) {
        return(fluidRow(column(
          width = 12,
          tags$li("You already provided a matrix/data.frame with the experimental covariates or a DESeqDataSet object as input. You can check your input data in the collapsible box here below."), offset = 2)))

      } else {
        return(fileInput(inputId = ns("uploadmetadatafile"),
                         label = "Upload a sample metadata matrix file",
                         accept = c("text/csv", "text/comma-separated-values",
                                    "text/tab-separated-values", "text/plain",
                                    ".csv", ".tsv"), multiple = FALSE))
      }
    })
    
    # output upload_metadata ----
    output$upload_metadata <- renderUI({
      shiny::validate(
        need(!is.null(input$uploadcmfile), "Please upload a count matrix file first.")
      )
      on.exit({
        if (!is.null(getDefaultReactiveDomain())) {
          removeNotification(id = "readCountmatrix")
        }
      })
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("readCountmatrix",
                         id = "readCountmatrix",
                         duration = NULL
        )
      }
      
      # Use the new utility function
      cm <- read_countmatrix(input$uploadcmfile$datapath[1])
      if (nrow(input$uploadcmfile) > 1){
        for (r in 2:nrow(input$uploadcmfile)) {
          # Use the new utility function
          cmt <- read_countmatrix(input$uploadcmfile$datapath[r])
          cm2 = base::merge(cm, cmt, by = 0,  all=T)
          rownames(cm2) = cm2$Row.names
          cm2$Row.names = NULL
          cm = cm2
        }
      }
      cat(file = stderr(), paste("read count data: rows:", nrow(cm), "ncol:", ncol(cm), "\n", "sample names:", colnames(cm), "\n"))
      # browser()
      return(cm)
    })
    
    readCountmatrix <- reactive({
      shiny::validate(
        need(!is.null(input$uploadcmfile), "Please upload a count matrix file first.")
      )
      on.exit({
        if (!is.null(getDefaultReactiveDomain())) {
          removeNotification(id = "readCountmatrix")
        }
      })
      if (!is.null(getDefaultReactiveDomain())) {
        showNotification("readCountmatrix",
                         id = "readCountmatrix",
                         duration = NULL
        )
      }
      
      guessed_sep <- sepguesser(input$uploadcmfile$datapath[1])
      cm <- utils::read.delim(input$uploadcmfile$datapath[1], header = TRUE,
                              as.is = TRUE, sep = guessed_sep, quote = "",
                              row.names = 1, # https://github.com/federicomarini/pcaExplorer/issues/1
                              ## TODO: tell the user to use tsv, or use heuristics
                              ## to check what is most frequently occurring separation character? -> see sepGuesser.R
                              check.names = FALSE)
      if (nrow(input$uploadcmfile) > 1){
        for (r in 2:nrow(input$uploadcmfile)) {
          cmt <- utils::read.delim(input$uploadcmfile$datapath[r], header = TRUE,
                                   as.is = TRUE, sep = guessed_sep, quote = "",
                                   row.names = 1, # https://github.com/federicomarini/pcaExplorer/issues/1
                                   ## TODO: tell the user to use tsv, or use heuristics
                                   ## to check what is most frequently occurring separation character? -> see sepGuesser.R
                                   check.names = FALSE)
          cm2 = base::merge(cm, cmt, by = 0,  all=T)
          rownames(cm2) = cm2$Row.names
          cm2$Row.names = NULL
          cm = cm2
        }
      }
      cm[is.na(cm)] <- 0
      cat(file = stderr(), paste("read count data: rows:", nrow(cm), "ncol:", ncol(cm), "\n", "sample names:", colnames(cm), "\n"))
      # browser()
      return(cm)
    })
    
    
    readMetadata <- reactive({
      shiny::validate(
        need(!is.null(input$uploadmetadatafile), "Please upload a metadata file first.")
      )
      # browser()
      
      # Use the new utility function
      expdesign <- read_metadata(input$uploadmetadatafile$datapath)
      cat(file = stderr(), paste("read metadata: rows:", nrow(expdesign), "colnames:", colnames(expdesign), "\n"))
      # browser()
      return(expdesign)
    })
    
    # as in http://stackoverflow.com/questions/29716868/r-shiny-how-to-get-an-reactive-data-frame-updated-each-time-pressing-an-actionb
    observeEvent(input$uploadcmfile,
                 {
                   # browser()
                   values$countmatrix <- readCountmatrix()
                   values$dds_obj <- NULL
                   values$res_obj <- NULL
                 })
    
    observeEvent(input$uploadmetadatafile,{
      cat(file = stderr(), paste("input$uploadmetadatafile\n"))
      
      # browser()
      if("restoreBookmark" %in% names(values))
        if(values$restoreBookmark & !is.null(values$expdesign)){
          return()
        }
      
      values$expdesign <- readMetadata()
      values$dds_obj <- NULL
      values$res_obj <- NULL
    })
    
    observeEvent(input$dds_intercept, {
      values$dds_intercept = input$dds_intercept
    })
    
    # server outliers --------------------------------------------------------
    # output ui_selectoutliers ----
    output$ui_selectoutliers <- renderUI({
      shiny::validate(
        need(!is.null(values$dds_obj), "DESeqDataSet object is not available. Please create it first.")
      )
      # cat(file = stderr(), "ui_selectoutliers UI\n")
      md = colnames(values$dds_obj)
      cm = colnames(values$dds_obj)
      if (!is.null(readCountmatrix()))
        cm <- readCountmatrix()
      if (!is.null(readMetadata()))
        md <- readMetadata()
      comSamples <- intersect(colnames(cm), rownames(md))
      sele <- values$removedsamples
      
      
      checkboxGroupInput("selectoutliers","Select the samples to remove - candidate outliers",
                         choices = comSamples, selected = sele)
      # selectInput("selectoutliers","Select the samples to remove - candidate outliers",
      #             choices = colnames(values$dds_obj), selected = NULL,multiple = TRUE
      # )
    })
    
    
    diyDDS <- reactive({
      cat(file = stderr(), paste("diyDDS\n"))
      
      shiny::validate(
        need(!is.null(values$countmatrix) & !is.null(values$expdesign) & !is.null(input$dds_design),
             "Please provide count matrix, experimental design, and select a design formula to generate the dds object.")
      )
      # browser()
      if (!is.null(readCountmatrix()))
        values$countmatrix <- readCountmatrix()
      if (!is.null(readMetadata()))
        values$expdesign <- readMetadata()
      comSamples <- intersect(colnames(values$countmatrix), rownames(values$expdesign))
      if(length(comSamples)==0){
        showNotification("count and meta data don't have any samples in common. Please reload data.", type = "error")
      }
      metaData <- readMetadata()
      values$expdesign <- values$expdesign[comSamples, ]
      dsgString = paste(input$dds_design, collapse = " + ")
      dsgString = str_replace(dsgString, fixed(" + * + "), " * ")
      shiny::validate(
        need(!is.null(tryCatch( as.formula(paste0("~", dsgString)),
                             error = function(e) return(NULL) )), "Invalid design formula. Please check your design selection.")
      )
      dsgn <- as.formula(paste0("~", dsgString))
      
      # dsgn <- input$dds_design
      filterExp <- input$geneFilter
      values$geneFilter = input$geneFilter
      values$dds_design = input$dds_design
      values$res_obj = NULL
      if (!is(dsgn,"formula")) {
        dStr <- paste0("~", paste(input$dds_design, collapse = " + "))
        # if (input$multiplyDesign) {
        #   dStr <- sub('\\+', '*', dStr)
        # }
        dsgn <- as.formula(dStr)
      }
      locfunc <- stats::median
      counts <- values$countmatrix[, comSamples]
      # save(file = "~/SCHNAPPsDebug/WEIN.RData", list = ls())
      # cp = load("~/SCHNAPPsDebug/WEIN.RData")
      if (nchar(filterExp)>0){
        counts <- counts[grep(filterExp, rownames(counts), invert = TRUE), ]
      }
      values$countmatrix <- counts
      cGenes = 1:nrow(counts)
      # browser()
      md = readMetadata()
      rc = readCountmatrix()
      colData = values$expdesign[comSamples, ]
      c1 = colnames(values$countmatrix)
      c2= rownames(values$expdesign)
      design = dsgn
      # browser()
      # save(file = "~/SCHNAPPsDebug/WEINDDS.RData", list = c('counts', "colData", "design", "comSamples", "c1", "c2", "md", "rc"))
      
      # Use the new utility function
      dds <- tryCatch(
        {
          dds <- create_dds(
            countmatrix = counts,
            expdesign = colData,
            design_formula = design,
            gene_filter = filterExp,
            dds_intercept = values$dds_intercept,
            design_factor = input$dds_design[1]
          )
        },
        error = function(e) {
          cat(file = stderr(), paste("error during creation of dds object:", e))
          # save(file = "~/SCHNAPPsDebug/WEINDDS.RData", list = c('counts', "colData", "design", "comSamples", "c1", "c2", "md", "rc"))
          # cp =load("/Users/bernd/SCHNAPPsDebug/WEINDDS.RData")
          showNotification(
            paste(
              "Error during creation of DDS object",
              "-----", e
            ),
            type = "error"
          )
          shiny::validate(
            need(FALSE, paste("Error during creation of DDS object:", e))
          )
        }
      )
      return(dds)
    })
    
    # output debugdiy ----
    output$debugdiy <- renderPrint({
      shiny::validate(
        need(!is.null(values$dds_obj), "DESeqDataSet object is not available. Please create it first.")
      )
      print(values$dds_obj)
      print("Design:")
      print(design(values$dds_obj))
      print("Gene filter:")
      print(values$geneFilter)
    })
    
    return(list(
      genefilter = reactive(input$geneFilter),
      design   = reactive(input$dds_design),
      intercept= reactive(input$dds_intercept)
      # ,
      # dds_obj  = dds_obj,
      # diyDDS   = diyDDS
    ))
  })
}

