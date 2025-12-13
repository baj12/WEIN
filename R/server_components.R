#' Server Components for WEIN
#'
#' Contains all server-related functions for the WEIN Shiny application.
#'
#' @name WEIN-server
#' @docType package
#' @keywords internal
"_PACKAGE"


#' Create the server for WEIN
#'
#' Defines the complete server logic for the WEIN Shiny application.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param dds_obj DESeqDataSet object
#' @param res_obj DESeqResults object
#' @param annotation_obj Annotation data frame
#' @param countmatrix Count matrix
#' @param expdesign Experimental design data frame
#' @param gene_signatures Gene signatures list
#' @param dds_design Design formula for DESeq2
#' @param cur_species Current species for annotation
#' @param cur_type Current ID type for annotation
#'
#' @return A Shiny server definition
#' @export
WEIN_server <- function(
    app_data = list(
      dds_obj = NULL,
      res_obj = NULL,
      annotation_obj = NULL,
      countmatrix = NULL,
      expdesign = NULL,
      gene_signatures = NULL,
      dds_design = NULL,
      cur_species = NULL,
      cur_type = NULL))
{
  dds_obj <- app_data$dds_obj
  res_obj <- app_data$res_obj
  annotation_obj <- app_data$annotation_obj
  countmatrix <- app_data$countmatrix
  expdesign <- app_data$expdesign
  gene_signatures <- app_data$gene_signatures
  dds_design <- app_data$dds_design
  cur_species <- app_data$cur_species
  cur_type <- app_data$cur_type
  
  function(input, output, session){
    require(ashr)
    require(airway)
    require(shinyjqui)
    require(apeglm)
    require(DESeq2) 
    require(SummarizedExperiment) 
    require(GenomicRanges) 
    require(IRanges)
    require(S4Vectors) 
    require(ggplot2) 
    require(heatmaply) 
    require(plotly)
    require(pcaExplorer) 
    require(IHW) 
    require(gplots) 
    require(UpSetR) 
    require(goseq) 
    require(stringr) 
    require(dplyr)
    require(limma) 
    require(GOstats) 
    require(GO.db) 
    require(AnnotationDbi) 
    require(shiny)
    require(shinydashboard) 
    require(shinyBS) 
    require(DT) 
    require(rentrez) 
    require(rintrojs) 
    require(ggrepel) 
    require(knitr)
    require(rmarkdown) 
    require(shinyAce) 
    require(BiocParallel) 
    require(grDevices) 
    require(base64enc)
    require(methods)
    require(testthat) 
    require(BiocStyle) 
    require(airway) 
    require(org.Hs.eg.db)
    require(TxDb.Hsapiens.UCSC.hg38.knownGene)
    require(DEFormats)
    require(edgeR)
    require(reshape2)
    require(cowplot)
    
    # server setup reactivevalues -----------------------------------------------------------
    ## placeholder for the figures to export
    exportPlots <- reactiveValues()
    # will store all the reactive values relevant to the app
    values <- reactiveValues(
      countmatrix = countmatrix,
      expdesign = expdesign,
      dds_obj = dds_obj,
      res_obj = res_obj,
      annotation_obj = annotation_obj,
      gene_signatures = gene_signatures,
      dds_design = dds_design,
      cur_species = cur_species,
      cur_type = cur_type,
      color_by = "",
      avail_symbols = "",
      FDR = 0.05,
      avail_ids = "",
      export_width = 16,
      export_height = 10,
      genelistUP = c(),
      genelistDOWN = c(),
      genelistUPDOWN = c(),
      choose_expfac = NULL,
      choose_expfac2 = NULL,
    )
    
    # Check if we're restoring from a saved state file
    # If dds_obj is provided but countmatrix/expdesign are not, populate them
    if(!is.null(dds_obj) && (is.null(countmatrix) || is.null(expdesign))) {
      values$countmatrix <- counts(dds_obj, normalized = FALSE)
      values$expdesign <- as.data.frame(colData(dds_obj))
    }
    
    # this part sets the "matching" objects if something is provided that is depending on these
    if(!is.null(dds_obj)){
      values$countmatrix <- counts(dds_obj, normalized = FALSE)
      values$expdesign <- as.data.frame(colData(dds_obj))
    }
    
    ## Update directory
    userdir <- tempfile()
    dir.create(userdir, recursive = TRUE)
    
    # server tours setup -----------------------------------------------------------
    tourModuleServer("tour_manager", input, session)
    
    ### bookmarking --------------------------------------------------------------------------
    onBookmark(function(state) {
      # browser()
      state$values$ValList = list()
      valList = reactiveValuesToList(values)
      for(na in names(valList)){
        state$values$ValList[[na]]<- valList[[na]]
      }
    })
    
    onRestore(function(state) {
      # browser()
      vList = state$values$ValList
      for(na in names(vList)){
        values[[na]] = state$values$ValList[[na]]
      }
      values$restoreBookmark = TRUE
      dds_obj <<- values$dds_obj
      countmatrix <<- values$countmatrix 
      expdesign <<- values$expdesign 
      
      res_obj <<- values$res_obj
      annotation_obj <<- values$annotation_obj
      gene_signatures <<- values$gene_signatures
      dds_design <<- values$dds_design
      
    })
    
    onRestored(function(state){
      # browser()
      # Update UI inputs to match restored values
      if (!is.null(values$cur_species)) {
        updateSelectInput(session, "ui_outputs_manager-speciesSelect", selected = values$cur_species)
      }
      if (!is.null(values$cur_type)) {
        updateSelectInput(session, "ui_outputs_manager-idtype", selected = values$cur_type)
      }
      if (!is.null(values$dds_intercept)) {
        updateSelectInput(session, "ui_outputs_manager-dds_intercept", selected = values$dds_intercept)
      }
      if (!is.null(values$dds_design)) {
        updateSelectInput(session, "ui_outputs_manager-dds_design", selected = values$dds_design)
      }
      if (!is.null(values$geneFilter)) {
        updateTextInput(session, "ui_outputs_manager-geneFilter", value = values$geneFilter)
      }
      if (!is.null(values$FDR)) {
        updateNumericInput(session, "FDR", value = values$FDR)
      }
      if (!is.null(values$export_width)) {
        updateNumericInput(session, "export_width", value = values$export_width)
      }
      if (!is.null(values$export_height)) {
        updateNumericInput(session, "export_height", value = values$export_height)
      }
      # Try to update color_by if it exists in values
      if (!is.null(values$color_by) && length(values$color_by) > 0) {
        tryCatch({
          updateSelectInput(session, "color_by", selected = values$color_by)
        }, error = function(e) {
          # Ignore errors if the input doesn't exist or isn't ready yet
        })
      }
      values$restoreBookmark = FALSE
    })
    
    
    # if i want to focus a little more on the ihw object
    values$ihwres <- NULL
    
    
    
    # server ui steps -----------------------------------------------------------
    
    # server retrieving anno --------------------------------------------------
    annoSpecies_df <- 
      data.frame(species=c("","Anopheles","Arabidopsis","Bovine","Worm",
                           "Canine","Fly","Zebrafish","E coli strain K12",
                           "E coli strain Sakai","Chicken","Human","Mouse",
                           "Rhesus","Malaria","Chimp","Rat",
                           "Yeast","Streptomyces coelicolor", "Pig","Toxoplasma gondii",
                           "Xenopus"),
                 pkg=c("","org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db",
                       "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db",
                       "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db",
                       "org.Mmu.eg.db", "org.Pf.plasmo.db", "org.Pt.eg.db", "org.Rn.eg.db",
                       "org.Sc.sgd.db", "org.Sco.eg.db", "org.Ss.eg.db", "org.Tgondii.eg.db",
                       "org.Xl.eg.db"),
                 stringsAsFactors = FALSE)
    
    annoSpecies_df <- annoSpecies_df[order(annoSpecies_df$species),]
    # this one is relevant for creating links to the genes
    annoSpecies_df$ensembl_db <- c("","","","Bos_taurus","Canis_familiaris","Gallus_gallus","Pan_troglodytes",
                                   "","","Drosophila_melanogaster","Homo_sapiens","","Mus_musculus",
                                   "Sus_scrofa","Rattus_norvegicus","Macaca_mulatta","","","Caenorhabditis_elegans",
                                   "Xenopus_tropicalis","Saccharomyces_cerevisiae","Danio_rerio"
    )
    # this one is the shortcut for the limma::goana function
    annoSpecies_df$species_short[grep(pattern = "eg.db",annoSpecies_df$pkg)] <- gsub(".eg.db","",gsub("org.","",annoSpecies_df$pkg))[grep(pattern = "eg.db",annoSpecies_df$pkg) ]
    # to match to the goseq genome setting
    annoSpecies_df$goseq_shortcut <- c("","anoGam1","Arabidopsis","bosTau8","canFam3","galGal4","panTro4","E. coli K12","E. coli Sakai",
                                       "dm6","hg19","Malaria","mm10","susScr3","rn6","rheMac","","","ce11","xenTro","sacCer3","danRer10")
    rownames(annoSpecies_df) <- annoSpecies_df$species # easier to access afterwards
    # annoSpecies_df <- annoSpecies_df[annoSpecies_df$species %in% c("","Human", "Mouse", "Rat", "Fly", "Chimp"),]
    
    
    
    
    # not implemented
    output$diagno_dispests <- renderPlot({
          # Suppress par() warnings that commonly occur in Shiny reactive contexts
          suppressWarnings({
            BiocGenerics::plotDispEsts(values$dds_obj)
          })
        })
    
    
    output$checkdds <- reactive({
      is.null(values$dds_obj)
    })
    output$checkresu<-reactive({
      is.null(values$res_obj)
    })
    
    outputOptions(output, 'checkresu', suspendWhenHidden = FALSE)
    outputOptions(output, 'checkdds', suspendWhenHidden = FALSE)
    
    # manager calls ----
    ui_setup_server("ui_setup", values)  
    data_setup_server("ui_outputs_manager", values, annoSpecies_df)
    count_overview_server("count_overview_manager", values)
    extract_results_server("extract_results_manager", values, annoSpecies_df, exportPlots)
    summary_plots_server("summary_plots", values, annoSpecies_df, exportPlots) 
    gene_finder_server("gene_finder", values, annoSpecies_df, exportPlots) 
    functional_analysis_server("fun_ana", values, annoSpecies_df, exportPlots)
    signature_explorer_server("sig_expl", values, annoSpecies_df, exportPlots) 
      
    # observers for values reactive
    observeEvent(input$color_by,
                 {values$color_by = input$color_by})
    
    observeEvent(input$avail_symbols,
                 {values$avail_symbols = input$avail_symbols}
    )
    
    observeEvent(input$FDR,
                 {values$FDR <- input$FDR}
    )
    observeEvent(input$export_width,
                 {values$export_width <- input$export_width}
    )
    observeEvent(input$export_height,
                 {values$export_height <- input$export_height}
    )
    observeEvent(input$avail_ids,
                 {values$avail_ids <- input$avail_ids}
    )
    
    # server managing gene lists --------------------------------------------------------
 
   
    
     
   
    
    
    available_orgdb <- rownames(installed.packages())[
      grep(pattern = "^org.*db$",rownames(installed.packages()))]
    

    
    
    # server ui update/observers --------------------------------------------------------
    output$color_by <- renderUI({
      # browser()
      if(is.null(values$dds_obj))
        return(NULL)
      poss_covars <- names(colData(values$dds_obj))
      selectInput('color_by', label = 'Group/color by: ',
                  choices = c(NULL, poss_covars), selected = "STIMULUS",multiple = TRUE)
    })
    
    # this trick speeds up the populating of the select(ize) input widgets,
    # see http://stackoverflow.com/questions/38438920/shiny-selectinput-very-slow-on-larger-data-15-000-entries-in-browser
    observe({
      updateSelectizeInput(
        session = session,
        inputId = 'avail_ids',
        choices = c(Choose = '', rownames(values$dds_obj)),
        server = TRUE)
    })
    
    
    observe({
      avVals = values$annotation_obj$gene_name[match(rownames(values$dds_obj), values$annotation_obj$gene_id)]
      updateSelectizeInput(
        session = session,
        # selected = oldSelected,
        inputId = 'avail_symbols',
        choices = c(Choose = '', avVals),
        server = TRUE)
    })
    
    output$available_genes <- renderUI({
      # browser()
      if(!is.null(values$annotation_obj)) {
        avVals = values$annotation_obj$gene_name[match(rownames(values$dds_obj), values$annotation_obj$gene_id)]
        # oldSelected = oldSelected[oldSelected %in% avVals]
        selectizeInput("avail_symbols", label = "Select the gene(s) of interest",
                       choices = avVals,
                       # selected = oldSelected,
                       multiple = TRUE)
      } else { # else use the rownames as identifiers
        selectizeInput("avail_ids", label = "Select the gene(s) of interest - ids",
                       choices = NULL, selected = NULL, multiple = TRUE)
      }
    })
    
    
    
    output$dds_design <- renderPrint({
      design(values$dds_obj)
    })
    
    output$res_names <- renderPrint({
      resultsNames(values$dds_obj)
    })
    
    # summary plots and outputs --------------------------------------------------------
    
    
    ### loading report template
    # update aceEditor module
    observe({
      # loading rmd report from disk
      inFile <- system.file("extdata", "irt.Rmd",package = "WEIN")
      
      isolate({
        if(!is.null(inFile) && !is.na(inFile)) {
          
          rmdfilecontent <- paste0(readLines(inFile),collapse="\n")
          
          shinyAce::updateAceEditor(session, "acereport_rmd", value = rmdfilecontent)
        }
      })
    })
    
    
    
    output$ui_iSEEexport <- renderUI({
      validate(
        need(((!is.null(values$dds_obj)) & (!is.null(values$res_obj))),
             message = "Please build and compute the dds and res object to export as 
             SummarizedExperiment for use in iSEE")
      )
      return(
        tagList(
          textInput(
            "se_export_name",label = "Choose a filename for the serialized .rds object",
            value = "se_WEIN_toiSEE.rds"),
          downloadButton(
            "button_iSEEexport",
            label = "Export as serialized SummarizedExperiment"
          )
        )
      )
    })
    
    output$button_iSEEexport <- downloadHandler(
      filename = function() {
        # paste0("se_WEIN_toiSEE_",gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time())))),".rds")
        input$se_export_name
      }, content = function(file) {
        se <- wrapup_for_iSEE(values$dds_obj, values$res_obj)
        saveRDS(se, file = file)
      }
    )
    
    # server state saving --------------------------------------------------------
    ### to environment
    observe({
      if(is.null(input$task_exit_and_save) || input$task_exit_and_save ==0 ) return()
      
      # quit R, unless you are running an interactive session
      if(interactive()) {
        # flush input and values to the environment in two distinct objects (to be reused later?)
        isolate({
          
          # WEIN_env <<- new.env(parent = emptyenv())
          cur_inputs <- reactiveValuesToList(input)
          cur_values <- reactiveValuesToList(values)
          tstamp <- gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time()))))
          
          # myvar <- "frfr"
          # assign("test", myvar, WEIN_env)
          
          # better practice rather than assigning to global env - notify users of this
          assign(paste0("WEIN_inputs_", tstamp),cur_inputs, envir = WEIN_env)
          assign(paste0("WEIN_values_", tstamp),cur_values, envir = WEIN_env)
          stopApp("WEIN closed, state successfully saved to global R environment.")
          
        })
      } else {
        stopApp("WEIN closed")
        q("no")
      }
    })
    
    ### to binary data
    saveState <- function(filename) {
      isolate({
        LiveInputs <- reactiveValuesToList(input)
        # values[names(LiveInputs)] <- LiveInputs
        r_data <- reactiveValuesToList(values)
        save(LiveInputs, r_data , file = filename)
      })
    }
    
    output$task_state_save <- downloadHandler(
      filename = function() {
        paste0("WEINState_",gsub(" ","_",gsub("-","",gsub(":","-",as.character(Sys.time())))),".RData")
      },
      content = function(file) {
        saveState(file)
      }
    )
    
    output$sessioninfo <- renderPrint({
      sessionInfo()
    })
    
    # server export plots and tables --------------------------------------------------------
    
    ## here, all export of plots and tables
    
    
    
    
    # Periodic memory cleanup
    observe({
      # Run garbage collection every 5 minutes
      invalidateLater(300000, session)  # 5 minutes in milliseconds
      gc()
    })
    
    # Cleanup temporary files when session ends
    session$onSessionEnded(function() {
      # Clean up any temporary files created during the session
      if (exists("userdir") && dir.exists(userdir)) {
        unlink(userdir, recursive = TRUE)
      }
    })
   }
 }
