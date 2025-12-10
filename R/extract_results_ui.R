extract_results_ui <-  function(){
  nsER <- NS("extract_results_manager")
  tabPanel(
    "Extract Results", icon = icon("table"),
    # see: http://stackoverflow.com/questions/21609436/r-shiny-conditionalpanel-output-value?noredirect=1&lq=1
    conditionalPanel(
      condition="!output.checkdds",
      headerPanel("Extract and inspect the DE results"),
      fluidRow(
        column(
          width = 8,
          shinyBS::bsCollapse(
            id = "help_extractresults",open = NULL,
            shinyBS::bsCollapsePanel(
              "Help",
              includeMarkdown(system.file("extdata", "help_results.md",package = "idealImmunoTP")))
          )
        )
      ),
      
      actionButton("tour_results", "Click me for a quick tour of the section", icon("info"),
                   style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"),
      br(),
      fluidRow(
        column(
          width = 6,
          uiOutput(nsER("choose_fac")),
          uiOutput(nsER("choose_fac2"))
        )
      ),
      fluidRow(
        column(
          width = 4,
          # factor as covariate
          wellPanel(
            width = 4, id = "factor_opts",
            uiOutput(nsER("fac1")),
            uiOutput(nsER("fac2"))
          )
          
        )
      ),
      
      ## general options for result function
      # alpha is set via FDR on the left side
      fluidRow(
        column(
          width = 4,
          wellPanel(id = "resu_opts",
                    selectInput(nsER("resu_indfil"),label = "Apply independent filtering automatically",
                                choices = c(TRUE,FALSE), selected = TRUE),
                   selectInput(nsER("resu_ihw"), "Use Independent Hypothesis Weighting (IHW) as a filtering function",
                                choices = c(TRUE, FALSE), selected = FALSE)
          )
        )
      ),
      #, evtl also the *filter* parameter of the function, i.e. baseMean if not specified
      fluidRow(
        column(
          width = 6,
          uiOutput(nsER("runresults")),
          uiOutput(nsER("store_result")),
          verbatimTextOutput(nsER("diyres_summary"))
        )
      ),
      
      DT::dataTableOutput(nsER("table_res")),
      downloadButton(nsER("downloadTblResu"),"Download", class = "btn btn-success"),
      fluidRow(
        column(
          width = 6,
          h3("Diagnostic plots"),
          plotOutput(nsER("pvals_hist")),
          div(align = "right", style = "margin-right:15px; margin-bottom:10px",
              downloadButton(nsER("download_plot_pvals_hist"), "Download Plot"),
              textInput(nsER("filename_plot_pvals_hist"),label = "Save as...",value = "plot_pvals_hist.pdf"))
        ),
        column(
          width = 6,
          plotOutput(nsER("pvals_hist_strat"))
        ),
         column(
          width = 6,
          plotOutput(nsER("logfc_hist")),
          div(align = "right", style = "margin-right:15px; margin-bottom:10px",
              downloadButton(nsER("download_plot_logfc_hist"), "Download Plot"),
              textInput(nsER("filename_plot_logfc_hist"),label = "Save as...",value = "plot_logfc_hist.pdf"))
        )
      )
    ),
    conditionalPanel(
      condition="output.checkdds",
      h2("You did not create the dds object yet. Please go the main tab and generate it")
    )
    
    
  )
}
