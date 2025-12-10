
panel_count_overview_ui <- function(){
  nsCO <- NS("count_overview_manager")
  tabPanel(
    "Counts Overview",
    icon = icon("eye"),
    conditionalPanel(
      condition="!output.checkdds",
      headerPanel("Get an overview on your data"),
      fluidRow(
        column(
          width = 8,
          shinyBS::bsCollapse(
            id = "help_countsoverview",open = NULL, 
            shinyBS::bsCollapsePanel(
              "Help",
              includeMarkdown(system.file("extdata", "help_overview.md",package = "idealImmunoTP")))
          )
        )
      ),
      
      actionButton("tour_countsoverview", "Click me for a quick tour of the section", 
                   icon("info"),
                   style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"),
      br(),
      selectInput(nsCO("countstable_unit"), label = "Data scale in the table",
                  choices = list("Counts (raw)" = "raw_counts",
                                 "Counts (normalized)" = "normalized_counts",
                                 # "Regularized logarithm transformed" = "rlog_counts",
                                 "Log10 (pseudocount of 1 added)" = "log10_counts"
                                 # ,
                                 # "TPM (Transcripts Per Million)" = "tpm_counts"
                  )),
      
      DT::dataTableOutput(nsCO("showcountmat")),
      downloadButton(nsCO("downloadData"),"Download", class = "btn btn-success"),
      hr(),
      fluidRow(
        column(
          width = 8,
          h3("Basic summary for the counts"),
          p("Number of uniquely aligned reads assigned to each sample"),
          # verbatimTextOutput("reads_summary"),
          wellPanel(
            fluidRow(
              column(
                width = 6,
                numericInput(nsCO("threshold_rowsums"),"Threshold on the row sums of the counts",value = 0, min = 0)),
              column(
                width = 6,
                numericInput(nsCO("threshold_rowmeans"),"Threshold on the row means of the normalized counts",value = 0, min = 0))
            )),
          p("According to the selected filtering criteria, this is an overview on the provided count data"),
          verbatimTextOutput(nsCO("detected_genes")),
          
          selectInput(nsCO("filter_crit"),label = "Choose the filtering criterium",
                      choices = c("row means", "row sums"), selected = "row means"),
          
          actionButton(nsCO("featfilt_dds"), "Filter the DDS object",class = "btn btn-primary")
        )
      ),
      
      h3("Sample to sample scatter plots"),
      selectInput(nsCO("corr_method"),"Correlation method",choices = list("pearson","spearman", "kendall")),
      checkboxInput(inputId = nsCO("corr_uselogs"),
                    label = "Use log2 values for plot axes and values",
                    value = TRUE),
      checkboxInput(inputId = nsCO("corr_usesubset"),
                    label = "Use a subset of max 1000 genes (quicker to plot)",
                    value = TRUE),
      p("Compute sample to sample correlations on the normalized counts - warning, it can take a while to plot all points (depending mostly on the number of samples you provided)."),
      actionButton(nsCO("compute_pairwisecorr"), "Run", class = "btn btn-primary"),
      uiOutput(nsCO("pairwise_plotUI")),
      br(),
      shinyjqui::jqui_resizable(uiOutput(nsCO("heatcorr_plotUI"))),
      br(),
      box(width = 12,
          title = "PCA/Size Factors", solidHeader = TRUE,
          fluidRow(
            column(
              width=3,
              uiOutput(nsCO("pca_Dim1"))
            ),
            column(
              width=3,
              uiOutput(nsCO("pca_Dim2"))
            )
          ),
          fluidRow(
            column(
              width = 6,
              shinyjqui::jqui_resizable(uiOutput(nsCO("pca_plotUI")))
            ),
            column(
              width = 6,
              shinyjqui::jqui_resizable(uiOutput(nsCO("pcaEV_plotUI")))
            )
          ),
          fluidRow(
            column(
              width = 6,
              shinyjqui::jqui_resizable(uiOutput(nsCO("pca_plotUI34")))
            ),
            column(
              width = 6,
              shinyjqui::jqui_resizable(uiOutput(nsCO("pca_plotUI56")))
            )
          ),
          fluidRow(
            
            column(
              width = 6,
              shinyjqui::jqui_resizable(uiOutput(nsCO("sizeFactors_plotUI")))
            )
          )
      ),
      br(),
      box(width = 12,
          title = "Highest expressing genes", solidHeader = TRUE,
          fluidRow(
            column(
              width = 12,
              uiOutput(nsCO("geneHeatmap_plotUI")) %>% shinyjqui::jqui_resizable()
            )
          ),
          fluidRow(
            column(
              width = 12,
              uiOutput(nsCO("geneHeatmap_genesUI")) %>% shinyjqui::jqui_resizable()
            )
            
          )
          
      ),
      conditionalPanel(
        condition="output.checkdds",
        h2("You did not create the dds object yet. Please go the main tab and generate it")
      )
    )
  )
}