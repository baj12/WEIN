data_setup_ui <- function(){
  nsOM    <- NS("ui_outputs_manager")
  nsSetup <- NS("ui_setup")
  nsTR <- NS("tour_manager")
  
  tabPanel(
    "Data Setup",icon = icon("upload"), # value="tab-ds",
    value = "tab-datasetup",
    headerPanel("Setup your data for the analysis"),
    fluidRow(
      column(
        width = 8,
        shinyBS::bsCollapse(
          id = nsTR("help_datasetup"),open = NULL, 
          shinyBS::bsCollapsePanel(
            "Help",
            includeMarkdown(system.file("extdata", "help_datasetup.md",package = "WEIN"))
          )
        )
      )
    ),
    
    actionButton(nsTR("tour_datasetup"), "Click me for a quick tour of the section", icon("info"),
                 style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),
    
    box(
      width = 12, 
      title = "Step 1", status = "danger", solidHeader = TRUE,
      h2("Upload your count matrix and the info on the experimental design"),
      
      fluidRow(
        column(
          width = 4, uiOutput(nsOM("upload_count_matrix")),uiOutput(nsOM("upload_metadata_ui"))
        ),
        column(
          width = 4,
          br(),
          actionButton("help_format",label = "",icon = icon("question-circle", verify_fa = FALSE),
                       style="color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"),
          shinyBS::bsTooltip(
            "help_format", 
            "How to provide your input data to WEIN",
            "bottom", options = list(container = "body")
          )
        )
      ),
      
      fluidRow(
        column(
          width = 6,
          box(width = NULL, title = "Count matrix preview",status = "primary",
              solidHeader = TRUE,collapsible = TRUE, collapsed = TRUE,
              fluidRow(
                column(
                  width = 12,
                  offset = 0.5,
                  DT::dataTableOutput(nsOM("dt_cm")))
              )
          )
        ),
        column(
          width = 6,
          box(width = NULL, title = "Experimental design preview",status = "primary",
              solidHeader = TRUE,collapsible = TRUE, collapsed = TRUE,
              fluidRow(
                column(
                  width = 12,
                  offset = 0.5,
                  DT::dataTableOutput(nsOM("dt_ed")))
              )
          )
        )
      ),
      fluidRow(
        column(
          width =6,
          box(width = NULL, title = "Non zero genes per sample",status = "primary",
              solidHeader = TRUE,collapsible = TRUE, collapsed = TRUE,
              fluidRow(
                column(
                  width = 12,
                  offset = 0.5,
                  shinyjqui::jqui_resizable(plotOutput(nsSetup("nonZeroCountsPlot"))))
              )
              
          )
        ),
        
        column(
          width =6,
          box(width = NULL, title = "Aligned sequences per sample",status = "primary",
              solidHeader = TRUE,collapsible = TRUE, collapsed = TRUE,
              fluidRow(
                column(
                  width = 12,
                  offset = 0.5,
                  shinyjqui::jqui_resizable(plotOutput(nsSetup("alignedSequencesPlot"))))
              )
              
          )
        )
      )
    ),
    uiOutput(nsOM("ui_step2")),
    fluidRow(
      column(
        width = 6, uiOutput(nsOM("ui_stepanno"))
      ),
      column(
        width = 6,uiOutput(nsOM("ui_stepoutlier"))
      )
    ),
    uiOutput(nsOM("ui_step3"))
  ) # end of Data Setup panel
}